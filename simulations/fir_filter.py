# ===============================================================================
# @file:    fir_filter.py
# @note:    This script is evaluation of FIR filter algorithm
# @author:  Ziga Miklosic
# @date:    02.01.2021
# @brief:   Evaluation of FIR filter design. FIR1 filter coefficients are 
#           calculated using firwin function while FIR2 filter coefficients
#           are defined via webpage T-filter.
# ===============================================================================

# ===============================================================================
#       IMPORTS  
# ===============================================================================
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import firwin, freqz, kaiserord

# ===============================================================================
#       CONSTANTS
# ===============================================================================

## ****** USER CONFIGURATIONS ******

## Sample frequency
#   Sample frequency of real system   
#
# Unit: Hz
SAMPLE_FREQ = 1000.0

# Ideal sample frequency
#   As a reference to sample rate constrained embedded system
#
# Unit: Hz
IDEAL_SAMPLE_FREQ = 20000.0

## Time window
#
# Unit: second
TIME_WINDOW = 2.5

## Number of samples in time window
SAMPLE_NUM = int(( IDEAL_SAMPLE_FREQ * TIME_WINDOW ) + 1.0 )

## Select input filter signal type
INPUT_SINE = 0
INPUT_RECT = 1

## Mux input signal
INPUT_SIGNAL_SELECTION = INPUT_SINE

## Input signal frequency
#
# Unit: Hz
INPUT_SIGNAL_FREQ = 50

## Number of FIR taps
FIR_TAP_NUM     = 8
FIR_TAP_NUM_2   = 19

## FIR cutoff frequency
#
#   Unit: Hz
FIR_CUT_OFF     = 10.0

# Coefficient calculated by T-filter web page
# Link: http://t-filter.engineerjs.com/
FIR_2_COEF = [
    0.00983679251710613,
    0.016906009403821313,
    0.029012511729722044,
    0.044185823706193596,
    0.06146104679182853,
    0.07937592793186861,
    0.09611093303268682,
    0.10977246099391554,
    0.11872311326692045,
    0.12184110110169752,
    0.11872311326692045,
    0.10977246099391554,
    0.09611093303268682,
    0.07937592793186861,
    0.06146104679182853,
    0.044185823706193596,
    0.029012511729722044,
    0.016906009403821313,
    0.00983679251710613
]

## ****** END OF USER CONFIGURATIONS ******

# ===============================================================================
#       FUNCTIONS
# ===============================================================================

# ===============================================================================
# @brief: Input signal mux
#
# @param[in]:    sel     - Multiplexor selector  
# @param[in]:    in_1    - Input 1 
# @param[in]:    in_2    - Input 2 
# @return:       Either in_1 or in_2
# ===============================================================================
def input_signal_mux(sel, in_1, in_2):
    if ( INPUT_SINE == sel ):
        return in_1
    elif ( INPUT_RECT == sel ):
        return in_2
    else:
        pass

# ===============================================================================
# @brief: Generate sine as injected signal
#
# @param[in]:    time    - Linear time  
# @param[in]:    amp     - Amplitude of sine
# @param[in]:    off     - DC offset of sine
# @param[in]:    phase   - Phase of sine
# @return:       Generated signal
# ===============================================================================
def generate_sine(time, freq, amp, off, phase):
    return (( amp * np.sin((2*np.pi*freq*time) + phase )) + off )


# ===============================================================================
# @brief: Generate rectangle signal
#
# @param[in]:    time    - Linear time  
# @param[in]:    amp     - Amplitude of rectange
# @param[in]:    off     - DC offset of rectangle
# @param[in]:    phase   - Phase of rectangle
# @return:       Generated signal
# ===============================================================================
def generate_rect(time, freq, amp, off, phase):
    _carier = generate_sine(time, freq, 1.0, 0.0, phase)
    _sig = 0

    if ( _carier > 0 ):
        _sig = amp + off
    else:
        _sig = off

    return _sig 


# ===============================================================================
#       CLASSES
# ===============================================================================

## Circular buffer
class CircBuffer:

    def __init__(self, size):
        self.buf = [0.0] * size
        self.idx = 0
        self.size = size
    
    def __manage_index(self):
        if self.idx >= (self.size - 1):
            self.idx = 0
        else:
            self.idx = self.idx + 1

    def set(self, val):
        if self.idx < self.size:
            self.buf[self.idx] = val
            self.__manage_index()
        else:
            raise AssertionError

    def get(self, idx):
        if idx < self.size:
            return self.buf[idx]
        else:
            raise AssertionError
    
    def get_tail(self):
        return self.idx

    def get_whole_buffer(self):
        return self.buf

## FIR Filter
class FIR:

    def __init__(self, tap, coef):

        # Store tap number and coefficient
        self.tap = tap
        self.coef = coef

        # Create circular buffer
        self.buf = CircBuffer(size=tap)


    def update(self, x):
        
        # Fill buffer
        self.buf.set( x )

        # Get tail index of circ buffer
        tail = self.buf.get_tail()

        # Offset by one
        if tail >= ( self.tap - 1 ):
            tail = 0
        else:
            tail = tail + 1

        # Convolve
        y = 0
        for j in range(self.tap):
            #y += ( self.coef[j] * self.buf.get( self.tap - j ))
            
            y += ( self.coef[j] * self.buf.get( tail ))

            if tail >= ( self.tap - 1 ):
                tail = 0
            else:
                tail = tail + 1

        return y


# ===============================================================================
#       MAIN ENTRY
# ===============================================================================
if __name__ == "__main__":

    # Time array
    _time, _dt = np.linspace( 0.0, TIME_WINDOW, num=SAMPLE_NUM, retstep=True )

    # Calculate FIR coefficient
    _fir_coef = firwin( numtaps=FIR_TAP_NUM, cutoff=FIR_CUT_OFF, window="hamming", fs=SAMPLE_FREQ )

    # Calculate frequency characteristics 
    w, h = freqz( _fir_coef )
    w2, h2 = freqz( FIR_2_COEF )

    # Filter object
    _filter_FIR = FIR( FIR_TAP_NUM, _fir_coef )
    _filter_FIR_2 = FIR( FIR_TAP_NUM_2, FIR_2_COEF )

    # Filter input/output
    _x = [ 0 ] * SAMPLE_NUM
    _x_d = [0]

    _y_d_fir = [0]
    _y_d_fir_2 = [0]

    # Generate inputs
    _sin_x = []
    _rect_x = []

    # Down sample
    _downsamp_cnt = 0
    _downsamp_samp = [0]
    _d_time = [0]
    
    for n in range(SAMPLE_NUM):
        _sin_x.append( generate_sine( _time[n], INPUT_SIGNAL_FREQ, 1.0, 0.0, 0.0 ))  
        _rect_x.append( generate_rect( _time[n], INPUT_SIGNAL_FREQ, 1.0, 0.0, 0.0 ))  
 
    # Apply filter
    for n in range(SAMPLE_NUM):
        
        # Mux input signals
        _x[n] = input_signal_mux( INPUT_SIGNAL_SELECTION, _sin_x[n], _rect_x[n] )

        # Down sample to SAMPLE_FREQ
        if _downsamp_cnt >= (( 1 / ( _dt * SAMPLE_FREQ )) - 1 ):
            _downsamp_cnt = 0

            # Utils
            _downsamp_samp.append(0)
            _d_time.append( _time[n])
            _x_d.append( _x[n] )

            # FIR 
            _y_d_fir.append( _filter_FIR.update( _x[n] ) )
            _y_d_fir_2.append( _filter_FIR_2.update( _x[n] ) )
        else:
            _downsamp_cnt += 1
    
    # Plot results
    fig, ax = plt.subplots(2, 1)
    fig.suptitle("FIR Filter\nInput signal freq: " + str(INPUT_SIGNAL_FREQ) + "Hz, Sample freq (fs): " +str(SAMPLE_FREQ) + "Hz", fontsize=20)
    ax[0].plot( _time, _x,                  "b", label="Input-generated" )
    ax[0].plot( _d_time, _downsamp_samp,    "r.", label="Sample points")
    ax[0].plot( _d_time, _y_d_fir,          "g", label=str(FIR_CUT_OFF) + "Hz/" + str(FIR_TAP_NUM))
    ax[0].plot( _d_time, _y_d_fir_2,        "y", label="xHz/" + str(FIR_TAP_NUM_2))
    ax[0].grid()
    ax[0].title.set_text("Time domain")
    ax[0].set_ylabel("Amplitude")
    ax[0].legend(loc="upper right")

    # Convert to Hz unit
    w = ( w / np.pi * SAMPLE_FREQ / 2)
    w2 = ( w2 / np.pi * SAMPLE_FREQ / 2)

    ax[1].plot(w, 20 * np.log10(abs(h)), 'b', label=str(FIR_CUT_OFF) + "Hz/" + str(FIR_TAP_NUM))
    ax[1].plot(w2, 20 * np.log10(abs(h2)), 'y', label="xHz/" + str(FIR_TAP_NUM_2))
    ax[1].set_ylabel('Amplitude [dB]', color='b')
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].grid()
    ax[1].title.set_text("Amplitude and phase characteristics")
    ax[1].legend(loc="upper right")

    ax_11 = ax[1].twinx()
    angles = np.unwrap( np.angle(h) )
    angles_2 = np.unwrap( np.angle(h2) )
    ax_11.plot(w, (angles*180/np.pi), 'g', label="phase 1")
    ax_11.plot(w2, (angles_2*180/np.pi), 'r', label="phase 2")
    ax_11.set_ylabel('Angle [degrees]', color='g')
    ax_11.axis('tight')
    ax_11.legend(loc="lower right")

    plt.show()
    


# ===============================================================================
#       END OF FILE
# ===============================================================================
