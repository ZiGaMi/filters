# ===============================================================================
# @file:    fir_filter.py
# @note:    This script is evaluation of FIR filter algorithm
# @author:  Ziga Miklosic
# @date:    02.01.2021
# @brief:   Evaluation of FIR filter design
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

## Sample frequency
#   Sample frequency of real system   
#
# Unit: Hz
SAMPLE_FREQ = 2000.0

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
INPUT_SIGNAL_FREQ = 100



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

FIR_TAP_NUM = 16
FIR_COEFFICIENT = [
993.7115243977110590E-6,
 0.011898750525239039,
 0.026093294471711985,
 0.048239818572596682,
 0.074979302516632151,
 0.102206493068131094,
 0.124584980431458908,
 0.137240329406428413,
 0.137240329406428413,
 0.124584980431458908,
 0.102206493068131094,
 0.074979302516632151,
 0.048239818572596682,
 0.026093294471711985,
 0.011898750525239039,
 993.7115243977110590E-6,
]

if __name__ == "__main__":

    # Time array
    _time, _dt = np.linspace( 0.0, TIME_WINDOW, num=SAMPLE_NUM, retstep=True )

    FIR_TAPS = 8
    #_fir_coef = firwin( numtaps=FIR_TAPS, cutoff=10.0, window="hamming", fs=SAMPLE_FREQ )
    _fir_coef = firwin( numtaps=FIR_TAPS, cutoff=[0.01, 100], window="hamming", fs=SAMPLE_FREQ, pass_zero=False )
    _fir_coef = FIR_COEFFICIENT

    # Filter object
    _filter_FIR     = FIR( FIR_TAP_NUM, FIR_COEFFICIENT )
    _filter_FIR_2   = FIR( FIR_TAPS, _fir_coef )

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
        #_sin_x.append( generate_sine( _time[n], INPUT_SIGNAL_FREQ, 1.0, 0.0, 0.0 ) + generate_sine( _time[n], 1500.0, 0.05, 0.0, 0.0 )  )  
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
    fig.suptitle("Input signal freq: " + str(INPUT_SIGNAL_FREQ) + "Hz", fontsize=20)
    ax[0].plot( _time, _x,                  "b", label="Input-generated" )
    ax[0].plot( _d_time, _downsamp_samp,    "r.", label="Sample points")
    #ax[0].plot( _d_time, _y_d_fir,          "g", label="FIR")
    ax[0].plot( _d_time, _y_d_fir_2,        "y", label="FIR2")
    ax[0].grid()
    #ax[0].title.set_text("RC Low Pass Filter - Ideal (fs=" + str( 1 / _dt ) + "Hz)")
    ax[0].set_ylabel("Amplitude")
    ax[0].legend(loc="upper right")


    w, h = freqz( _fir_coef )

    # Convert to Hz unit
    w = ( w / np.pi * SAMPLE_FREQ / 2)

    ax[1].plot(w, 20 * np.log10(abs(h)), 'b')
    ax[1].set_ylabel('Amplitude [dB]', color='b')
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].grid()

    ax_11 = ax[1].twinx()
    angles = np.unwrap( np.angle(h) )
    ax_11.plot(w, (angles*180/np.pi), 'g')
    ax_11.set_ylabel('Angle [degrees]', color='g')
    ax_11.axis('tight')


    """
    # Nyquist rate.
    nyq_rate = SAMPLE_FREQ / 2

    # Width of the roll-off region.
    width = 500 / nyq_rate

    # Attenuation in the stop band.
    ripple_db = 12.0

    num_of_taps, beta = kaiserord(ripple_db, width)
    if num_of_taps % 2 == 0:
        num_of_taps = num_of_taps + 1

    # Cut-off frequency.
    cutoff_hz = 5000.0

    # Estimate the filter coefficients.
    taps = firwin(num_of_taps, cutoff_hz/nyq_rate, window=('kaiser', beta), pass_zero=False)

    w, h = freqz(taps, worN=4000)


    fig, ax = plt.subplots(1, 1)
    ax.plot((w/np.pi)*nyq_rate, 20*np.log10(np.abs(h)), linewidth=2)

    ax.axvline(cutoff_hz + width*nyq_rate, linestyle='--', linewidth=1, color='g')
    ax.axvline(cutoff_hz - width*nyq_rate, linestyle='--', linewidth=1, color='g')
    ax.axhline(-ripple_db, linestyle='--', linewidth=1, color='c')
    delta = 10**(-ripple_db/20)
    ax.axhline(20*np.log10(1 + delta), linestyle='--', linewidth=1, color='r')
    ax.axhline(20*np.log10(1 - delta), linestyle='--', linewidth=1, color='r')
    ax.grid()
    """


    #plt.subplots_adjust(left=0.05, right=0.98, bottom=0.05, wspace=0.08)

    plt.show()
    


# ===============================================================================
#       END OF FILE
# ===============================================================================
