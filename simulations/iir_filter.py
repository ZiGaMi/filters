# ===============================================================================
# @file:    iir_filter.py
# @note:    This script is evaluation of IIR filter algorithm
# @author:  Ziga Miklosic
# @date:    06.01.2021
# @brief:   Evaluation of IIR filter design. 
# ===============================================================================

# ===============================================================================
#       IMPORTS  
# ===============================================================================
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import freqz, butter, cheby1

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

## IIR Filter
class IIR:

    def __init__(self, a, b, order=1):

        # Store tap number and coefficient
        self.order = order
        self.a = a
        self.b = b

        # Create circular buffer
        self.x = CircBuffer(order)
        self.y = CircBuffer(order)


    def update(self, x):
        
        # Fill input
        self.x.set( x )

        """
        # Get tail index of circ buffer
        tail = self.x.get_tail()

        # Offset by one
        if tail >= ( self.order - 1 ):
            tail = 0
        else:
            tail = tail + 1

        # Convolve
        y = 0
        for j in range(self.order + 1):
            
            y += ( self.coef[j] * self.x.get( tail ))

            if tail >= ( self.order - 1 ):
                tail = 0
            else:
                tail = tail - 1
        """


        # Fill output
        self.y.set(y)

        return y


# ===============================================================================
#       MAIN ENTRY
# ===============================================================================
if __name__ == "__main__":

    # Time array
    _time, _dt = np.linspace( 0.0, TIME_WINDOW, num=SAMPLE_NUM, retstep=True )

    # Natural frequency
    _w = 1.0

    # Damping factor
    _z = .707

    # Cutoff frequency
    _fc = 100.0

    # Sample time
    _ts = 2 * np.tan( _fc / SAMPLE_FREQ * np.pi ) 

    # Calculate frequency characteristics 
    a2 = ( _ts**2 *_w**2 / 4 ) - ( _ts/2 * _z * _w ) + 1
    a1 = ( _ts**2 * _w**2 / 2 ) - 2
    a0 = ( _ts**2 *_w**2 / 4 ) + ( _ts/2 * _z * _w ) + 1
    
    b2 = 1
    b1 = -2
    b0 = 1
    
    a = [ a0, a1, a2 ] 
    b = [ b0, b1, b2 ] 


   
    _omega_b = 2 * _fc / SAMPLE_FREQ

    #b_b, a_b = butter( N=1, Wn=_omega_b, btype="lowpass", analog=False, fs=SAMPLE_FREQ )
    b_b, a_b = butter( N=1, Wn=_fc, btype="highpass", analog=False, fs=SAMPLE_FREQ )

    print("a: %s" % a)
    print("b: %s" % b)
    print("a_b: %s" % a_b)
    print("b_b: %s" % b_b)

    w, h = freqz( b, a, 4096 )
    w_b, h_b = freqz( b_b, a_b, 4096 )
    
    
    """
    # Filter object
    _filter_IIR = IIR( 2, a, b ) 

    # Filter input/output
    _x = [ 0 ] * SAMPLE_NUM
    _x_d = [0]

    _y_d_iir = [0]

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

            # IIR 
            _y_d_iir.append( _filter_IIR.update( _x[n] ))

        else:
            _downsamp_cnt += 1
    
    """
    
    # Plot results
    fig, ax = plt.subplots(2, 1)
    fig.suptitle("IIR Filter\nInput signal freq: " + str(INPUT_SIGNAL_FREQ) + "Hz, Sample freq (fs): " +str(SAMPLE_FREQ) + "Hz", fontsize=20)
    
    """
    ax[0].plot( _time, _x,                  "b", label="Input-generated" )
    ax[0].plot( _d_time, _downsamp_samp,    "r.", label="Sample points")
    ax[0].plot( _d_time, _y_d_iir,          "g" )
    ax[0].grid()
    ax[0].title.set_text("Time domain")
    ax[0].set_ylabel("Amplitude")
    ax[0].legend(loc="upper right")


    """

    
    # Convert to Hz unit
    w = ( w / np.pi * SAMPLE_FREQ / 2)
    w_b = ( w_b / np.pi * SAMPLE_FREQ / 2)

    ax[1].plot(w, 20 * np.log10(abs(h)), 'b')
    ax[1].plot(w_b, 20 * np.log10(abs(h_b)), 'g')
    ax[1].set_ylabel('Amplitude [dB]', color='b')
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].grid()
    ax[1].title.set_text("Amplitude and phase characteristics")
    #ax[1].legend(loc="upper right")

    """
    ax_11 = ax[1].twinx()
    angles = np.unwrap( np.angle(h) )
    angles_b = np.unwrap( np.angle(h_b) )
    ax_11.plot(w, (angles*180/np.pi), 'b')
    ax_11.plot(w_b, (angles_b*180/np.pi), 'g')
    ax_11.set_ylabel('Angle [degrees]', color='g')
    ax_11.axis('tight')
    #ax_11.legend(loc="lower right")
    """
    

    plt.show()
    


# ===============================================================================
#       END OF FILE
# ===============================================================================
