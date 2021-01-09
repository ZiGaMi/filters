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
from scipy.signal import freqz, butter, cheby1, lfilter, filtfilt

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
TIME_WINDOW = 1

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
INPUT_SIGNAL_FREQ = 10.0

## Cutoff freqeuncy of filter
#
# Unit: Hz
HPF_FC_1 = 10.0
HPF_FC_2 = 10.0
HPF_FC_BUTTER = 10.0
HPF_FC_CHEBY = 10.0

## Damping factor (zeta) for 2nd order filter
#
#   z = 0       -> underdamped
#   z = 0.7071  -> criticaly damped, sweep spot
#   z = 1       -> overdamped
HPF_Z_1 = .25
HPF_Z_2 = 1.75


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
# @brief:   calculate 2nd order high pass filter based on following 
#           transfer function:
#           
#               h(s) = s^2 / ( s^2 + z*s + 1 )
#
# @param[in]:    fc     - Corner frequenc
# @param[in]:    z      - Damping factor
# @param[in]:    fs     - Sample frequency
# @return:       b,a    - Array of b,a IIR coefficients
# ===============================================================================
def calculate_2nd_order_HPF_coeff(fc, z, fs):
    
    _ts = 2 * np.tan( fc / SAMPLE_FREQ * np.pi ) 
    
    # Calculate coefficient
    # NOTE: This for of coefficient is result of bi-linear transform
    a2 = ( _ts**2 / 4 ) - ( _ts/2 * z ) + 1
    a1 = ( _ts**2 / 2 ) - 2
    a0 = ( _ts**2 / 4 ) + ( _ts/2 * z ) + 1
    b2 = 1
    b1 = -2
    b0 = 1

    # Fill array
    a = [ a0, a1, a2 ] 
    b = [ b0, b1, b2 ] 

    return b, a


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

    """
        Returns array of time sorted samples in buffer
        [ n, n-1, n-2, ... n - size - 1 ]
    """
    def get_time_ordered_samples(self):

        _ordered = [0.0] * self.size
        _start_idx = 0

        _start_idx = self.idx - 1
        if _start_idx < 0:
            _start_idx += self.size

        # Sort samples per time
        for n in range( self.size ):
            _last_idx = _start_idx - n
            if _last_idx < 0:
                _last_idx += self.size 
            _ordered[n] = self.buf[_last_idx]

        return _ordered

    def get_tail(self):
        return self.idx

    def get_size(self):
        return self.size

    def get_whole_buffer(self):
        return self.buf

## IIR Filter
class IIR:

    def __init__(self, a, b, order):

        # Store tap number and coefficient
        self.order = order
        self.a = a
        self.b = b

        # Create circular buffer
        self.x = CircBuffer(order+1)
        self.y = CircBuffer(order+1)


    def update(self, x):
        
        # Fill input
        self.x.set(float( x ))

        # Get input/outputs history
        _x = self.x.get_time_ordered_samples()
        _y = self.y.get_time_ordered_samples()

        # Calculate new value
        y = 0.0
        for j in range(self.order+1):

            y = y + float(float(self.b[j]) * _x[j])

            if j > 0:
                y = y - float(float(self.a[j]) * _y[j-1])

        y = float( y * ( 1 / float(self.a[0] )))

        # Fill output
        self.y.set(y)

        return y


# ===============================================================================
#       MAIN ENTRY
# ===============================================================================
if __name__ == "__main__":

    # Time array
    _time, _dt = np.linspace( 0.0, TIME_WINDOW, num=SAMPLE_NUM, retstep=True )

    # Calculate coefficient for 2nd order HPF
    b, a        = calculate_2nd_order_HPF_coeff( HPF_FC_1, HPF_Z_1, SAMPLE_FREQ )
    b_2, a_2    = calculate_2nd_order_HPF_coeff( HPF_FC_2, HPF_Z_2, SAMPLE_FREQ )

    # Calculate coefficient for 2nd order Butterworth & Chebyshev filter
    b_b, a_b = butter( N=2, Wn=HPF_FC_BUTTER, btype="highpass", analog=False, fs=SAMPLE_FREQ )
    b_c, a_c = cheby1( N=2, Wn=HPF_FC_CHEBY, btype="highpass", analog=False, fs=SAMPLE_FREQ, rp=0.99 )

    # Get frequency characteristics
    w, h = freqz( b, a, 4096 )
    w_2, h_2 = freqz( b_2, a_2, 4096 )
    w_b, h_b = freqz( b_b, a_b, 4096 )
    w_c, h_c = freqz( b_c, a_c, 4096 )
    
    # Filter object
    _filter_IIR     = IIR( a, b, order=2 ) 
    _filter_IIR_2   = IIR( a_2, b_2, order=2 ) 

    # Filter input/output
    _x = [ 0 ] * SAMPLE_NUM
    _x_d = [0]

    _y_d_iir = [0]
    _y_d_iir_2 = [0]

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
            _y_d_iir_2.append( _filter_IIR_2.update( _x[n] ))

        else:
            _downsamp_cnt += 1
    
    # Apply filter on signal
    #__y = lfilter(b, a, _x_d)
    
    # Plot results
    fig, ax = plt.subplots(2, 1)
    fig.suptitle("Highpass 2nd order IIR Filter Design\nInput signal freq: " + str(INPUT_SIGNAL_FREQ) + "Hz, Sample freq (fs): " +str(SAMPLE_FREQ) + "Hz", fontsize=20)
    
    ax[0].plot( _time, _x,                  "b",    label="Input-generated" )
    ax[0].plot( _d_time, _downsamp_samp,    "r.",   label="Sample points")
    ax[0].plot( _d_time, _y_d_iir,          ".-g",    label=str(HPF_FC_1) + "Hz/" + str(HPF_Z_1))
    ax[0].plot( _d_time, _y_d_iir_2,        ".-y",    label=str(HPF_FC_2) + "Hz/" + str(HPF_Z_2))
    ax[0].grid()
    ax[0].title.set_text("Time domain")
    ax[0].set_ylabel("Amplitude")
    ax[0].legend(loc="upper right")

    
    # Convert to Hz unit
    w = ( w / np.pi * SAMPLE_FREQ / 2)
    w_2 = ( w_2 / np.pi * SAMPLE_FREQ / 2)

    w_b = ( w_b / np.pi * SAMPLE_FREQ / 2)
    w_c = ( w_c / np.pi * SAMPLE_FREQ / 2)

    ax[1].plot(w,   20 * np.log10(abs(h)),   'g', label=str(HPF_FC_1) + "Hz/" + str(HPF_Z_1))
    ax[1].plot(w_2, 20 * np.log10(abs(h_2)), 'y', label=str(HPF_FC_2) + "Hz/" + str(HPF_Z_2))
    
    ax[1].plot(w_b, 20 * np.log10(abs(h_b)), 'r', label="butterworth" )
    ax[1].plot(w_c, 20 * np.log10(abs(h_c)), 'b', label="chebysev" )
    
    ax[1].set_ylabel('Amplitude [dB]')
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].grid()
    ax[1].title.set_text("Amplitude and phase characteristics")
    ax[1].set_xscale("log")
    ax[1].legend(loc="upper right")
    
    """
    ax_11 = ax[1].twinx()
    angles = np.unwrap( np.angle(h) )
    angles_b = np.unwrap( np.angle(h_b) )
    ax_11.plot(w, (angles*180/np.pi), 'b')
    ax_11.plot(w_b, (angles_b*180/np.pi), 'r')
    ax_11.set_ylabel('Angle [degrees]', color='g')
    ax_11.axis('tight')
    """

    plt.show()
    


# ===============================================================================
#       END OF FILE
# ===============================================================================
