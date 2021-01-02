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

# ===============================================================================
#       CONSTANTS
# ===============================================================================

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
TIME_WINDOW = 0.5

## Number of samples in time window
SAMPLE_NUM = int(( IDEAL_SAMPLE_FREQ * TIME_WINDOW ) + 1.0 )

## Select input filter signal type
INPUT_SINE = 0
INPUT_RECT = 1

## Mux input signal
INPUT_SIGNAL_SELECTION = INPUT_RECT

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


class CircularBuffer:

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
        if self.idx < self.size:
            return self.buf[self.idx]
        else:
            raise AssertionError

class FIR:

    def __init__(self, tap, coef):

        # Store tap number and coefficient
        self.tap = tap
        self.coef = coef

        # Create circular buffer
        self.buf = CircularBuffer(size=tap)


    def update(self, x):
        
        # Fill buffer
        self.buf.set( x )

        # Convolve
        y = 0
        for j in range(self.tap):
            y += ( self.coef[j] * self.buf.get( self.tap - j ))

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

    # Filter object
    _filter_FIR     = FIR( FIR_TAP_NUM, FIR_COEFFICIENT )

    # Filter input/output
    _x = [ 0 ] * SAMPLE_NUM
    _x_d = [0]

    _y_fir = []
    _time_fir = []

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

            # LPF
            _y_d_lpf_1.append( _filter_D_LPF_1.update( _x_d[-1] ) )
            _y_d_lpf_2.append( _filter_D_LPF_2.update( _x_d[-1] ) )
            _y_d_lpf_3.append( _filter_D_LPF_3.update( _x_d[-1] ) )

            # FIR 
            _y_fir.append( _filter_FIR.update( _x[n] ) )
            _time_fir.append( _time[n] )
        else:
            _downsamp_cnt += 1
    
    # Plot results
    """
    fig, ax = plt.subplots(2, 1)
    fig.suptitle("Input signal freq: " + str(INPUT_SIGNAL_FREQ) + "Hz", fontsize=20)
    ax[0,0].plot( _time, _x, "b", label="Input-generated" )
    ax[0,0].plot( _time, _y_lpf_1, "g", label="Ideal filter")
    ax[0,0].plot( _d_time, _downsamp_samp, "r.", label="Sample points")
    ax[0,0].plot( _time, _y_lpf_1, "g", label="RC1: " + str(LPF_FC_2) + "Hz/" + str(LPF_ORDER_2))
    ax[0,0].plot( _time, _y_lpf_2, "r", label="RC2: " + str(LPF_FC_2) + "Hz/" + str(LPF_ORDER_2))
    ax[0,0].plot( _time, _y_lpf_3, "y", label="RC3: " + str(LPF_FC_3) + "Hz/" + str(LPF_ORDER_3))
    ax[0,0].grid()
    ax[0,0].title.set_text("RC Low Pass Filter - Ideal (fs=" + str( 1 / _dt ) + "Hz)")
    ax[0,0].set_ylabel("Amplitude")
    ax[0,0].legend(loc="upper right")


    ax[0,1].plot( _d_time, _x_d, "b.-", label="Input-sampled")
    ax[0,1].plot( _d_time, _y_d_lpf_1, "g.-", label="RC1: " + str(LPF_FC_1) + "Hz/" + str(LPF_ORDER_1))
    ax[0,1].plot( _d_time, _y_d_lpf_2, "r.-", label="RC2: " + str(LPF_FC_2) + "Hz/" + str(LPF_ORDER_2))
    ax[0,1].plot( _d_time, _y_d_lpf_3, "y.-", label="RC3: " + str(LPF_FC_3) + "Hz/" + str(LPF_ORDER_3))
    ax[0,1].grid()
    ax[0,1].title.set_text("RC Low Pass Filter - Real (fs=" + str(SAMPLE_FREQ) + "Hz)")
    ax[0,1].legend(loc="upper right")
    """
    
    plt.subplots_adjust(left=0.05, right=0.98, bottom=0.05, wspace=0.08)

    plt.show()
    


# ===============================================================================
#       END OF FILE
# ===============================================================================
