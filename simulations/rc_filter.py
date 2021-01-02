# ===============================================================================
# @file:    rc_filter.py
# @note:    This script is evaluation of RC filter algorithm
# @author:  Ziga Miklosic
# @date:    31.12.2020
# @brief:   Evaluation of RC/CR filter design. Implemented are HPF and LPF
#           with configurable fc and order. 
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
#
# Unit: Hz
SAMPLE_FREQ = 100

## Time window
#
# Unit: second
TIME_WINDOW = 0.5

## Number of samples in time window
SAMPLE_NUM = 10000

## Select input filter signal type
INPUT_SINE = 0
INPUT_RECT = 1

## Mux input signal
INPUT_SIGNAL_SELECTION = INPUT_SINE

## Input signal frequency
#
# Unit: Hz
INPUT_SIGNAL_FREQ = 10

## LPF fc
#
# Unit: Hz
LPF_FC_1 = 10.0
LPF_FC_2 = 100.0
LPF_FC_3 = 100.0

## LPF order
LPF_ORDER_1 = 1
LPF_ORDER_2 = 10
LPF_ORDER_3 = 5

## HPF fc
#
# Unit: Hz
HPF_FC_1 = 10.0
HPF_FC_2 = 20.0
HPF_FC_3 = 30.0

## HPF order
HPF_ORDER_1 = 1
HPF_ORDER_2 = 10
HPF_ORDER_3 = 5


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

## RC Low Pass Filter
class RC_LPF:

    def __init__(self, fc, dt, order=1, init_val=0):

        # Calculate filter coefficient
        time_const = 1/(2*np.pi*fc)
        self.k = dt / ( dt + time_const )

        # Store order
        self.order = order

        # Init output of filter
        self.y = [init_val] * order

    def update(self, x):
        for n in range( self.order ):
            if n == 0:
                self.y[0] = ( self.y[0] + self.k * ( x - self.y[0] ))
            else:
                self.y[n] = ( self.y[n] + self.k * ( self.y[n-1] - self.y[n] ))

        return self.y[ self.order - 1 ]
        


## CR High Pass Filter
class CR_HPF:

    def __init__(self, fc, dt, order=1):

        # Calculate filter coefficient
        time_const = 1/(2*np.pi*fc)
        self.k = time_const / ( dt + time_const )

        # Store order
        self.order = order

        # Init inputs & outputs
        self.y = [0] * order
        self.x = [0] * order

    def update(self, x):

        for n in range( self.order ):
            if n == 0:
                self.y[n] = ( self.k * self.y[n] + self.k * ( x - self.x[0] ))
                self.x[0] = x
            else:
                self.y[n] = ( self.k * self.y[n] + self.k * ( self.y[n-1] - self.x[n] ))
                self.x[n] = self.y[n-1]

        return self.y[ self.order - 1 ]



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
    _filter_LPF_1   = RC_LPF( fc=LPF_FC_1, dt=_dt, order=LPF_ORDER_1, init_val=0)
    #_filter_LPF_2   = RC_LPF( fc=LPF_FC_2, dt=_dt, order=LPF_ORDER_2, init_val=1.0)
    #_filter_LPF_3   = RC_LPF( fc=LPF_FC_3, dt=_dt, order=LPF_ORDER_3, init_val=1.5)
    
    _filter_HPF_1   = CR_HPF( fc=HPF_FC_1, dt=_dt, order=HPF_ORDER_1)
    _filter_HPF_2   = CR_HPF( fc=HPF_FC_2, dt=_dt, order=HPF_ORDER_2)
    _filter_HPF_3   = CR_HPF( fc=HPF_FC_3, dt=_dt, order=HPF_ORDER_3)

    _filter_FIR     = FIR( FIR_TAP_NUM, FIR_COEFFICIENT )

    # Filter input/output
    _x = [ 0 ] * SAMPLE_NUM

    _y_lpf_1 = []
    _y_lpf_2 = []
    _y_lpf_3 = []

    _y_hpf_1 = []
    _y_hpf_2 = []
    _y_hpf_3 = []

    _y_fir = []
    _time_fir = []

    # Generate inputs
    _sin_x = []
    _rect_x = []

    # Down sample conunter
    _downsamp_cnt = 0
    
    for n in range(SAMPLE_NUM):
        _sin_x.append( generate_sine( _time[n], INPUT_SIGNAL_FREQ, 1.0, 0.0, 0.0 ))  
        _rect_x.append( generate_rect( _time[n], INPUT_SIGNAL_FREQ, 1.0, 0.0, 0.0 ))  
 
    # Apply filter
    for n in range(SAMPLE_NUM):
        
        # Mux input signals
        _x[n] = input_signal_mux( INPUT_SIGNAL_SELECTION, _sin_x[n], _rect_x[n] )

        # LPF
        _y_lpf_1.append( _filter_LPF_1.update( _x[n] ) )
        #_y_lpf_2.append( _filter_LPF_2.update( _x[n] ) )
        #_y_lpf_3.append( _filter_LPF_3.update( _x[n] ) )

        # HPF
        _y_hpf_1.append( _filter_HPF_1.update( _x[n] ) )
        _y_hpf_2.append( _filter_HPF_2.update( _x[n] ) )
        _y_hpf_3.append( _filter_HPF_3.update( _x[n] ) )




        # Down sample to SAMPLE_FREQ
        if _downsamp_cnt >= (( 1 / ( _dt * SAMPLE_FREQ )) - 1 ):
            _downsamp_cnt = 0

            # FIR 
            _y_fir.append( _filter_FIR.update( _x[n] ) )
            _time_fir.append( _time[n] )
        else:
            _downsamp_cnt += 1
    
    # Plot results
    fig, (ax_1, ax_2, ax_3) = plt.subplots(3, 1)
    fig.suptitle("Input signal freq: " + str(INPUT_SIGNAL_FREQ) + "Hz", fontsize=20)
    ax_1.plot( _time, _x, "b", label="input" )
    ax_1.plot( _time, _y_lpf_1, "g", label=  str(LPF_FC_1) + "Hz/" + str(LPF_ORDER_1))
    #ax_1.plot( _time, _y_lpf_2, "r", label=  str(LPF_FC_2) + "Hz/" + str(LPF_ORDER_2) )
    #ax_1.plot( _time, _y_lpf_3, "y", label=  str(LPF_FC_3) + "Hz/" + str(LPF_ORDER_3) )
    ax_1.grid()
    ax_1.title.set_text("RC Low Pass Filter")
    #ax_1.set_xlabel("Time [s]")
    ax_1.set_ylabel("Amplitude")
    ax_1.legend(loc="upper right")

    ax_2.plot( _time, _x, "b" )
    ax_2.plot( _time, _y_hpf_1, "g", label=  str(HPF_FC_1) + "Hz/" + str(HPF_ORDER_1))
    ax_2.plot( _time, _y_hpf_2, "r", label=  str(HPF_FC_2) + "Hz/" + str(HPF_ORDER_2))
    ax_2.plot( _time, _y_hpf_3, "y", label=  str(HPF_FC_3) + "Hz/" + str(HPF_ORDER_3))
    ax_2.grid()
    ax_2.title.set_text("CR High Pass Filter")
    #ax_2.set_xlabel("Time [s]")
    ax_2.set_ylabel("Amplitude")
    ax_2.legend(loc="upper right")

    
    ax_3.plot( _time, _x, "b" )
    ax_3.plot( _time_fir, _y_fir, "g")
    ax_3.grid()
    ax_3.title.set_text("FIR filter")
    ax_3.set_xlabel("Time [s]")
    ax_3.set_ylabel("Amplitude")
    #ax_3.legend(loc="upper right")
    

    plt.show()
    


# ===============================================================================
#       END OF FILE
# ===============================================================================
