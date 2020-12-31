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

## Time window
#
# Unit: second
TIME_WINDOW = 20.0

## Number of samples in time window
SAMPLE_NUM = 10000

## Select input filter signal type
INPUT_SINE = 0
INPUT_RECT = 1

## Mux input signal
INPUT_SIGNAL_SELECTION = INPUT_RECT

## Input signal frequency
#
# Unit: Hz
INPUT_SIGNAL_FREQ = 0.1

## LPF fc
#
# Unit: Hz
LPF_FC_1 = 1.0
LPF_FC_2 = 2.0
LPF_FC_3 = 3.0

## LPF order
LPF_ORDER_1 = 1
LPF_ORDER_2 = 2
LPF_ORDER_3 = 3

## HPF fc
#
# Unit: Hz
HPF_FC_1 = 1.0
HPF_FC_2 = 2.0
HPF_FC_3 = 3.0

## HPF order
HPF_ORDER_1 = 1
HPF_ORDER_2 = 2
HPF_ORDER_3 = 3


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
        self.k = dt / ( dt + 1/fc )

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


## RC High Pass Filter
class RC_HPF:

    def __init__(self, fc, dt, order=1, init_val=0):

        # Calculate filter coefficient
        self.k = (1/fc) / ( dt + 1/fc )

        # Store order
        self.order = order

        # Init inputs & outputs
        self.y = [init_val] * order
        self.x = [init_val] * order

    def update(self, x):

        for n in range( self.order ):
            if n == 0:
                self.y[n] = ( self.k * self.y[n] + self.k * ( x - self.x[0] ))
                self.x[0] = x
            else:
                self.y[n] = ( self.k * self.y[n] + self.k * ( self.y[n-1] - self.x[n] ))
                self.x[n] = self.y[n-1]

        return self.y[ self.order - 1 ]




# ===============================================================================
#       MAIN ENTRY
# ===============================================================================


if __name__ == "__main__":

    # Time array
    _time, _dt = np.linspace( 0.0, TIME_WINDOW, num=SAMPLE_NUM, retstep=True )

    # Filter object
    _filter_LPF_1   = RC_LPF( fc=LPF_FC_1, dt=_dt, order=LPF_ORDER_1, init_val=0)
    _filter_LPF_2   = RC_LPF( fc=LPF_FC_2, dt=_dt, order=LPF_ORDER_2, init_val=0)
    _filter_LPF_3   = RC_LPF( fc=LPF_FC_3, dt=_dt, order=LPF_ORDER_3, init_val=0)
    _filter_HPF_1   = RC_HPF( fc=HPF_FC_1, dt=_dt, order=HPF_ORDER_1, init_val=0)
    _filter_HPF_2   = RC_HPF( fc=HPF_FC_2, dt=_dt, order=HPF_ORDER_2, init_val=0)
    _filter_HPF_3   = RC_HPF( fc=HPF_FC_3, dt=_dt, order=HPF_ORDER_3, init_val=0)

    # Filter input/output
    _x = [ 0 ] * SAMPLE_NUM
    _y_lpf_1 = []
    _y_lpf_2 = []
    _y_lpf_3 = []
    _y_hpf_1 = []
    _y_hpf_2 = []
    _y_hpf_3 = []

    # Generate inputs
    _sin_x = []
    _rect_x = []
    
    for n in range(SAMPLE_NUM):
        _sin_x.append( generate_sine( _time[n], INPUT_SIGNAL_FREQ, 1.0, 0.0, 0.0 ))  
        _rect_x.append( generate_rect( _time[n], INPUT_SIGNAL_FREQ, 1.0, 0.0, 0.0 ))  
 
    # Apply filter
    for n in range(SAMPLE_NUM):
        
        # Mux input signals
        _x[n] = input_signal_mux( INPUT_SIGNAL_SELECTION, _sin_x[n], _rect_x[n] )

        # LPF
        _y_lpf_1.append( _filter_LPF_1.update( _x[n] ) )
        _y_lpf_2.append( _filter_LPF_2.update( _x[n] ) )
        _y_lpf_3.append( _filter_LPF_3.update( _x[n] ) )

        # HPF
        _y_hpf_1.append( _filter_HPF_1.update( _x[n] ) )
        _y_hpf_2.append( _filter_HPF_2.update( _x[n] ) )
        _y_hpf_3.append( _filter_HPF_3.update( _x[n] ) )

    
    # Plot results
    fig, (ax_1, ax_2) = plt.subplots(2, 1)
    ax_1.plot( _time, _x, "b", label="input" )
    ax_1.plot( _time, _y_lpf_1, "g", label=  str(LPF_FC_1) + "Hz/" + str(LPF_ORDER_1))
    ax_1.plot( _time, _y_lpf_2, "r", label=  str(LPF_FC_2) + "Hz/" + str(LPF_ORDER_2) )
    ax_1.plot( _time, _y_lpf_3, "y", label=  str(LPF_FC_3) + "Hz/" + str(LPF_ORDER_3) )
    ax_1.grid()
    ax_1.title.set_text("RC Low Pass Filter (input freq:" + str(INPUT_SIGNAL_FREQ) + "Hz)")
    ax_1.legend(loc="upper right")

    ax_2.plot( _time, _x, "b" )
    ax_2.plot( _time, _y_hpf_1, "g", label=  str(HPF_FC_1) + "Hz/" + str(HPF_ORDER_1))
    ax_2.plot( _time, _y_hpf_2, "r", label=  str(HPF_FC_2) + "Hz/" + str(HPF_ORDER_2))
    ax_2.plot( _time, _y_hpf_3, "g", label=  str(HPF_FC_3) + "Hz/" + str(HPF_ORDER_3))
    ax_2.grid()
    ax_2.title.set_text("CR High Pass Filter (input freq:" + str(INPUT_SIGNAL_FREQ) + "Hz)")
    ax_2.legend(loc="upper right")

    plt.show()
    


# ===============================================================================
#       END OF FILE
# ===============================================================================
