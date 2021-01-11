# ===============================================================================
# @file:    rc_filter.py
# @note:    This script is evaluation of RC/CR filter algorithm
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

from filter_utils import FunctionGenerator, SignalMux

# ===============================================================================
#       CONSTANTS
# ===============================================================================

## ****** USER CONFIGURATIONS ******

## Sample frequency
#   Sample frequency of real system   
#
# Unit: Hz
SAMPLE_FREQ = 200.0

# Ideal sample frequency
#   As a reference to sample rate constrained embedded system
#
# Unit: Hz
IDEAL_SAMPLE_FREQ = 20000.0

## Time window
#
# Unit: second
TIME_WINDOW = 2.5

## Input signal shape
INPUT_SIGNAL_AMPLITUDE = 1.0
INPUT_SIGNAL_OFFSET = 0.0
INPUT_SIGNAL_PHASE = 0.0

## Mux input signal
INPUT_SIGNAL_SELECTION = SignalMux.MUX_CTRL_SINE

## Input signal frequency
#
# Unit: Hz
INPUT_SIGNAL_FREQ = 5

## LPF fc
#
# Unit: Hz
LPF_FC_1 = 5.0
LPF_FC_2 = 5.0
LPF_FC_3 = 5.0

## LPF order
LPF_ORDER_1 = 1
LPF_ORDER_2 = 2
LPF_ORDER_3 = 3

## HPF fc
#
# Unit: Hz
HPF_FC_1 = 10
HPF_FC_2 = 10
HPF_FC_3 = 10

## HPF order
HPF_ORDER_1 = 1
HPF_ORDER_2 = 2
HPF_ORDER_3 = 3

## ****** END OF USER CONFIGURATIONS ******

## Number of samples in time window
SAMPLE_NUM = int(( IDEAL_SAMPLE_FREQ * TIME_WINDOW ) + 1.0 )


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



# ===============================================================================
#       MAIN ENTRY
# ===============================================================================
if __name__ == "__main__":

    # Time array
    _time, _dt = np.linspace( 0.0, TIME_WINDOW, num=SAMPLE_NUM, retstep=True )

    # Filter object
    _filter_LPF_1   = RC_LPF( fc=LPF_FC_1, dt=_dt,              order=LPF_ORDER_1, init_val=0)
    _filter_D_LPF_1 = RC_LPF( fc=LPF_FC_1, dt=(1/SAMPLE_FREQ),  order=LPF_ORDER_1, init_val=0)
    _filter_LPF_2   = RC_LPF( fc=LPF_FC_2, dt=_dt,              order=LPF_ORDER_2, init_val=0)
    _filter_D_LPF_2 = RC_LPF( fc=LPF_FC_2, dt=(1/SAMPLE_FREQ),  order=LPF_ORDER_2, init_val=0)
    _filter_LPF_3   = RC_LPF( fc=LPF_FC_3, dt=_dt,              order=LPF_ORDER_3, init_val=0)
    _filter_D_LPF_3 = RC_LPF( fc=LPF_FC_3, dt=(1/SAMPLE_FREQ),  order=LPF_ORDER_3, init_val=0)
    
    _filter_HPF_1   = CR_HPF( fc=HPF_FC_1, dt=_dt,              order=HPF_ORDER_1)
    _filter_D_HPF_1 = CR_HPF( fc=HPF_FC_1, dt=(1/SAMPLE_FREQ),  order=HPF_ORDER_1)
    _filter_HPF_2   = CR_HPF( fc=HPF_FC_2, dt=_dt,              order=HPF_ORDER_2)
    _filter_D_HPF_2 = CR_HPF( fc=HPF_FC_2, dt=(1/SAMPLE_FREQ),  order=HPF_ORDER_2)
    _filter_HPF_3   = CR_HPF( fc=HPF_FC_3, dt=_dt,              order=HPF_ORDER_3)
    _filter_D_HPF_3 = CR_HPF( fc=HPF_FC_3, dt=(1/SAMPLE_FREQ),  order=HPF_ORDER_2)

    # Filter input/output
    _x = [ 0 ] * SAMPLE_NUM
    _x_d = [0]

    _y_lpf_1 = []
    _y_d_lpf_1 = [0]
    _y_lpf_2 = []
    _y_d_lpf_2 = [0]
    _y_lpf_3 = []
    _y_d_lpf_3 = [0]

    _y_hpf_1 = []
    _y_d_hpf_1 = [0]
    _y_hpf_2 = []
    _y_d_hpf_2 = [0]
    _y_hpf_3 = []
    _y_d_hpf_3 = [0]

    # Generate inputs
    _fg_sine = FunctionGenerator( INPUT_SIGNAL_FREQ, INPUT_SIGNAL_AMPLITUDE, INPUT_SIGNAL_OFFSET, INPUT_SIGNAL_PHASE, "sine" )
    _fg_rect = FunctionGenerator( INPUT_SIGNAL_FREQ, INPUT_SIGNAL_AMPLITUDE, INPUT_SIGNAL_OFFSET, INPUT_SIGNAL_PHASE, "rect" )
    _sin_x = []
    _rect_x = []

    # Signal mux
    # NOTE: Support only sine and rectange
    _signa_mux = SignalMux( 2 )

    # Down sample
    _downsamp_cnt = 0
    _downsamp_samp = [0]
    _d_time = [0]
    
    # Generate stimuli signals
    for n in range(SAMPLE_NUM):
        _sin_x.append( _fg_sine.generate( _time[n] ))
        _rect_x.append( _fg_rect.generate( _time[n] ))
 
    # Apply filter
    for n in range(SAMPLE_NUM):
        
        # Mux input signals
        _x[n] = _signa_mux.out( INPUT_SIGNAL_SELECTION, [ _sin_x[n], _rect_x[n] ] )

        # LPF
        _y_lpf_1.append( _filter_LPF_1.update( _x[n] ) )
        _y_lpf_2.append( _filter_LPF_2.update( _x[n] ) )
        _y_lpf_3.append( _filter_LPF_3.update( _x[n] ) )

        # HPF
        _y_hpf_1.append( _filter_HPF_1.update( _x[n] ) )
        _y_hpf_2.append( _filter_HPF_2.update( _x[n] ) )
        _y_hpf_3.append( _filter_HPF_3.update( _x[n] ) )

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

            # HPF
            _y_d_hpf_1.append( _filter_D_HPF_1.update( _x_d[-1] ) )
            _y_d_hpf_2.append( _filter_D_HPF_2.update( _x_d[-1] ) )
            _y_d_hpf_3.append( _filter_D_HPF_3.update( _x_d[-1] ) )

        else:
            _downsamp_cnt += 1
    
    # Plot results
    fig, ax = plt.subplots(2, 2, sharex="row", sharey="row")
    fig.suptitle("Input signal freq: " + str(INPUT_SIGNAL_FREQ) + "Hz", fontsize=20)
    
    ax[0,0].plot( _time, _x, "b", label="Input-generated" )
    ax[0,0].plot( _time, _y_lpf_1, "g", label="Ideal filter")
    ax[0,0].plot( _d_time, _downsamp_samp, "r.", label="Sample points")
    ax[0,0].plot( _time, _y_lpf_1, "g", label="RC1: " + str(LPF_FC_1) + "Hz/" + str(LPF_ORDER_1))
    ax[0,0].plot( _time, _y_lpf_2, "r", label="RC2: " + str(LPF_FC_2) + "Hz/" + str(LPF_ORDER_2))
    ax[0,0].plot( _time, _y_lpf_3, "y", label="RC3: " + str(LPF_FC_3) + "Hz/" + str(LPF_ORDER_3))
    ax[0,0].grid()
    ax[0,0].title.set_text( "RC Low Pass Filter \n'CONTINOUS' TIME DOMAIN \nfs=" + str( 1 / _dt ) + "Hz\n \
    fc1/fs=" + str(LPF_FC_1/IDEAL_SAMPLE_FREQ) + " fc2/fs=" + str(LPF_FC_2/IDEAL_SAMPLE_FREQ) + " fc3/fs=" + str(LPF_FC_3/IDEAL_SAMPLE_FREQ))
    ax[0,0].set_ylabel("Amplitude")
    ax[0,0].legend(loc="upper right")

    ax[0,1].plot( _d_time, _x_d, "b.-", label="Input-sampled")
    ax[0,1].plot( _d_time, _y_d_lpf_1, "g.-", label="RC1: " + str(LPF_FC_1) + "Hz/" + str(LPF_ORDER_1))
    ax[0,1].plot( _d_time, _y_d_lpf_2, "r.-", label="RC2: " + str(LPF_FC_2) + "Hz/" + str(LPF_ORDER_2))
    ax[0,1].plot( _d_time, _y_d_lpf_3, "y.-", label="RC3: " + str(LPF_FC_3) + "Hz/" + str(LPF_ORDER_3))
    ax[0,1].grid()
    ax[0,1].title.set_text( "RC Low Pass Filter \nDISCRETE TIME DOMAIN (on embedded system)\n fs=" + str(SAMPLE_FREQ) + "Hz\n \
    fc1/fs=" + str(LPF_FC_1/SAMPLE_FREQ) + " fc2/fs=" + str(LPF_FC_2/SAMPLE_FREQ) + " fc3/fs=" + str(LPF_FC_3/SAMPLE_FREQ))
    ax[0,1].legend(loc="upper right")

    ax[1,0].plot( _time, _x, "b" )
    ax[1,0].plot( _time, _y_hpf_1, "g", label="CR1: " + str(HPF_FC_1) + "Hz/" + str(HPF_ORDER_1))
    ax[1,0].plot( _time, _y_hpf_2, "r", label="CR2: " + str(HPF_FC_2) + "Hz/" + str(HPF_ORDER_2))
    ax[1,0].plot( _time, _y_hpf_3, "y", label="CR2: " + str(HPF_FC_3) + "Hz/" + str(HPF_ORDER_3))
    ax[1,0].grid()
    ax[1,0].title.set_text("CR High Pass Filter")
    ax[1,0].set_ylabel("Amplitude")
    ax[1,0].legend(loc="upper right")

    ax[1,1].plot( _d_time, _x_d, "b.-" )
    ax[1,1].plot( _d_time, _y_d_hpf_1, "g.-", label="CR1: " + str(HPF_FC_1) + "Hz/" + str(HPF_ORDER_1))
    ax[1,1].plot( _d_time, _y_d_hpf_2, "r.-", label="CR2: " + str(HPF_FC_2) + "Hz/" + str(HPF_ORDER_2))
    ax[1,1].plot( _d_time, _y_d_hpf_3, "y.-", label="CR2: " + str(HPF_FC_3) + "Hz/" + str(HPF_ORDER_3))
    ax[1,1].grid()
    ax[1,1].title.set_text("CR High Pass Filter")
    ax[1,1].set_ylabel("Amplitude")
    ax[1,1].legend(loc="upper right")
    ax[1,1].set_xlabel("Time [s]")
    
    plt.subplots_adjust(left=0.05, right=0.98, bottom=0.05, wspace=0.08)

    plt.show()
    


# ===============================================================================
#       END OF FILE
# ===============================================================================
