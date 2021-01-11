# ===============================================================================
# @file:    washout_filter.py
# @note:    This script is evaluation of washout filter algorithm
# @author:  Ziga Miklosic
# @date:    11.01.2021
# @brief:   Evaluation of washout filter design. This evaluation is for designing
#           a working washout filter used for steward platform. 
# ===============================================================================

# ===============================================================================
#       IMPORTS  
# ===============================================================================
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import freqz, butter, cheby1, lfilter, filtfilt, bilinear

from filter_utils import FunctionGenerator, SignalMux
from iir_filter import IIR, calculate_2nd_order_HPF_coeff, calculate_1nd_order_HPF_coeff, calculate_2nd_order_LPF_coeff

# ===============================================================================
#       CONSTANTS
# ===============================================================================

## ****** USER CONFIGURATIONS ******

## Sample frequency
#   Sample frequency of real system   
#
# Unit: Hz
SAMPLE_FREQ = 100.0

# Ideal sample frequency
#   As a reference to sample rate constrained embedded system
#
# Unit: Hz
IDEAL_SAMPLE_FREQ = 20000.0

## Time window
#
# Unit: second
TIME_WINDOW = 2

## Input signal shape
INPUT_SIGNAL_AMPLITUDE = 2.0
INPUT_SIGNAL_OFFSET = 0.0
INPUT_SIGNAL_PHASE = -1.57

## Mux input signal
INPUT_SIGNAL_SELECTION = SignalMux.MUX_CTRL_RECT

## Input signal frequency
#
# Unit: Hz
INPUT_SIGNAL_FREQ = 0.5

## Number of samples in time window
SAMPLE_NUM = int(( IDEAL_SAMPLE_FREQ * TIME_WINDOW ) + 1.0 )


# =====================================================
## TRANSLATION CHANNEL SETTINGS

# HPF Wht 2nd order filter
WASHOUT_HPF_WHT_FC  = 1.0
WASHOUT_HPF_WHT_Z   = 0.7071

# HPF Wrtzt 1st order filter
WASHOUT_HPF_WRTZT_FC  = 1.0

# SCALE AND LIMIT
# [ x, y, z]
WASHOUT_SCALE_A_T = [ 1.0, 1.0, 1.0 ]
WASHOUT_LIMIT_A_T = [ 2.0, 2.0, 2.0 ] # Limits are symetrical


# =====================================================
## COORDINATION CHANNEL SETTINGS

# LPF W12 2nd order filter
WASHOUT_LPF_W12_FC  = 10.0
WASHOUT_LPF_W12_Z   = 0.7071

# SCALE AND LIMIT 
# [ x, y, z]
WASHOUT_SCALE_A_C = [ 1.0, 1.0, 1.0 ]
WASHOUT_LIMIT_A_C = [ 2.0, 2.0, 2.0 ] # Limits are symetrical

# Gravity constant
G = 9.18
G_INV = 1.0 / G

# TILT MATRIX
WASHOUT_TILT_MATRIX = [ [0,     -G_INV,  0],
                        [G_INV, 0,       0],
                        [0,     0,       0] ]

# TILT LIMIT
WASHOUT_TILT_LIMIT = [ 0.2, 0.2, 0.2 ]


# =====================================================
## ROTATION CHANNEL SETTINGS

# HPF W11 1st order filter
WASHOUT_HPF_W11_FC  = 1.0

# SCALE AND LIMIT
# [ rool, pitch, yaw]
WASHOUT_SCALE_BETA = [ 1.0, 1.0, 1.0 ]
WASHOUT_LIMIT_BETA = [ 0.2, 0.2, 0.2 ]


## ****** END OF USER CONFIGURATIONS ******

# ===============================================================================
#       FUNCTIONS
# ===============================================================================

# ===============================================================================
#       CLASSES
# ===============================================================================    

## IIR Filter
class Washout:

    def __init__(self, Wht, Wrtzt, W11, W12, fs):

        # Translation channel filters
        b, a = calculate_2nd_order_HPF_coeff( Wht[0], Wht[1], fs )
        self._hpf_wht = [0] * 3
        self._hpf_wht[0] = IIR( a, b, 2 )
        self._hpf_wht[1] = IIR( a, b, 2 )
        self._hpf_wht[2] = IIR( a, b, 2 )

        b, a = calculate_1nd_order_HPF_coeff( Wrtzt, fs )
        self._hpf_wrtzt = [0] * 3
        self._hpf_wrtzt[0] = IIR( a, b, 1 )
        self._hpf_wrtzt[1] = IIR( a, b, 1 )
        self._hpf_wrtzt[2] = IIR( a, b, 1 )

        # Coordination channel filters
        b, a = calculate_2nd_order_LPF_coeff( W12[0], W12[1], fs )
        self._lpf_w12 = [0] * 3
        self._lpf_w12[0] = IIR( a, b, 2 )
        self._lpf_w12[1] = IIR( a, b, 2 )
        self._lpf_w12[2] = IIR( a, b, 2 )

        # Rotation channel filters
        b, a = calculate_1nd_order_HPF_coeff( W11, fs )
        self._hpf_w22 = [0] * 3
        self._hpf_w22[0] = IIR( a, b, 1 )
        self._hpf_w22[1] = IIR( a, b, 1 )
        self._hpf_w22[2] = IIR( a, b, 1 )


    # ===============================================================================
    # @brief: Update washout filter
    #
    # @param[in]:    a      - Vector of accelerations
    # @param[in]:    beta   - Vector of angular velocities
    # @return:       p, r   - Positions and rotation of a steward platform
    # ===============================================================================
    def update(self, a, beta):

        # Translation channel
        a_t = [0] * 3

        # Coordination channel
        a_c = [0] * 3

        # Rotations channel
        beta_r = [0] * 3

        # Loop all axis
        for n in range(3):

            # Apply scaling and limitations
            a_t[n]  = self.__scale_limit( a[n], WASHOUT_SCALE_A_T[n], WASHOUT_LIMIT_A_T[n] )
            a_c[n]  = self.__scale_limit( a[n], WASHOUT_SCALE_A_C[n], WASHOUT_LIMIT_A_C[n] )
            beta_r[n] = self.__scale_limit( beta[n], WASHOUT_SCALE_BETA[n], WASHOUT_LIMIT_BETA[n] )

            # Translation filters
            a_t[n] = self._hpf_wht[n].update( a_t[n] )
            a_t[n] = self._hpf_wrtzt[n].update( a_t[n] )

            # Coordiation filters
            a_c[n] = self._lpf_w12[n].update( a_c[n] )

            # Rotational filter
            beta_r[n] = self._hpf_w22[n].update( beta_r[n] )

            # Tilt coordination matrix
            # NOTE: This operation convert acceleration to angle
            _a_c = a_c
            for j in range(3):
                #a_c[n] = WASHOUT_TILT_MATRIX[n][j] * _a_c[j]
                pass

            # Tilt rate limiter
            # TODO: Implement rate limiter not limiter!!!
            a_c[n] = self.__scale_limit( a_c[n], 1.0,  WASHOUT_TILT_LIMIT[n] )

            # Add tilt to rotation channel
            beta_r[n] += a_c[n]

        return a_t, beta_r


    # ===============================================================================
    # @brief: Scale and limit input signals
    #
    #   NOTE: For know limiting is very simple, shall be change to poly order 3!!!
    #
    # @param[in]:    x      - Input signal value
    # @param[in]:    scale  - Scale factor
    # @param[in]:    lim    - Limit factor
    # @return:       y      - Scaled and limited value
    # ===============================================================================
    def __scale_limit(self, x, scale, lim):

        y = scale * x

        if y > lim:
            y = lim
        elif y < -lim:
            y = -lim
        else:
            pass

        return y


# ===============================================================================
#       MAIN ENTRY
# ===============================================================================
if __name__ == "__main__":

    # Time array
    _time, _dt = np.linspace( 0.0, TIME_WINDOW, num=SAMPLE_NUM, retstep=True )

    
    # Filter object
    _filter_washout = Washout(  Wht=[WASHOUT_HPF_WHT_FC, WASHOUT_HPF_WHT_Z], Wrtzt=WASHOUT_HPF_WRTZT_FC, \
                                W11=WASHOUT_HPF_W11_FC, W12=[WASHOUT_LPF_W12_FC, WASHOUT_LPF_W12_Z], fs=SAMPLE_FREQ )

    # Filter input/output
    _x = [ 0 ] * SAMPLE_NUM
    _x_d = [0]

    # Position
    _y_d_p = [[0], [0], [0]] * 3
    
    # Rotation
    _y_d_r = [[0], [0], [0]] * 3


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

        # Down sample to SAMPLE_FREQ
        if _downsamp_cnt >= (( 1 / ( _dt * SAMPLE_FREQ )) - 1 ):
            _downsamp_cnt = 0

            # Utils
            _downsamp_samp.append(0)
            _d_time.append( _time[n])
            _x_d.append( _x[n] )


            p, r = _filter_washout.update( [ _x[n], 0, 0 ], [ 0, 0, 0 ] )
            _y_d_p[0].append( p[0] )
            _y_d_p[1].append( p[1] )
            _y_d_p[2].append( p[2] )
            _y_d_r[0].append( r[0] )
            _y_d_r[1].append( r[1] )
            _y_d_r[2].append( r[2] )


        else:
            _downsamp_cnt += 1
    
    # Plot results
    fig, ax = plt.subplots(2, 1)
    fig.suptitle( "WASHOUT FILTERS", fontsize=20 )

    ax[0].set_title("Translations", fontsize=16)
    ax[0].plot( _time, _x,                  "b",    label="Input-generated" )
    #ax[0].plot( _d_time, _downsamp_samp,    "r.",   label="Sample points")
    ax[0].plot( _d_time, _y_d_p[0],         "g.-",    label="ax")
    ax[0].plot( _d_time, _y_d_p[1],         "r.-",    label="ay")
    ax[0].plot( _d_time, _y_d_p[2],         "y.-",    label="az")
    
    ax[0].set_ylabel('Amplitude')
    ax[0].set_xlabel('Time [s]')
    ax[0].grid()
    ax[0].legend(loc="upper right")

    ax[1].set_title("Rotations", fontsize=16)
    ax[1].plot( _time, _x,                  "b",    label="Input-generated" )
    #ax[1].plot( _d_time, _downsamp_samp,    "r.",   label="Sample points")
    ax[1].plot( _d_time, _y_d_r[0],         "g.-",    label="roll")
    ax[1].plot( _d_time, _y_d_r[1],         "r.-",    label="pitch")
    ax[1].plot( _d_time, _y_d_r[2],         "y.-",    label="yaw")

    ax[1].set_ylabel('Amplitude')
    ax[1].set_xlabel('Time [s]')
    ax[1].grid()
    ax[1].legend(loc="upper right")

    plt.show()
    


# ===============================================================================
#       END OF FILE
# ===============================================================================
