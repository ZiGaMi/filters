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
TIME_WINDOW = 6

## Input signal shape
INPUT_SIGNAL_AMPLITUDE = 1.0
INPUT_SIGNAL_OFFSET = 0.0
INPUT_SIGNAL_PHASE = 0.0

## Mux input signal
INPUT_SIGNAL_SELECTION = SignalMux.MUX_CTRL_SINE

## Input signal frequency
#
# Unit: Hz
INPUT_SIGNAL_FREQ = 1

## Number of samples in time window
SAMPLE_NUM = int(( IDEAL_SAMPLE_FREQ * TIME_WINDOW ) + 1.0 )

## Cutoff freqeuncy of filter
#
# Unit: Hz
HPF_FC_1 = 1.0
HPF_FC_2 = 1.0

HPF_FC_BUTTER = 1.0
HPF_FC_CHEBY = 1.0

LPF_FC_1 = 1.0
LPF_FC_2 = 1.0

LPF_FC_BUTTER = 1.0
LPF_FC_CHEBY = 1.0

## Damping factor (zeta) for 2nd order filter
#
#   z = 0       -> underdamped
#   z = 0.7071  -> criticaly damped, sweep spot
#   z = 1       -> overdamped
HPF_Z_1 = .25
HPF_Z_2 = 0.701

LPF_Z_1 = .25
LPF_Z_2 = 1.0


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
WASHOUT_LIMIT_A_T = [ 4.0, 4.0, 4.0 ]


# =====================================================
## COORDINATION CHANNEL SETTINGS

# LPF W12 2nd order filter
WASHOUT_LPF_W12_FC  = 1.0
WASHOUT_LPF_W12_Z   = 0.7071

# SCALE AND LIMIT
# [ x, y, z]
WASHOUT_SCALE_A_C = [ 1.0, 1.0, 1.0 ]
WASHOUT_LIMIT_A_C = [ 4.0, 4.0, 4.0 ]


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
        self._hpf_w12 = [0] * 3
        self._hpf_w12[0] = IIR( a, b, 2 )
        self._hpf_w12[1] = IIR( a, b, 2 )
        self._hpf_w12[2] = IIR( a, b, 2 )

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
        p = []
        r = []

        # Translation accelerations
        a_t = [0] * 3

        # Coordination accelerations
        a_c = [0] * 3

        # Rotations
        beta = [0] * 3

        # Apply scaling and limitations
        for n in range(3):
            a_t[n]  = self.__scale_limit( a[n], WASHOUT_SCALE_A_T[n], WASHOUT_LIMIT_A_T[n] )
            a_c[n]  = self.__scale_limit( a[n], WASHOUT_SCALE_A_C[n], WASHOUT_LIMIT_A_C[n] )
            beta[n] = self.__scale_limit( a[n], WASHOUT_SCALE_BETA[n], WASHOUT_LIMIT_BETA[n] )

        # Translation filters



        return p, r


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

    print("lpf_a: %s" % lpf_a)
    print("lpf_b: %s" % lpf_b)
    print("a_b_lpf: %s" % a_b_lpf)
    print("b_b_lpf: %s" % b_b_lpf)

    # Filter input/output
    _x = [ 0 ] * SAMPLE_NUM
    _x_d = [0]

    _y_d_iir = [0]
    _y_d_iir_2 = [0]

    _y_d_iir_lpf = [0]
    _y_d_iir_lpf_2 = [0]
    
    _y_d_iir_notch = [0]
    _y_d_iir_notch_2 = [0]

    # Generate inputs
    _fg_sine = FunctionGenerator( INPUT_SIGNAL_FREQ, INPUT_SIGNAL_AMPLITUDE, INPUT_SIGNAL_OFFSET, INPUT_SIGNAL_PHASE, "sine" )
    _fg_ac_noise = FunctionGenerator( 50.0, 0.2, 0, 0, "sine" )
    _fg_rect = FunctionGenerator( INPUT_SIGNAL_FREQ, INPUT_SIGNAL_AMPLITUDE, INPUT_SIGNAL_OFFSET, INPUT_SIGNAL_PHASE, "rect" )
    _sin_x = []
    _rect_x = []
    _ac_power_noise = []

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
        _ac_power_noise.append( _sin_x[-1] +  _fg_ac_noise.generate( _time[n] ) )
 
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

            # HPF
            _y_d_iir.append( _filter_IIR.update( _x[n] ))
            _y_d_iir_2.append( _filter_IIR_2.update( _x[n] ))

            # LPF
            _y_d_iir_lpf.append( _filter_IIR_LPF.update( _x[n] ))
            _y_d_iir_lpf_2.append( _filter_IIR_LPF_2.update( _x[n] ))

            # Notch 
            _y_d_iir_notch.append( _filter_IIR_NOTCH.update( _ac_power_noise[n] ))
            _y_d_iir_notch_2.append( _filter_IIR_NOTCH_2.update( _ac_power_noise[n] ))

        else:
            _downsamp_cnt += 1
    
    # Plot results
    fig, ax = plt.subplots(2, 1)
    fig.suptitle("Highpass 2nd order IIR Filter Design\nInput signal freq: " + str(INPUT_SIGNAL_FREQ) + "Hz, Sample freq (fs): " +str(SAMPLE_FREQ) + "Hz", fontsize=20)
    
    ax[0].plot( _time, _x,                  "b",    label="Input-generated" )
    ax[0].plot( _d_time, _downsamp_samp,    "r.",   label="Sample points")
    ax[0].plot( _d_time, _y_d_iir,          ".-g",    label=str(HPF_FC_1) + "Hz/" + str(HPF_Z_1))
    ax[0].plot( _d_time, _y_d_iir_2,        ".-y",    label=str(HPF_FC_2) + "Hz/" + str(HPF_Z_2))
    
    ax[0].set_ylabel('Amplitude')
    ax[0].set_xlabel('Time [s]')
    ax[0].grid()
    ax[0].title.set_text("Time domain")
    ax[0].set_ylabel("Amplitude")
    ax[0].legend(loc="upper right")
    ax[0].set_ylabel("Amplitude")

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
    
    # uncomment if neede phase delay
    """
    ax_11 = ax[1].twinx()
    angles = np.unwrap( np.angle(h) )
    angles_b = np.unwrap( np.angle(h_b) )
    ax_11.plot(w, (angles*180/np.pi), 'b')
    ax_11.plot(w_b, (angles_b*180/np.pi), 'r')
    ax_11.set_ylabel('Angle [degrees]', color='g')
    ax_11.axis('tight')
    """

    fig2, ax2 = plt.subplots(2, 1)
    fig2.suptitle("Lowpass 2nd order IIR Filter Design\nInput signal freq: " + str(INPUT_SIGNAL_FREQ) + "Hz, Sample freq (fs): " +str(SAMPLE_FREQ) + "Hz", fontsize=20)
    
    ax2[0].plot( _time, _x,                  "b",    label="Input-generated" )
    ax2[0].plot( _d_time, _downsamp_samp,    "r.",   label="Sample points")
    ax2[0].plot( _d_time, _y_d_iir_lpf,      ".-g",  label=str(LPF_FC_1) + "Hz/" + str(LPF_Z_1))
    ax2[0].plot( _d_time, _y_d_iir_lpf_2,    ".-y",  label=str(LPF_FC_2) + "Hz/" + str(LPF_Z_2))
    
    ax2[0].set_ylabel('Amplitude')
    ax2[0].set_xlabel('Time [s]')
    ax2[0].grid()
    ax2[0].title.set_text("Time domain")
    ax2[0].set_ylabel("Amplitude")
    ax2[0].legend(loc="upper right")

    w_b_lpf = ( w_b_lpf / np.pi * SAMPLE_FREQ / 2)
    w_c_lpf = ( w_c_lpf / np.pi * SAMPLE_FREQ / 2)
    
    w_lpf = ( w_lpf / np.pi * SAMPLE_FREQ / 2)
    w_lpf_2 = ( w_lpf_2 / np.pi * SAMPLE_FREQ / 2)

    ax2[1].plot(w_lpf, 20 * np.log10(abs(h_lpf)),       'g', label=str(LPF_FC_1) + "Hz/" + str(LPF_Z_1) )
    ax2[1].plot(w_lpf_2, 20 * np.log10(abs(h_lpf_2)),   'y', label=str(LPF_FC_2) + "Hz/" + str(LPF_Z_2) )
    ax2[1].plot(w_b_lpf, 20 * np.log10(abs(h_b_lpf)),   'r', label="butterworth" )
    ax2[1].plot(w_c_lpf, 20 * np.log10(abs(h_c_lpf)),   'b', label="chebysev" )

    ax2[1].set_ylabel('Amplitude [dB]')
    ax2[1].set_xlabel('Frequency [Hz]')
    ax2[1].grid()
    ax2[1].set_ylim(-80, 20)
    ax2[1].set_xlim(0.001, SAMPLE_FREQ/4)
    ax2[1].set_xscale("log")
    ax2[1].legend(loc="upper right")


    fig3, ax3 = plt.subplots(2, 1)
    fig3.suptitle("Notch 2nd order IIR Filter Design\nInput signal freq: " + str(INPUT_SIGNAL_FREQ) + "Hz, Sample freq (fs): " +str(SAMPLE_FREQ) + "Hz", fontsize=20)
   
    w_notch = ( w_notch / np.pi * SAMPLE_FREQ / 2)
    w_notch_2 = ( w_notch_2 / np.pi * SAMPLE_FREQ / 2)


    ax3[0].plot( _time, _ac_power_noise,     "b",    label="Input-generated" )
    ax3[0].plot( _d_time, _downsamp_samp,    "r.",   label="Sample points")
    ax3[0].plot( _d_time, _y_d_iir_notch,    ".-g",  label="notch1")
    #ax3[0].plot( _d_time, _y_d_iir_notch_2,  ".-y",  label="notch2")
    
    ax3[0].set_ylabel('Amplitude')
    ax3[0].set_xlabel('Time [s]')
    ax3[0].grid()
    ax3[0].legend(loc="upper right")

    ax3[1].plot(w_notch, 20 * np.log10(abs(h_notch)),  'g')
    ax3[1].plot(w_notch_2, 20 * np.log10(abs(h_notch_2)),  'r')

    ax3[1].set_ylabel('Amplitude [dB]')
    ax3[1].set_xlabel('Frequency [Hz]')
    ax3[1].grid()
    #ax3[1].set_ylim(-80, 20)
    #ax3[1].set_xlim(0.001, SAMPLE_FREQ/4)
    #ax3[1].set_xscale("log")
    ax3[1].legend(loc="upper right")

    plt.show()
    


# ===============================================================================
#       END OF FILE
# ===============================================================================
