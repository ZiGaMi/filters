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
from scipy.signal import freqz, butter, cheby1, lfilter, filtfilt, bilinear

from filter_utils import FunctionGenerator, SignalMux, CircBuffer

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

LPF_FC_1 = .5
LPF_FC_2 = 1.0

LPF_FC_BUTTER = 1.0
LPF_FC_CHEBY = 1.0

## Damping factor (zeta) for 2nd order filter
#
#   z = 0       -> underdamped
#   z = 0.7071  -> criticaly damped, sweep spot
#   z = 1       -> overdamped
HPF_Z_1 = 3
HPF_Z_2 = 0.701

LPF_Z_1 = 3.0
LPF_Z_2 = 10.0

## ****** END OF USER CONFIGURATIONS ******

# ===============================================================================
#       FUNCTIONS
# ===============================================================================

# ===============================================================================
# @brief:   calculate 2nd order high pass filter based on following 
#           transfer function:
#           
#               h(s) = s^2 / ( s^2 + 2*z*w*s + w^2 )  --- bilinear ---> h(z)
#
# @param[in]:    fc     - Corner frequenc
# @param[in]:    z      - Damping factor
# @param[in]:    fs     - Sample frequency
# @return:       b,a    - Array of b,a IIR coefficients
# ===============================================================================
def calculate_2nd_order_HPF_coeff(fc, z, fs):
    
    # Calculate omega
    w = 2*np.pi*fc

    # Make bilinear transformation
    b, a = bilinear( [1,0,0], [1,2*z*w,w**2], fs )

    return b, a


# ===============================================================================
# @brief:   calculate 2nd order high pass filter based on following 
#           transfer function:
#           
#   Equations were taken from: https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
#
# @param[in]:    fc     - Corner frequenc
# @param[in]:    z      - Damping factor
# @param[in]:    fs     - Sample frequency
# @return:       b,a    - Array of b,a IIR coefficients
# ===============================================================================
def calculate_biquad_hpf_coeff(fc, z, fs):

    # Calculate omega
    w = 2*np.pi*fc/fs

    # Calculate quality factor 
    Q = 1 / ( 2.0 * z )

    # Calculate alpha
    alpha = np.sin( w ) / ( 2.0 * Q )

    # Calculate factors
    cos_w = np.cos(w)

    b0 = ( 1 + cos_w) / 2.0
    b1 = -( 1 + cos_w )
    b2 = ( 1 + cos_w ) / 2.0

    a0 = 1 + alpha
    a1 = -2.0 * cos_w
    a2 = 1 - alpha

    return [b0,b1,b2], [a0,a1,a2]

# ===============================================================================
# @brief:   calculate 1nd order high pass filter based on following 
#           transfer function:
#           
#               h(s) = s / ( s + w )  --- bilinear ---> h(z)
#
# @param[in]:    fc     - Corner frequenc
# @param[in]:    fs     - Sample frequency
# @return:       b,a    - Array of b,a IIR coefficients
# ===============================================================================
def calculate_1nd_order_HPF_coeff(fc, fs):
    
    # Calculate omega
    w = 2*np.pi*fc

    # Make bilinear transformation
    b, a = bilinear( [1,0], [1,w], fs )

    return b, a

# ===============================================================================
# @brief:   calculate 2nd order low pass filter based on following 
#           transfer function:
#           
#               h(s) = w^2 / ( s^2 + 2*z*w*s + w^2 ) --- bilinear ---> h(z)
#
# @param[in]:    fc     - Corner frequenc
# @param[in]:    z      - Damping factor
# @param[in]:    fs     - Sample frequency
# @return:       b,a    - Array of b,a IIR coefficients
# ===============================================================================
def calculate_2nd_order_LPF_coeff(fc, z, fs):

    # Calculate omega
    w = 2*np.pi*fc

    # Using bilinear transformation
    b, a = bilinear( [0,0,w**2], [1,2*z*w,w**2], fs )

    return b, a

# ===============================================================================
# @brief:   calculate 2nd order low pass filter based on following 
#           transfer function:
#          
#   Equations were taken from: https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
#
# @param[in]:    fc     - Corner frequenc
# @param[in]:    z      - Damping factor
# @param[in]:    fs     - Sample frequency
# @return:       b,a    - Array of b,a IIR coefficients
# ===============================================================================
def calculate_biquad_lpf_coeff(fc, z, fs):

    # Calculate omega
    w = 2*np.pi*fc/fs

    # Calculate quality factor 
    Q = 1 / ( 2.0 * z )

    # Calculate alpha
    alpha = np.sin( w ) / ( 2.0 * Q )

    # Calculate factors
    cos_w = np.cos(w)

    b0 = ( 1 - cos_w) / 2.0
    b1 = 1 - cos_w
    b2 = ( 1 - cos_w ) / 2.0

    a0 = 1 + alpha
    a1 = -2.0 * cos_w
    a2 = 1 - alpha

    return [b0,b1,b2], [a0,a1,a2]


# ===============================================================================
# @brief:   calculate 2nd order notch filter. This code is from 
#           "Second-order IIR Notch Filter Design and implementation of digital
#           signal processing system" acticle
#
# @param[in]:    fc     - Corner frequenc
# @param[in]:    fs     - Sample frequency
# @return:       b,a    - Array of b,a IIR coefficients
# ===============================================================================
def calculate_2nd_order_notch_coeff(fc, fs, r):
   
    _w = 2 * np.pi * fc / fs

    # Calculate coefficient
    a2 = r*r
    a1 = -2*r*np.cos( _w ) 
    a0 = 1
    
    b2 = 1
    b1 = -2*np.cos( _w )
    b0 = 1
    
    # Fill array
    a = [ a0, a1, a2 ] 
    b = [ b0, b1, b2 ] 

    return b, a

# ===============================================================================
#       CLASSES
# ===============================================================================    

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
        self.x.set( x )

        # Get input/outputs history
        _x = self.x.get_time_ordered_samples()
        _y = self.y.get_time_ordered_samples()

        # Calculate new value
        y = 0.0
        for j in range(self.order+1):

            y = y + (self.b[j] * _x[j])

            if j > 0:
                y = y - (self.a[j] * _y[j-1])

        y = ( y * ( 1 / self.a[0] ))

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
    
    # Calculate coefficient for 2nd order LPF
    lpf_b, lpf_a        = calculate_2nd_order_LPF_coeff( LPF_FC_1, LPF_Z_1, SAMPLE_FREQ )
    lpf_b_2, lpf_a_2    = calculate_2nd_order_LPF_coeff( LPF_FC_2, LPF_Z_2, SAMPLE_FREQ )

    # Calculate coefficient for 2nd order notch filter
    notch_b, notch_a   = calculate_2nd_order_notch_coeff( 10.0, SAMPLE_FREQ, r=0.925 )
    notch_b_2, notch_a_2   = calculate_2nd_order_notch_coeff( 10.0, SAMPLE_FREQ, r=0.90 )

    # Calculate coefficient for 2nd order Butterworth & Chebyshev filter
    b_b, a_b = butter( N=2, Wn=HPF_FC_BUTTER, btype="highpass", analog=False, fs=SAMPLE_FREQ )
    b_c, a_c = cheby1( N=2, Wn=HPF_FC_CHEBY, btype="highpass", analog=False, fs=SAMPLE_FREQ, rp=0.99 )

    b_b_lpf, a_b_lpf = butter( N=2, Wn=LPF_FC_BUTTER, btype="lowpass", analog=False, fs=SAMPLE_FREQ )
    b_c_lpf, a_c_lpf = cheby1( N=2, Wn=LPF_FC_CHEBY, btype="lowpass", analog=False, fs=SAMPLE_FREQ, rp=0.99 )

    # Get frequency characteristics
    w, h = freqz( b, a, 4096 )
    w_2, h_2 = freqz( b_2, a_2, 4096 )
    
    w_b, h_b = freqz( b_b, a_b, 4096 )
    w_c, h_c = freqz( b_c, a_c, 4096 )

    w_b_lpf, h_b_lpf = freqz( b_b_lpf, a_b_lpf, 4096 )
    w_c_lpf, h_c_lpf = freqz( b_c_lpf, a_c_lpf, 4096 )
    
    w_lpf, h_lpf = freqz( lpf_b, lpf_a, 4096 )
    w_lpf_2, h_lpf_2 = freqz( lpf_b_2, lpf_a_2, 4096 )

    w_notch, h_notch = freqz( notch_b, notch_a, 4096 )
    w_notch_2, h_notch_2 = freqz( notch_b_2, notch_a_2, 4096 )
    
    # Filter object
    _filter_IIR     = IIR( a, b, order=2 ) 
    _filter_IIR_2   = IIR( a_2, b_2, order=2 ) 

    _filter_IIR_LPF   = IIR( lpf_a, lpf_b, order=2 ) 
    _filter_IIR_LPF_2 = IIR( lpf_a_2, lpf_b_2, order=2 ) 

    _filter_IIR_NOTCH = IIR( notch_a, notch_b, order=2 ) 
    _filter_IIR_NOTCH_2 = IIR( notch_a_2, notch_b_2, order=2 ) 

    print("lpf_a: %s" % lpf_a)
    print("lpf_b: %s" % lpf_b)

    biq_b, biq_a = calculate_biquad_lpf_coeff( LPF_FC_1, LPF_Z_1, SAMPLE_FREQ )
    w_lpf_2, h_lpf_2 = freqz( biq_b, biq_a, 4096 )

    biq_b_hpf, biq_a_hpf = calculate_biquad_hpf_coeff( HPF_FC_1, HPF_Z_1, SAMPLE_FREQ )
    w_2, h_2 = freqz( biq_b_hpf, biq_a_hpf, 4096 )

    print("biq_a: %s" % biq_a)
    print("biq_b: %s" % biq_b)

    print("a_b_lpf: %s" % a_b_lpf)
    print("b_b_lpf: %s" % b_b_lpf)
    print("notch_a: %s" % notch_a)
    print("notch_b: %s" % notch_b)

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
    _fg_ac_noise = FunctionGenerator( 10.0, 0.2, 0, 0, "sine" )
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
    
    #ax[1].plot(w_b, 20 * np.log10(abs(h_b)), 'r', label="butterworth" )
    #ax[1].plot(w_c, 20 * np.log10(abs(h_c)), 'b', label="chebysev" )
    
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
