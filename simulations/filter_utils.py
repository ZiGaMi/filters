# ===============================================================================
# @file:    filter_utils.py
# @note:    This script is has various utility classes/function
# @author:  Ziga Miklosic
# @date:    10.01.2021
# @brief:   Classes/Functions for general usage to evaluate filter
# ===============================================================================

# ===============================================================================
#       IMPORTS  
# ===============================================================================
import sys
import numpy as np

# ===============================================================================
#       CONSTANTS
# ===============================================================================

# ===============================================================================
#       CLASSES
# ===============================================================================

# Signal fucntion generator
class FunctionGenerator:

    # ===============================================================================
    # @brief: Initialize function generator object
    #
    # @param[in]:    freq    - Frequency of signal  
    # @param[in]:    amp     - Amplitude of signal
    # @param[in]:    off     - DC offset of signal
    # @param[in]:    phase   - Phase of signal
    # @return:       Generated signal
    # ===============================================================================
    def __init__(self, freq, amp, off, phase, kind):

        self.freq = freq
        self.amp = amp
        self.off = off
        self.phase = phase
        self.kind = kind

    # ===============================================================================
    # @brief: Generate signal selected at init time
    #
    # @param[in]:    time    - Linear time  
    # @return:       Generated signal
    # ===============================================================================
    def generate(self, time):
        
        _sig = 0

        if self.kind == "sine":
            _sig =  self.__generate_sine(time, self.freq, self.amp, self.off, self.phase)
        elif self.kind == "rect":
            _sig = self.__generate_rect(time, self.freq, self.amp, self.off, self.phase)
        else:
            raise AssertionError

        return _sig

    # ===============================================================================
    # @brief: Generate sine signal
    #
    # @param[in]:    time    - Linear time  
    # @param[in]:    amp     - Amplitude of sine
    # @param[in]:    off     - DC offset of sine
    # @param[in]:    phase   - Phase of sine
    # @return:       Generated signal
    # ===============================================================================
    def __generate_sine(self, time, freq, amp, off, phase):
        _sig = (( amp * np.sin((2*np.pi*freq*time) + phase )) + off )
        return _sig

    # ===============================================================================
    # @brief: Generate rectangle signal
    #
    # @param[in]:    time    - Linear time  
    # @param[in]:    amp     - Amplitude of rectange
    # @param[in]:    off     - DC offset of rectangle
    # @param[in]:    phase   - Phase of rectangle
    # @return:       Generated signal
    # ===============================================================================
    def __generate_rect(self, time, freq, amp, off, phase):
        _carier = self.__generate_sine(time, freq, 1.0, 0.0, phase)
        _sig = 0

        if ( _carier > 0 ):
            _sig = amp + off
        else:
            _sig = off

        return _sig 


# Signal multiplexor
class SignalMux:

    # Signal mux options
    MUX_CTRL_SINE = 0   # Sine
    MUX_CTRL_RECT = 1   # Rectange

    # ===============================================================================
    # @brief: Initialize multiplexor
    #
    # @param[in]:    size - Size of input signals
    # @return:       void
    # ===============================================================================
    def __init__(self, size=1):
        self.size = size

    # ===============================================================================
    # @brief: Route output of mulitplexor based on control lines
    #
    # @param[in]:    ctrl       - Control line  
    # @param[in]:    sig_in     - Array of signals (expected of size specified at init)
    # @return:       out        - Output of mux
    # ===============================================================================
    def out(self, ctrl, sig_in):
        
        _out = 0
        
        for n in range( self.size ):
            if n == ctrl:
                _out = sig_in[n]
                break

        return _out


# ===============================================================================
#       END OF FILE
# ===============================================================================
