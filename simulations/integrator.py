# ===============================================================================
# @file:    integrator.py
# @note:    This script is evaluation of discrete integrator
# @author:  Ziga Miklosic
# @date:    20.02.2021
# @brief:   Evaluation of digital integrator
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

## ****** USER CONFIGURATIONS ******

## Sample frequency
#   Sample frequency of real system   
#
# Unit: Hz
SAMPLE_FREQ = 50.0

# Ideal sample frequency
#   As a reference to sample rate constrained embedded system
#
# Unit: Hz
IDEAL_SAMPLE_FREQ = 20000.0

## Time window
#
# Unit: second
TIME_WINDOW = 3.14


## ****** END OF USER CONFIGURATIONS ******

## Number of samples in time window
SAMPLE_NUM = int(( IDEAL_SAMPLE_FREQ * TIME_WINDOW ) + 1.0 )


# ===============================================================================
#       CLASSES
# ===============================================================================

class SimpleIntegrator:

    def __init__(self, fs):
        self.y = 0
        self.dt = 1 / fs

    def update(self, x):
        self.y += x * self.dt
        return self.y

    def get(self):
        return self.y


class TrapezoidIntegrator:

    def __init__(self, fs):
        self.y_prev = 0
        self.x_prev = 0
        self.y = 0
        self.dt = 1 / fs

    def update(self, x):
        self.y += ( x + self.x_prev ) * self.dt / 2.0
        self.x_prev = x
        return self.y

    def get(self):
        return self.y


class KeplerIntegrator:
    
    def __init__(self, fs):
        self.y_prev = 0
        self.x_prev_1 = 0
        self.x_prev_2 = 0
        self.y = 0
        self.dt = 1 / fs

    def update(self, x):
        self.y += ( x + 4*self.x_prev_1 + self.x_prev_2 ) * self.dt/6
        self.x_prev_2 = self.x_prev_1
        self.x_prev_1 = x
        return self.y 
    
    def get(self):
        return self.y


# ===============================================================================
#       MAIN ENTRY
# ===============================================================================
if __name__ == "__main__":

    # Time array
    _time, _dt = np.linspace( 0.0, TIME_WINDOW, num=SAMPLE_NUM, retstep=True )


    # Create integrators
    _simple_int = SimpleIntegrator( fs=SAMPLE_FREQ )
    _trap_int = TrapezoidIntegrator( fs=SAMPLE_FREQ )
    _kepler_int = KeplerIntegrator( fs=SAMPLE_FREQ )    

    _ideal_simple_int = SimpleIntegrator( fs=IDEAL_SAMPLE_FREQ )
    _ideal_trap_int = TrapezoidIntegrator( fs=IDEAL_SAMPLE_FREQ )
    _ideal_kepler_int = KeplerIntegrator( fs=IDEAL_SAMPLE_FREQ )    

    # Integrals to plot
    _d_simple_int = [0]
    _d_trap_int = [0]
    _d_kepler_int = [0]

    _i_simple_int = []
    _i_trap_int = []
    _i_kepler_int = []

    # Signal to integrate
    _x_d = [0]

    K = -1
    N = 1
    input = 0

    # Down sample
    _downsamp_cnt = 0
    _downsamp_samp = [0]
    _d_time = [0]
 
    # Apply filter
    for n in range(SAMPLE_NUM):
        
        input = np.sin( 2*_time[n] )

        # Down sample to SAMPLE_FREQ
        if _downsamp_cnt >= (( 1 / ( _dt * SAMPLE_FREQ )) - 1 ):
            _downsamp_cnt = 0

            # Utils
            _downsamp_samp.append(0)
            _d_time.append( _time[n])    

            # Perform integration
            _simple_int.update( input )
            _trap_int.update( input )
            _kepler_int.update( input )

            # Store data            
            _d_simple_int.append( _simple_int.get())
            _d_trap_int.append( _trap_int.get())
            _d_kepler_int.append( _kepler_int.get())

            # Discrete input
            _x_d.append( input )

        else:
            _downsamp_cnt += 1

        # ============================================
        # IDEAL SAMPLING FREQUENCY
        # ============================================

        _ideal_simple_int.update( input )
        _ideal_trap_int.update( input )
        _ideal_kepler_int.update( input )

        _i_simple_int.append( _ideal_simple_int.get())
        _i_trap_int.append( _ideal_trap_int.get())
        _i_kepler_int.append( _ideal_kepler_int.get())
    

    print("Simple integral: %s" % _simple_int.get() )
    print("Trap integral: %s" % _trap_int.get() )
    print("Kepler integral: %s" % _kepler_int.get() )


    # =============================================================================================
    ## PLOT CONFIGURATIONS
    # =============================================================================================
    plt.style.use(['dark_background'])
    PLOT_MAIN_TITLE_SIZE    = 18
    PLOT_MAIN_TITLE         = "DISCRETE INTEGRATOR SIMULATIONS | fs: " + str(SAMPLE_FREQ) + "Hz" 
    PLOT_TITLE_SIZE         = 16
    PLOT_AXIS_LABEL_SIZE    = 12
    PLOT_ADJUST_LEFT        = 0.06
    PLOT_ADJUST_RIGHT       = 0.98
    PLOT_ADJUST_TOP         = 0.91
    PLOT_ADJUST_BOTTOM      = 0.05


    ## ==============================================================================================
    #     
    ## ==============================================================================================
    fig, ax = plt.subplots(1, 2, sharex=True)
    fig.suptitle( PLOT_MAIN_TITLE , fontsize=PLOT_MAIN_TITLE_SIZE )

    # Subplot 0
    ax[0].plot(_d_time, _x_d, "w")

    # real
    ax[1].plot(_d_time, _d_simple_int, ".y", label="simple")
    ax[1].plot(_d_time, _d_trap_int, ".b", label="trap")
    ax[1].plot(_d_time, _d_kepler_int, ".r", label="simp")

    # ideal
    #ax[1].plot(_time, _i_simple_int, "w", label="simple")
    #ax[1].plot(_time, _i_trap_int, "w", label="trap")
    ax[1].plot(_time, _i_kepler_int, "w", label="simp")
    



    #ax[0][0].set_title("Input acceleration & rotation", fontsize=PLOT_TITLE_SIZE)
    ax[0].grid(alpha=0.25)
    ax[1].grid(alpha=0.25)
    ax[1].legend(loc="upper right")
    #ax[0][0].set_ylabel('Acceleration [m/s^2],\nRotation [rad]', fontsize=PLOT_AXIS_LABEL_SIZE)


    plt.show()


# ===============================================================================
#       END OF FILE
# ===============================================================================