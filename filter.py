import sys
import matplotlib.pyplot as plt
import numpy as np


# Time window
TIME_WINDOW = 10.0 # seconds

# Number of samples in time window
SAMPLE_NUM = 10000

# Select input filter signal type
INPUT_SINE = 0
INPUT_RECT = 1

# Mux input signal
INPUT_SIGNAL_SELECTION = INPUT_SINE



def input_signal_mux(sel, in_1, in_2):
    if ( INPUT_SINE == sel ):
        return in_1
    elif ( INPUT_RECT == sel ):
        return in_2
    else:
        pass

def generate_sine(time, freq, amp, off, phase):
    return (( amp * np.sin((2*np.pi*freq*time) + phase )) + off )

def generate_rect(time, freq, amp, off, phase):
    _carier = generate_sine(time, freq, 1.0, 0.0, phase)
    _sig = 0

    if ( _carier > 0 ):
        _sig = amp + off
    else:
        _sig = off

    return _sig 

class RC_LPF:

    def __init__(self, fc, dt, order=1, init_val=0):

        # Calculate filter coefficient
        self.k = dt / ( dt + 1/fc )

        # Current value
        self.val = init_val

    def update(self, x):
        self.val = ( self.val + self.k * ( x - self.val ))
        return self.val


class RC_HPF:

    def __init__(self, fc, dt, order=1, init_val=0):

        # Calculate filter coefficient
        self.k = (1/fc) / ( dt + 1/fc )

        # Current value
        self.val = init_val

        # Prev input value
        self.x_1 = init_val

    def update(self, x):
        self.val = ( self.k * self.val + self.k * ( x - self.x_1 ))
        self.x_1 = x

        return self.val


        




if __name__ == "__main__":

    # Time array
    _time, _dt = np.linspace( 0.0, TIME_WINDOW, num=SAMPLE_NUM, retstep=True )

    # Filter object
    _filter_1Hz = RC_LPF( fc=1.0, dt=_dt, order=1, init_val=0)
    _filter_10Hz = RC_HPF( fc=10.0, dt=_dt, order=1, init_val=0)

    # Filter input/output
    _x = [ 0 ] * SAMPLE_NUM
    _y = []
    _y_1 = []

    # Generate inputs
    _sin_x = []
    _rect_x = []
    
    for n in range(SAMPLE_NUM):
        _sin_x.append( generate_sine( _time[n], 1.0, 1.0, 0.0, 0.0 ))  
        _rect_x.append( generate_rect( _time[n], 1.0, 1.0, 0.0, 0.0 ))  
 

    # Apply filter
    for n in range(SAMPLE_NUM):
        _x[n] = input_signal_mux( INPUT_SIGNAL_SELECTION, _sin_x[n], _rect_x[n] )
        _y.append( _filter_1Hz.update( _x[n] ) )
        _y_1.append( _filter_10Hz.update( _x[n] ) )


    fig, (ax_1, ax_2) = plt.subplots(2, 1)
    ax_1.plot( _time, _x, "b" )
    ax_1.plot( _time, _y, "g" )
    #ax_1.plot( _time, _y_1, ".-" )
    ax_1.grid()

    ax_2.plot( _time, _x, "b" )
    ax_2.plot( _time, _y_1, "g" )
    #ax_2.plot( _time, _rect_x, ".-g" )
    ax_2.grid()


    plt.show()
    
