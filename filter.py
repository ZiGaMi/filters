import sys
import matplotlib.pyplot as plt
import numpy as np


# Time window
TIME_WINDOW = 20.0 # seconds

# Number of samples in time window
SAMPLE_NUM = 1000



class RC_LPF:

    def __init__(self, fc, dt, order=1, init_val=0):

        # Calculate filter coefficient
        self.k = dt / ( dt + 1/fc )

        # Current value
        self.val = init_val

    def update(self, x):
        self.val = ( self.val + self.k * ( x - self.val ))
        return self.val


def generate_sine(time, freq, amp, off, phase):
    return (( amp * np.sin((2*np.pi*freq*time) + phase )) + off )

def generate_rect(time, freq, amp, off, phase):
    _carier = generate_sine(time, freq, amp, off, phase)
    _sig = 0

    if ( _carier > 0 ):
        _sig = amp + off
    else:
        _sin = off

    return _sig
        




if __name__ == "__main__":

    # Time array
    _time, _dt = np.linspace( 0.0, TIME_WINDOW, num=SAMPLE_NUM, retstep=True )

    # Filter object
    _filter_1Hz = RC_LPF( fc=1.0, dt=_dt, order=1, init_val=0)
    _filter_10Hz = RC_LPF( fc=10.0, dt=_dt, order=1, init_val=1)

    # Filter input/output
    _x = []
    _sin_x = []
    _rect_x = []

    # Generate inputs
    for n in range(SAMPLE_NUM):
        _sin_x.append( generate_sine( _time[n], .10, 1.0, 0.0, 0.0 ))  
        _rect_x.append( generate_rect( _time[n], .10, 1.0, 0.0, 0.0 ))  

    _y = []
    _y_1 = []

    # Apply filter
    for n in range(SAMPLE_NUM):
        _y.append( _filter_1Hz.update( _rect_x[n] ) )
        _y_1.append( _filter_10Hz.update( _rect_x[n] ) )


    fig, (ax_1, ax_2) = plt.subplots(2, 1)
   # ax_1.plot( _time, _x, ".-b" )
    ax_1.plot( _time, _y, ".-r" )
    ax_1.plot( _time, _y_1, ".-g" )
    ax_1.grid()

    #ax_2.plot( _time, _y, ".-r" )
    #ax_2.plot( _time, _y_1, ".-g" )
    ax_2.plot( _time, _sin_x, ".-r" )
    ax_2.plot( _time, _rect_x, ".-g" )
    ax_2.grid()


    plt.show()
    
