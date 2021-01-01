# ===============================================================================
# @file:    filter_csv.py
# @note:    Filter signals from CSV file
# @author:  Ziga Miklosic
# @date:    01.01.2021
# @brief:   This scripts read signal from CSV file and applies filter to it. 
# ===============================================================================

# ===============================================================================
#       IMPORTS  
# ===============================================================================
import sys
import matplotlib.pyplot as plt
import csv
import argparse
from scipy.fft import rfft, rfftfreq
import numpy as np

from rc_filter import RC_LPF

# ===============================================================================
#       CONSTANTS
# ===============================================================================

# Sample time of input file
SAMPLE_TIME = 0.01


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



# ===============================================================================
#       CLASSES
# ===============================================================================



# ===============================================================================
#       MAIN ENTRY
# ===============================================================================



# Arg parser
parser = argparse.ArgumentParser()

# File to plot
parser.add_argument("-f", "--file", help="File for analysis")

# Get args
args = parser.parse_args()

# Convert namespace to dict
args = vars(args)


if __name__ == "__main__":

	# Data files
	csv_file = 0
	
	# Time
	_time = []

	# Number of samples
	_n = 0

	# Signal
	_signal = []
	_signal_filter = []
	_signal_filter_1 = []

	# RC LPF
	_rc_lpf = RC_LPF( 40.0, SAMPLE_TIME, 2 )
	_rc_lpf_1 = RC_LPF( 50.0, SAMPLE_TIME, 3 )
	
	# Check for file
	if args["file"] != None:
		
		# Get file name
		_file = args["file"]
		
		# Print status
		print("Parsing... (%s)" % _file)

		# Open file for reading
		with open(_file, "r") as csvfile:
			
			# Read row
			spamreader = csv.reader(csvfile, delimiter=",")
		
			for idx, row in enumerate(spamreader):
				
				# Ignore first row
				if idx > 0:
                	
					_time.append( _time[-1] + SAMPLE_TIME )
					_signal.append( float(row[0]) )

					_signal_filter.append( _rc_lpf.update( _signal[-1] ) )
					_signal_filter_1.append( _rc_lpf_1.update( _signal[-1] ) )

				else:
					_time.append( 0.0 )
					_signal.append( 0.0 )

					_signal_filter.append( 0.0 )
					_signal_filter_1.append( 0.0 )

				# Sample number
				_n = _n + 1



	## ================================================================================
	## FFT
    ## ================================================================================	

	_signal_fft = rfft( _signal )
	_signal_filter_fft = rfft( _signal_filter )
	_signal_filter_1_fft = rfft( _signal_filter_1 )
	_freq = rfftfreq( _n, SAMPLE_TIME )

	## ================================================================================
	## PLOT DATA
    ## ================================================================================
	fig, (ax_1, ax_2) = plt.subplots(2, 1)

	fig.suptitle("Import file: " + str(_file), fontsize=20)

	ax_1.title.set_text("Time domain")
	ax_1.plot( _time, _signal, 			"g", label="raw" )
	ax_1.plot( _time, _signal_filter, 	"b", label="lpf" )
	ax_1.plot( _time, _signal_filter_1,	"r", label="lpf1" )
	ax_1.grid()
	ax_1.legend(loc="upper right")
	ax_1.set_xlabel("Time [s]")
	ax_1.set_ylabel("Acceleration [m/s^2]")

	ax_2.title.set_text("Frequency domain")
	ax_2.plot( _freq, np.abs(_signal_fft), 			"g", label="raw" )
	ax_2.plot( _freq, np.abs(_signal_filter_fft), 	"b", label="lpf" )
	ax_2.plot( _freq, np.abs(_signal_filter_1_fft), "r", label="lpf1" )
	ax_2.grid()
	ax_2.legend(loc="upper right")
	ax_2.set_xlabel("Freuqncy [Hz]")
	ax_2.set_ylabel("Power")

	plt.show()


# ===============================================================================
#       END OF FILE
# ===============================================================================
