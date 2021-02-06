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

from rc_filter import RC_LPF, CR_HPF
from fir_filter import FIR
from iir_filter import IIR

# ===============================================================================
#       CONSTANTS
# ===============================================================================

# ===============================================================================
#       FUNCTIONS
# ===============================================================================

# ===============================================================================
#       CLASSES
# ===============================================================================

# ===============================================================================
#       MAIN ENTRY
# ===============================================================================

# Arg parser
parser = argparse.ArgumentParser( 	description="This script is for evaluation of filtering raw data from CSV file", 
									epilog="Enjoy the program!\n")

# Add arguments
parser.add_argument("-f", "--file",			help="file to analyze", 						type=str, 	required=True)
parser.add_argument("-fs", "--sample_freq",	help="sample time in Hz", 						type=float, required=True)
parser.add_argument("-fc", "--cutoff_freq",	help="cutoff freq in Hz", 						type=float, required=False)
parser.add_argument("-fl", "--filter", 		help="type of filter. Options: 'RC', 'CR'", 	type=str, 	required=True)
parser.add_argument("-or", "--order", 		help="order of filter", 						type=int, 	required=False)

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

	# Sample frequecy & time
	_fs = args["sample_freq"]
	_ts = 1 / _fs

	# Signal
	_signal = []
	_signal_filter = []

	# Filter
	if args["filter"] == "RC":
		_filter = RC_LPF( args["cutoff_freq"], _ts, args["order"] )
	elif args["filter"] == "CR":
		_filter = CR_HPF( args["cutoff_freq"], _ts, args["order"] )
	elif args["filter"] == "FIR":
		_filter = FIR( 25, [-0.005392788920679466,
-0.0029050919204011926,
-0.001258909518015039,
0.0031062514135009396,
0.010866599224457303,
0.02226527409028689,
0.036982494521830354,
0.05407161051837529,
0.07200917783384833,
0.0889018940568695,
0.10277548127558198,
0.11189425867893746,
0.11507335279221934,
0.11189425867893746,
0.10277548127558198,
0.0889018940568695,
0.07200917783384833,
0.05407161051837529,
0.036982494521830354,
0.02226527409028689,
0.010866599224457303,
0.0031062514135009396,
-0.001258909518015039,
-0.0029050919204011926,
-0.005392788920679466 ] )
	elif args["filter"] == "IIR":
		_filter = IIR( [1.0, -1.04377111, 0.27236453], [0.05714836, 0.11429671, 0.05714836], order=2 )
	else:
		raise AssertionError

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
                	
					_time.append( _time[-1] + _ts )
					_signal.append( float(row[0]) + 1.0 )

					# Apply filter
					_signal_filter.append( _filter.update( _signal[-1] ) )

				else:
					_time.append( 0.0 )
					_signal.append( 0.0 )
					_signal_filter.append( 0.0 )

				# Sample number
				_n = _n + 1



	## ================================================================================
	## FFT
    ## ================================================================================	

	_signal_fft = rfft( _signal )
	_signal_filter_fft = rfft( _signal_filter )
	_freq = rfftfreq( _n, _ts )

	## ================================================================================
	## PLOT DATA
    ## ================================================================================
	fig, (ax_1, ax_2) = plt.subplots(2, 1)

	fig.suptitle("Import file: " + str(_file) + "\nFilter type: " + str(args["filter"]) + "\nfs: " + str(_fs) + "Hz, fc: " + str(args["cutoff_freq"]) + \
				"Hz, order: " + str(args["order"]), fontsize=20)

	ax_1.title.set_text("Time domain")
	ax_1.plot( _time, _signal, 			"g", label="raw" )
	ax_1.plot( _time, _signal_filter, 	"b", label="filtered" )
	ax_1.grid()
	ax_1.legend(loc="upper right")
	ax_1.set_xlabel("Time [s]")
	ax_1.set_ylabel("Acceleration [m/s^2]")

	ax_2.title.set_text("Frequency domain")
	ax_2.plot( _freq, np.abs(_signal_fft), 			"g", label="raw" )
	ax_2.plot( _freq, np.abs(_signal_filter_fft), 	"b", label="filtered" )
	ax_2.grid()
	ax_2.legend(loc="upper right")
	ax_2.set_xlabel("Freuqncy [Hz]")
	ax_2.set_ylabel("Power")

	plt.show()


# ===============================================================================
#       END OF FILE
# ===============================================================================
