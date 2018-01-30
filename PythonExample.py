# Open and read the text files from the digitizer.
# Get the number of samples from the file (Given by Record Length in the txt file) 
# and plot the data, waveform-by-waveform.
import os
import sys
import matplotlib.pyplot as plt
import numpy as np



# -----------------------------------------------------------


def OpenFile(inputfile):
	# Open the file passed as the first argument in the command line
	filename = inputfile
	f = open (filename, "r")
	# Get the Record Length from the file. 
	FirstLine = f.readline()
	print FirstLine
	
	return f, cols, rows

# -----------------------------------------------------------




OpenFile(1)