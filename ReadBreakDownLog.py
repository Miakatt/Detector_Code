# ReadBreakDownLog.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.optimize import curve_fit
import csv
import sys
import os
from matplotlib import gridspec
import math


def ReadCSV(arg):
	filename = sys.argv[arg]
	print filename

	with open(filename, 'r') as csvfile:
		read = csv.reader(csvfile, delimiter=',')
		for row in read:
			Temperature.append(float(row[0]))
			BD.append(float(row[1]))
			
	return Temperature, BD



def PlotData(Temperature, BD):
	plt.figure(figsize=(8,6), facecolor='white')
	plt.plot(Temperature, BD, 'ro')
	ax = plt.gca()
	ax.ticklabel_format(useOffset=False, style='plain')
	plt.draw()
	plt.pause(0.01)
#--------------------------------------------------
Temperature = []
BD = []


Temperature, BD = ReadCSV(1)
PlotData(Temperature, BD)
plt.show()



