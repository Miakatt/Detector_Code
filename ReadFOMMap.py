# READ FOM MAP DATA
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from numpy import unravel_index

def Open(arg):
    fullfilename = sys.argv[arg]
    filename, file_ext = os.path.splitext(sys.argv[arg])
    print(filename, file_ext)
    return fullfilename, file_ext


def Read(fullfilename):
	print fullfilename
	data = list(csv.reader(open(fullfilename)))
	print data
	dataArray  = np.array(data)
	floatArray = dataArray.astype(np.float)

	return floatArray

def PlotStuff(floatArray, W1Start, W1EndFirst, W1EndLast, W1StepSize, W2Start, W2EndFirst, W2EndLast, W2StepSize):
	
	plt.figure(facecolor='white', figsize=(10,8))
	plt.imshow(floatArray, interpolation='nearest') #, interpolation='none')
	#print 'shape : ', XSize, YSize
	plt.xlabel('Window 2 Integral (Starting at %1.2f us)' %  (W2Start*SampleRate*1e6))
	plt.ylabel('Window 1 Integral (Starting at %1.2f us)' %  (W1Start*SampleRate*1e6 ))

	# Integrating Sample Numbers
	XAxis = np.arange(W2Start+W2EndFirst, W2Start+W2EndLast+1, W2StepSize)
	YAxis = np.arange(W1Start+W1EndFirst, W1Start+W1EndLast+1, W1StepSize)	
	print YAxis

	ax = plt.gca()
	ax.set_aspect('auto')
	majorFormatter = FormatStrFormatter('%1.2f')
	majorXLocator  = MultipleLocator(2)
	minorXLocator  = MultipleLocator(1)
	majorYLocator  = MultipleLocator(1)
	minorYLocator  = MultipleLocator(1)

	ax.xaxis.set_major_locator(majorXLocator)
	ax.xaxis.set_minor_locator(minorXLocator)

	ax.yaxis.set_major_locator(majorYLocator)
	ax.yaxis.set_minor_locator(minorYLocator)

	# Format X Tick Marks
	XLabel = []
	x_axis_labels = ax.get_xticks().tolist()
	x_axis_labels.pop(-1)
	for i, X in enumerate(x_axis_labels):
		XLabel.append('{0:.2f}'.format(XAxis[X]*(SampleRate*1e6)))
	print 'XLabel  ', XLabel
	ax.set_xticklabels(XLabel)	
	

	# Format Y Tick Marks
	YLabel = []
	y_axis_labels = ax.get_yticks().tolist()
	

	print 'y_axis_labels ' , y_axis_labels
	for i, Y in enumerate(YAxis):
		YLabel.append('{0:.2f}'.format(YAxis[i-1]*(SampleRate*1e6)))
	print 'YLabel  ', YLabel	
	ax.set_yticklabels(YLabel)	

	print 'Maximum FOM Value : ' ,  np.amax(floatArray)
	MaxIndex =  np.unravel_index(floatArray.argmax(), floatArray.shape )
	print 'Max Index  :  ', MaxIndex	
	print 'W1Start :', W1Start*SampleRate*1e6
	print 'W1End for best FOM :' , YAxis[MaxIndex[0]]*SampleRate*1e6

	print 'W2Start :', W2Start*SampleRate*1e6
	print 'W2End for best FOM :' , XAxis[MaxIndex[1]]*SampleRate*1e6

	plt.colorbar()
	plt.draw()
	plt.pause(0.01)
	Title = raw_input('Enter Title   ')
	plt.title(Title)
	SaveTitle = Title.replace(" ", "_")
	plt.savefig(SaveTitle+'.pdf')
	plt.show()

#Sample Rate = 8ns/Sample
SampleRate = 8e-9

fullfilename, ext = Open(1)
floatArray = Read(fullfilename)
################### FOR CLLBC (and possibly CLLB) ###################
#(floatArray, W1Start, W1EndFirst, W1EndLast, W1StepSize, W2Start, W2EndFirst, W2EndLast, W2StepSize)
# From settings used in PSDOptimize.py
# 20 C Course Scan
#PlotStuff(floatArray, 20.0, 20.0, 240.0, 40,  150.0, 400.0, 1800.0, 40.0)
# 20C Scan
#PlotStuff(floatArray, 20,10,220,20, 100, 200, 600, 20)
# 40C and 50C
#PlotStuff(floatArray, 20.0, 10.0, 120.0, 10,  100.0, 80.0, 460.0, 20.0)
# -10C
#PlotStuff(floatArray, 5.0, 5.0, 80.0, 10,  30.0, 20.0, 420.0, 20.0)
# -10C 2nd Scan
#PlotStuff(floatArray, 5, 5, 80, 10, 30, 20, 400, 20)
# 0C
#PlotStuff(floatArray, 20, 20, 120, 10, 100, 80, 360, 20)
#-20C
#PlotStuff(floatArray, 5, 5, 80, 10, 30, 20, 400, 20)



######################### FOR CLLB Standard Method ############################
 #  Use ReadFOMMapDI.py for the delayed integration data.
#-20C
#PlotStuff(floatArray, 5, 5, 55, 5, 60, 20, 620, 20)
#-10C
PlotStuff(floatArray, 20, 20, 60, 5, 100, 100, 400, 10)

######################### FOR CLLB DelayedInt Method ############################
#-20C
#PlotStuff(floatArray, 10, 40, 100, 200, 10, 200, 800, 10)








# FOR CLYC
# 20C
#PlotStuff(floatArray, -5., 20., 50., 2., 70., 20., 200., 20.)
# 50C
#PlotStuff(floatArray, -5., 18., 70., 2., 40., 20., 260., 20.)
# 40C
#PlotStuff(floatArray, -5., 20., 50., 2., 60., 20., 160., 20.)
# 40C again. Larger scan - same settings -10C 
#PlotStuff(floatArray, -5., 16., 70., 2., 40., 20., 260., 20.)
# 0C
#PlotStuff(floatArray, -5., 20., 70., 2., 40., 20., 200., 20.)
#-10C
#PlotStuff(floatArray, -5., 10., 64., 2., 40., 20., 260., 20.)