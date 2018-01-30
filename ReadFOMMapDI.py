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
	dataArray  = np.array(data)
	floatArray = dataArray.astype(np.float)

	return floatArray

def PlotStuff(floatArray, W2Start , W2End):
    plt.figure(facecolor='white', figsize=(8,8))
    plt.imshow(floatArray, interpolation='nearest') #, interpolation='none')
    #print 'shape : ', XSize, YSize
    plt.xlabel('Window 2 Integral (us) ' )
    plt.ylabel('Window 2 Start delay after trigger (us) ' )

    
    XAxis = W2End
    YAxis = W2Start
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
    ax.set_xticklabels(XLabel)	
    
    # Format Y Tick Marks
    YLabel = []
    y_axis_labels = ax.get_yticks().tolist()
    print 'YAxis Labels ' , y_axis_labels
    y_axis_labels.pop(-1)
    for i, Y in enumerate(y_axis_labels):
        YLabel.append('{0:.2f}'.format(YAxis[Y]*(SampleRate*1e6)))
        ax.set_yticklabels(YLabel)	
    
    print 'Maximum FOM Value : ' ,  np.amax(floatArray)
    MaxIndex =  np.unravel_index(floatArray.argmax(), floatArray.shape )
    print 'Max Index  :  ', MaxIndex	
    
    #print 'W1End for best FOM :' , YAxis[MaxIndex[0]]*SampleRate*1e6
    
    print 'W2Start :', W2Start*SampleRate*1e6
    #print 'W2End for best FOM :' , XAxis[MaxIndex[1]]*SampleRate*1e6
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

# For CLLB
#20C
#PlotStuff(floatArray, 10, 40, 100, 200, 10, 200, 800, 10)
# 0C

# W2Start = np.arange(100, 200, 10)
# W2End   = np.arange(200, 400, 10)
# For all temperatures except 50C
#W2Start = np.arange(100, 200, 10)
#W2End   = np.arange(200, 400, 10)
W2Start = np.arange(100, 181, 10)
W2End   = np.arange(200, 601, 40)


# PlotStuff(floatArray, W2Start[0], W2Start[-1], W2Stepsize, W2StartEnd[0], W2End[-1], W2EndStepSize)
PlotStuff(floatArray, W2Start, W2End )

