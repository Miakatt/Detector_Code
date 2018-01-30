# Make FOM Map of mixed field data.

import numpy as np 
import sys
import time
import struct
import codecs
import matplotlib.pyplot as plt
import matplotlib as mpl
from operator import truediv
from scipy.optimize import curve_fit
from scipy.ndimage.interpolation import zoom
import scipy.signal as ss

def OpenFile(inputfile):
	filename = sys.argv[inputfile]
	print filename
	f = open (filename, "rb")
	binary_datapoint = f.read(4)
	cols = struct.unpack("i", binary_datapoint)[0]
	print "Columns ", cols

	binary_datapoint = f.read(4)
	rows = struct.unpack("i", binary_datapoint)[0]
	print "Rows ", rows
	
	return f, cols, rows

# -----------------------------------------------------------


def ReadBinaryFile(f, iterations, cols, rows):
	for it in range(iterations):
		chan = 0
		nData = 0
		del	WaveformA[:]	

		# Read 8 bytes = Double
		binary_datapoint = f.read(8)
		# 2400 = 600 samples/waveform x 4 channels in this case
		while binary_datapoint != b"" and nData<rows*cols:
			num_datapoint = struct.unpack("d", binary_datapoint)[0]
			# Fill datapoints for ChA
			if (nData%4 == 0): 
				WaveformA.append(num_datapoint)

			binary_datapoint = f.read(8)
			nData = nData + 1
			
		# Stack Waveforms in a nxm matrix	
		WaveformMatrix = StackWaveforms(WaveformA[0:cols])

	
# -----------------------------------------------------------

def StackWaveforms(WaveformA):
	WaveformMatrix.append(WaveformA)
	
	return WaveformMatrix

# -----------------------------------------------------------

def GetTriggerPosition(Matrix, Threshold):
	# Find the list entry that crosses the trigger threshold.
	TriggerIndexList = np.nanargmax(Matrix>Threshold, axis=1)
#	print TriggerIndexList
	return TriggerIndexList
# -----------------------------------------------------------

def BaselineSubtraction(Matrix, BaselineSamples):
	for b in range(Matrix.shape[0]):
		AvgBaseline = sum(Matrix[b,1:BaselineSamples] / BaselineSamples)
		#print "Average Baseline " , AvgBaseline
		Matrix[b,:] = Matrix[b,:]-AvgBaseline

	return Matrix
# -----------------------------------------------------------

def GetFullIntegral(Matrix):
	# Get Full Integrals for each waveform
	FullIntegral = Matrix.sum(1)
	#print "Full Int ", FullIntegral
	return FullIntegral

# -----------------------------------------------------------

def GetW1Integral(Matrix, TriggerIndexList, W1Start, W1End):
	W1Integral = []
	W1Start = TriggerIndexList+W1Start
	W1End = TriggerIndexList+W1End
	for S in range(len(W1Start)):
		# sum the values between W1Start and W1End and reduce to 3 decimal places
		W1Integral.append('{0:.3f}'.format(sum(Matrix[S, W1Start[S]:W1End[S]])))

	return W1Integral
# -----------------------------------------------------------


def GetW2Integral(Matrix, TriggerIndexList, W2Start, W2End):
	W2Integral = []
	W2Start = TriggerIndexList+W2Start
	W2End = TriggerIndexList+W2End
	for S in range(len(W2Start)):
		# sum the values between W1Start and W1End and reduce to 3 decimal places
		W2Integral.append('{0:.3f}'.format(sum(Matrix[S, W2Start[S]:W2End[S]])))

	return W2Integral
# -----------------------------------------------------------

def GetSignalPeak(Matrix):
	MaxSignal = []
	print len(Matrix)
	for M in range(len(Matrix)):
		MaxSignal.append(max(Matrix[M,:]))
	return MaxSignal


# -----------------------------------------------------------


def GetPSDFactor(SignalPeak, W1Integral, W2Integral):
	if (CLYC) or (Test):
		PSDFactor = map(truediv, W1Integral+W2Integral,W2Integral)
	elif (CLLB) or (CLLBC) or (Test):
		PSDFactor = map(truediv, W1Integral, W2Integral)
		# Try doing maximum of waveform/W2 Integral
		#PSDFactor = SignalPeak/W2Integral


	return PSDFactor

# -----------------------------------------------------------

def RemoveZeroValues(SignalPeak, W1Integral, W2Integral, FullIntegral):
# Remove any zero values in W2Integral, and corresponding values in W2Integral and FullIntegral
	SignalPeak    = np.array(SignalPeak)
	W1Integral    = np.array(W1Integral)
	W2Integral    = np.array(W2Integral)
	FullIntegral  = np.array(FullIntegral)
	#print len(SignalPeak), len(W1Integral), len(W2Integral), len(FullIntegral)

	NonZeroValues = W2Integral.nonzero()
	
	SignalPeak    = SignalPeak[NonZeroValues]
	W1Integral    = W1Integral[NonZeroValues]
	W2Integral    = W2Integral[NonZeroValues]	
	FullIntegral  = FullIntegral[NonZeroValues]
	
	return SignalPeak, W1Integral, W2Integral, FullIntegral

# -----------------------------------------------------------

def PlotWaveforms(Waveform, TriggerIndexList, W1Start, W1End, W2Start, W2End):
		
	for wf in range(Waveform.shape[0]):
		plt.plot(Waveform[wf],'r-')
		plt.axvline(TriggerIndexList[wf]+W1Start, color='r', linestyle='--')
		plt.axvline(TriggerIndexList[wf]+W1End  , color='b', linestyle='--')
		plt.axvline(TriggerIndexList[wf]+W2Start, color='g', linestyle='--')
		plt.axvline(TriggerIndexList[wf]+W2End  , color='c', linestyle='--')
		plt.axhline(Threshold, color='k', linestyle='--')
		plt.ylim(-0.01, 0.1)
		plt.hold(False)
		plt.draw()
		plt.pause(0.01)
		plt.cla()

# -----------------------------------------------------------

def MakeScatterPlot(figScat, PSDFactor, FullIntegral, nXBins, nYBins, XMin,XMax, YMin, YMax, First):
	# Create scatter plot of PSDFactor vs Full Int
	# and get user to set Y Projection limits.
	
	plt.figure(figScat.number)
	h = plt.hist2d(FullIntegral, PSDFactor, bins=[nXBins,nYBins], range=[[XMin,XMax],[YMin, YMax]])
	if (First):
		plt.colorbar()
#	plt.xlim(0,15)
	plt.ylim(YMin,YMax)
	plt.draw()
	plt.pause(0.01)
	return np.array(h[0])
# -----------------------------------------------------------

def SetProjectionLimits():
	print 'Please click lower limit: '
	XLow  = plt.ginput(1, timeout=0)[0]
	#print XLow[0]
	LowLine  = plt.axvline(XLow[0],  color='r', linestyle='--')
	plt.draw()
	plt.pause(0.01)

	print 'Please click upper limit: '
	XHigh = plt.ginput(1, timeout=0)[0]
	#print XHigh[0]
	HighLine = plt.axvline(XHigh[0], color='r', linestyle='--')
	plt.draw()
	plt.pause(0.01)

	return XLow[0], XHigh[0]

# -----------------------------------------------------------


def smooth(x,window_len=10,window='hanning'):
  
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)]

# -----------------------------------------------------------

def ProjectY(scatterhist, LowLimits, HighLimits, nXBins, XMax, XMin, GammaRow, NeutronRow):
	#Get Bin number from the projection limits.
	Low  = int(LowLimits* nXBins/(XMax-XMin))
	High = int(HighLimits* nXBins/(XMax-XMin))

	#print 'Low ' , Low 
	#print 'High ', High
	# Make Y Projection of data within the set limits
	ProjY = sum(scatterhist[Low:High,:])


	GammaPeak   = np.argmax(scatterhist[GammaRow,10:])
	print max(scatterhist[GammaRow])
	NeutronPeak = np.argmax(scatterhist[NeutronRow,10:])
	plt.plot(ProjY)
	plt.draw()
	plt.pause(1)
	print 'Gamma Peak:  ', GammaPeak
	print 'NeutronPeak:  ', NeutronPeak
	return ProjY, GammaPeak, NeutronPeak

# -----------------------------------------------------------

def gaus(x, *p0):
    amp, mean, sigma = p0
    return amp*np.exp(-(x-mean)**2/(2*sigma**2))
# -----------------------------------------------------------

def FitGaussian(ProjY, GammaPeak, NeutronPeak):
#	ProjY = smooth(ProjY,3)
	plt.plot(ProjY,'b-')
	plt.xlim([GammaPeak-100, NeutronPeak+100])
	plt.draw()
	plt.pause(0.1)
	Sigma = 1.1
	
	# Choose, either manual fitting method ...
#	print 'Click on the peaks...'
	# Fit Gaussian to Peaks
	# Find two highest bins in projY
#	XY= plt.ginput(2, timeout=0)
#	GammaProjPeak = int(round(XY[0][0]))
#	Amp1 = XY[0][1]
#	NeutronProjPeak = int(round(XY[1][0]))
#	Amp2 = XY[1][1]
#	PeakGammaWidth = int(round(0.1*GammaProjPeak))
#	PeakNeutronWidth = int(round(0.1*NeutronProjPeak))

	GammaProjPeak = int(GammaPeak)
	NeutronProjPeak = int(NeutronPeak)
	Amp1 = ProjY[GammaPeak]
	Amp2 = ProjY[NeutronPeak]
	PeakGammaWidth = int(round(0.2*GammaProjPeak))
	PeakNeutronWidth = int(round(0.2*NeutronProjPeak))

	#print 'Peaks: ', GammaProjPeak, NeutronProjPeak
	#print 'Amplitudes: ', Amp1,  Amp2
	#print 'Widths: ' ,PeakGammaWidth, PeakNeutronWidth

	p0=[Amp1, GammaPeak, Sigma]
	fitrange = range(GammaPeak-PeakGammaWidth, GammaPeak+PeakGammaWidth)
	fitdata = ProjY[GammaPeak-PeakGammaWidth:GammaPeak+PeakGammaWidth]
	try:
		popt, pcov = curve_fit(gaus, fitrange, fitdata,  p0=p0)
		mean1  = popt[1]
		sigma1 = popt[2]

	except RuntimeError:
		return 0,0,1,1

	except ValueError:
		print "Unable to fit first peak. Try again"	
		print "Select 1st peak"
		GammaPeak = int(round(plt.ginput(1)[0][0]))
		fitrange = range(GammaPeak-PeakGammaWidth, GammaPeak+PeakGammaWidth)
		fitdata = ProjY[GammaPeak-PeakGammaWidth:GammaPeak+PeakGammaWidth]

		plt.cla()
		p0=[Amp1, GammaPeak, Sigma]
		popt, pcov = curve_fit(gaus, fitrange, fitdata,  p0=p0)
		mean1  = popt[1]
		sigma1 = popt[2]

	if (PlotHistos):

		plt.plot(fitrange, gaus(fitrange, *popt), 'red', linewidth=2, )  
	#	plt.xlim([10, 50])
  
		plt.draw()	
		plt.pause(0.01)
		plt.hold(True)

	# FIT SECOND PEAK

	p0=[Amp2, NeutronPeak, Sigma]
	fitrange = range(NeutronPeak-PeakNeutronWidth, NeutronPeak+PeakNeutronWidth) 
	fitdata  = ProjY[NeutronPeak-PeakNeutronWidth:NeutronPeak+PeakNeutronWidth]
	try:
		popt, pcov = curve_fit(gaus, fitrange, fitdata,  p0=p0)
		mean2  = popt[1]
		sigma2 = popt[2]  

	except RuntimeError:
		return 0,0,1,1
	except ValueError:
		print "Didn't manage to fit second peak. Try again"
		print "Select 2nd peak"
		NeutronPeak = int(round(plt.ginput(1)[0][0]))
		fitrange = range(NeutronPeak-PeakNeutronWidth, NeutronPeak+PeakNeutronWidth) 
		fitdata  = ProjY[NeutronPeak-PeakNeutronWidth:NeutronPeak+PeakNeutronWidth]

		plt.cla()
		p0=[Amp2, NeutronPeak, Sigma]
		popt, pcov = curve_fit(gaus, fitrange, fitdata,  p0=p0)
		mean2  = popt[1]
		sigma2 = popt[2]  

	if (PlotHistos):
		plt.plot(ProjY,'b-')
		plt.plot(fitrange, gaus(fitrange, *popt), 'green', linewidth=2, )  
	#	plt.xlim([10, 50])

		plt.draw()
		plt.pause(0.1)

	return mean1, mean2, sigma1, sigma2

# -----------------------------------------------------------

def GetPeakRow(scatterhist, nYBins, YMin, YMax):
	print 'Select peak of gamma data points.'
	GammaRow = int(round(plt.ginput(1)[0][0]) * (nXBins/(XMax-XMin)) )
	print 'Select peak of neutron data points.'
	NeutronRow = int(round(plt.ginput(1)[0][0]) * (nXBins/(XMax-XMin)) )
	print 'Gamma Row:  ', GammaRow 
	print 'Neutron Row: ', NeutronRow
	return GammaRow, NeutronRow


# =======================================================================	
# =======================================================================	
# 								Start program here	
# =======================================================================
# =======================================================================	
W1Axis = []
W2Axis = []
CLYC  = 0
CLLBC = 1
CLLB  = 1
Test  = 0
NWaveforms  = 5000
PlotSignals = False
PlotHistos  = True
FitFOM 		= True
BaselineSamples = 50
Threshold = 0.008
XMin = 0
XMax = 35
YMin = 0
YMax = 3
# Bins for scatter plot
nXBins = 5*(XMax-XMin)
nYBins = 100*(YMax-YMin)
# FOR CLYC
# W1 and W2 start and end sample values are relative to the trigger
# NEEDS SETTING CORRECTLY
if (CLYC):
	W1Start = -2
	W1End   = range(W1Start+20,W1Start+100, 5)
	#print W1End
	W2Start = 120
	W2End   = range(W2Start+10, W2Start+100, 10)
	#print W2End
# FOR CLLB/CLLBC
if (CLLB) or (CLLBC):
	W1Start = 10
	W1End = range(W1Start+40, W1Start+60, 2)
	W2Start = 200
	W2End   = range(W2Start+100, W2Start+600, 40)

if (Test):
	W1Start = 20
	W1End = range(W1Start+26, W1Start+28,2)
	W2Start = 120
	W2End   = range(W2Start+320, W2Start+360, 10)
f, cols, rows = OpenFile(1)

WaveformA      = []
WaveformMatrix = []
#TriggerList    = []

WaveformA = ReadBinaryFile(f, NWaveforms, cols, rows)
# 'WaveformMatrix' is a nxm array of waveforms, where n = number of waverforms and m = number of samples/waverform
Matrix = np.array(WaveformMatrix)
# Subtract baseline	
Matrix = BaselineSubtraction(Matrix, BaselineSamples)
# Get trigger point for each waveform in Matrix
TriggerIndexList = GetTriggerPosition(Matrix, Threshold) 

if (PlotSignals):
	for P in range(len(W1End)):
		for Q in range(len(W2End)):
			# Plot the signals with the initial Window Integral Boundaries.
			PlotWaveforms(Matrix, TriggerIndexList, W1Start, W1End[P], W2Start, W2End[Q])



# Get W1 , W2 and Full Integrals. Convert to floats from strings.
SignalPeak = map(float, GetSignalPeak(Matrix))
W1Integral = map(float, GetW1Integral(Matrix, TriggerIndexList, W1Start, W1End[0]))
W2Integral = map(float, GetW2Integral(Matrix, TriggerIndexList, W2Start, W2End[0]))
FullIntegral = (map(float, GetFullIntegral(Matrix)) )

# Remove Zero Values
SignalPeak, W1Integral, W2Integral, FullIntegral = RemoveZeroValues(SignalPeak, W1Integral, W2Integral, FullIntegral)
# Calculate PSD Factor 

PSDFactor = GetPSDFactor(SignalPeak, W1Integral, W2Integral)

# Make Scatter plot. 1 time is to allow user to set projection limits.
FOMPlot = plt.figure(figsize=(8,8), facecolor='white')

figScat = plt.figure(figsize=(14,6), facecolor='white')
plt.figure(figScat.number)
plt.subplot(1,2,1)

scatterhist = MakeScatterPlot(figScat, PSDFactor, FullIntegral, nXBins, nYBins, XMin,XMax, YMin, YMax, 1)
# Set projection limits
LowLimits, HighLimits = SetProjectionLimits()
# GetPeakRow asks for the centre of the gamma and neutron peak in the scatter plot
# and uses this to get the column of data to calculate the maximum index in the column.

GammaRow, NeutronRow = GetPeakRow(scatterhist, nXBins, XMin, XMax)

# Now loop

print "Looping over Window Integrals... "
FirstFOMMap = True
FOMMap = np.zeros( (len(W1End), len(W2End)) )

for P in range(len(W1End)):	
	for Q in range(len(W2End)):
		Continue = 'R'
		print W1End[P], W2End[Q]
		#while Continue == 'R' or Continue == 'r':
			# Plot Waveforms if requested
		

		# Get W1 , W2 and Full Integrals. Convert to floats from strings.
		SignalPeak = map(float, GetSignalPeak(Matrix))
		W1Integral = map(float, GetW1Integral(Matrix, TriggerIndexList, W1Start, W1End[P]))
		W2Integral = map(float, GetW2Integral(Matrix, TriggerIndexList, W2Start, W2End[Q]))
		FullIntegral = (map(float, GetFullIntegral(Matrix)) )
		SignalPeak, W1Integral, W2Integral, FullIntegral = RemoveZeroValues(SignalPeak, W1Integral, W2Integral, FullIntegral)
		# Calculate PSD Factor 
		PSDFactor = GetPSDFactor(SignalPeak, W1Integral, W2Integral)
		# Make Scatter plot. 1 time is to allow user to set projection limits.
		plt.figure(figScat.number)
		plt.subplot(1,2,1)
		scatterhist = MakeScatterPlot(figScat, PSDFactor, FullIntegral, nXBins, nYBins, XMin,XMax, YMin, YMax, 0)
		plt.figure(figScat.number)
		plt.subplot(1,2,2)
		plt.cla()
		ProjY, GammaPeak, NeutronPeak = ProjectY(scatterhist, LowLimits, HighLimits, nXBins, XMax, XMin, GammaRow, NeutronRow)
		if (FitFOM):
			mean1, mean2, sigma1, sigma2 = FitGaussian(ProjY, GammaPeak, NeutronPeak) 
		#	print mean1, mean2, sigma1, sigma2   			
			FOM = float('{0:.3f}'.format( abs(mean2-mean1)/(2.35*(abs(sigma1)+abs(sigma2)) ) ) )
			print 'FOM: ', FOM
			prevFOM = FOM
			# Catch fitting screw-ups
			if (FOM < 0 or FOM>5):
				# Where there's a fitting screw-up, just assume the previous FOM value. Save fitting image.
				FOM = prevFOM
				plt.savefig('Window1_'+str(W1End[P])+'Window2_'+str(W2End[Q])+'_'+str(FOM)+'.pdf')
			FOMMap[P,Q] = FOM
			np.savetxt('FOM_Interim.txt', FOMMap, fmt='%1.3f',delimiter=',')
			FomFig = plt.figure(FOMPlot.number)
			ax = plt.gca()
			
			# Scale X and Y tick labels by steps
			ax.set_xticklabels(map(str, W1End[:]) )
			ax.set_yticklabels(map(str, W2End[:]) )
			plt.imshow(FOMMap, cmap='hot', interpolation='nearest')	

			ax.set_aspect('auto')
			plt.xlabel('W1 Window')
			plt.ylabel('W2 Window')
			plt.draw()	
			plt.pause(0.001)
			#Continue = raw_input('Repeat? [R/r to repeat.  Enter to Continue.] ')	


plt.imshow(FOMMap)
plt.colorbar()		
Title = raw_input('Enter Title for Figure: ')
plt.xlabel('W1 Window')
plt.ylabel('W2 Window')
plt.title(Title)
np.savetxt(Title+'.txt', FOMMap, fmt='%1.3f',delimiter=',')
plt.savefig(Title+'.pdf', bbox_inches='tight')
plt.show()
