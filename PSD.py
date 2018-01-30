import numpy as np
import sys
import time
import struct
import codecs
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from operator import truediv
from scipy.optimize import curve_fit
from scipy.ndimage.interpolation import zoom
import scipy.signal as ss

import gc



# -----------------------------------------------------------


def OpenFile(inputfile):
	filename = inputfile
	#print filename
	f = open (filename, "rb")
	binary_datapoint = f.read(4)
	cols = struct.unpack("i", binary_datapoint)[0]
	#print "Columns ", cols

	binary_datapoint = f.read(4)
	rows = struct.unpack("i", binary_datapoint)[0]
	#print "Rows ", rows

	return f, cols, rows


# -----------------------------------------------------------


def ReadBinaryFile(f, iterations, cols, rows):
	if 'WaveformA' in locals():
		del WaveformA[:]
	else:
		WaveformA = []
	if 'WaveformMatrix' in locals():
		del WaveformMatrix[:]
	else:
		WaveformMatrix = []

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
		WaveformMatrix = StackWaveforms(WaveformA[0:cols], WaveformMatrix)
	return WaveformMatrix

# -----------------------------------------------------------

def StackWaveforms(WaveformA, WaveformMatrix):
   # WaveformA = smooth(WaveformA, 10)
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
	for M in range(len(Matrix)):
		MaxSignal.append(float('{0:0.8f}'.format(max(Matrix[M,:]))))
	return MaxSignal


# -----------------------------------------------------------


def GetPSDFactor(Integrals=[]):
	SignalPeak = Integrals[0]
	W1Integral = Integrals[1]
	W2Integral = Integrals[2]

	if (CLYC) or (Test):
		PSDFactor = map(truediv, W1Integral,W2Integral)
	elif (CLLB) or (CLLBC) or (Test) or (CLLB_Delayed):
		PSDFactor = map(truediv, W1Integral, W2Integral)
		# Try doing maximum of waveform/W2 Integral
		#PSDFactor = SignalPeak/W2Integral


	return PSDFactor

# -----------------------------------------------------------

def RemoveZeroValues(Integrals = []):
	SignalPeak   = Integrals[0]
	W1Integral   = Integrals[1]
	W2Integral   = Integrals[2]
	FullIntegral = Integrals[3]
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

	Integrals = [SignalPeak, W1Integral, W2Integral, FullIntegral]
	return Integrals

# -----------------------------------------------------------

def PlotWaveforms( NeutronMatrix, GammaMatrix, NeutronTriggerIndexList, GammaTriggerIndexList, W, X, Y, Z ):

	for wf in range(NeutronMatrix.shape[0]):
		plt.plot(NeutronMatrix[wf],'r-')
		plt.axvline(NeutronTriggerIndexList[wf]+W, color='r', linestyle='--')
		plt.axvline(NeutronTriggerIndexList[wf]+X  , color='b', linestyle='--')
		plt.axvline(NeutronTriggerIndexList[wf]+Y, color='g', linestyle='--')
		plt.axvline(NeutronTriggerIndexList[wf]+Z  , color='c', linestyle='--')
		plt.axhline(Threshold, color='k', linestyle='--')
		plt.ylim(-0.01, 0.5)
		plt.hold(True)
		plt.draw()
		plt.pause(0.01)
		plt.cla()


	#for wf in range(GammaMatrix.shape[0]):
	#	plt.plot(GammaMatrix[wf],'r-')
	#	plt.axvline(GammaTriggerIndexList[wf]+W, color='r', linestyle='--')
	#	plt.axvline(GammaTriggerIndexList[wf]+X  ,color='b', linestyle='--')
	#	plt.axvline(GammaTriggerIndexList[wf]+Y, color='g', linestyle='--')
	#	plt.axvline(GammaTriggerIndexList[wf]+Z  , color='c', linestyle='--')
	#	plt.axhline(Threshold, color='k', linestyle='--')
	#	plt.ylim(-0.01, 0.05)
	#	plt.hold(False)
	#	plt.draw()
	#	plt.pause(0.01)
	#	plt.cla()


# -----------------------------------------------------------

def MakeScatterPlot(figScat, PSDFactor, FullIntegral, nXBins, nYBins, XMin,XMax, YMin, YMax, First):
	# Create scatter plot of PSDFactor vs Full Int
	# and get user to set Y Projection limits.
	plt.figure(figScat.number)
	plt.subplot(121)
	h = plt.hist2d(FullIntegral, PSDFactor, vmin=1, cmap=my_cmap, bins=[nXBins,nYBins], range=[[XMin,XMax],[YMin, YMax]])
#	if (First):
#		plt.colorbar()
	plt.ylim(YMin,YMax)
	plt.xlabel('Full Integral (arb. units)')
	plt.ylabel('PSD Factor')
	plt.draw()
	plt.pause(0.001)

	return np.array(h[0])
# -----------------------------------------------------------

def Make2DHist(PSDFactor, FullIntegral, nXBins, nYBins, XMin,XMax, YMin, YMax):

	hist = plt.hist2d(FullIntegral, PSDFactor, bins=[nXBins,nYBins], range=[[XMin,XMax],[YMin, YMax]])
	return np.array(hist[0])



# -----------------------------------------------------------

def SetProjectionLimits():
	print 'Please click lower limit: '
	XLow  = plt.ginput(1, timeout=0)[0]
	#print XLow[0]
	LowLine  = plt.axvline(XLow[0],  color='r', linestyle='--')
	plt.draw()
	plt.pause(0.001)

	print 'Please click upper limit: '
	XHigh = plt.ginput(1, timeout=0)[0]
	#print XHigh[0]
	HighLine = plt.axvline(XHigh[0], color='r', linestyle='--')
	plt.draw()
	plt.pause(0.001)

	print 'Please click lower neutron cut limit: '
	NeutronLimit = plt.ginput(1, timeout=0)[0]
	#print XHigh[0]
	NeutronLine = plt.axvline(NeutronLimit[0], color='r', linestyle='--')
	plt.draw()
	plt.pause(0.001)


	return XLow[0], XHigh[0], NeutronLimit[0]

# -----------------------------------------------------------


def smooth(x,window_len=10,window='hanning'):
    x = np.array(x)
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

def ProjectY(scatterhist, LowLimits, HighLimits, NeutronLimit, nXBins, XMax, XMin, Neutrons):
	# if Neutrons is True, Only project the scatterhist beterrn NeutronLimit and HighLimits
	# This is because we're getting more gammas from the neutron source than we get neutrons!

	#Get Bin number from the projection limits.
	Low  = int(LowLimits* nXBins/(XMax-XMin))
	High = int(HighLimits* nXBins/(XMax-XMin))
	NLimit = int(NeutronLimit*nXBins/(XMax-XMin))
	#print 'Low ' , Low
	#print 'High ', High
	# Make Y Projection of data within the set limits
	if not Neutrons:
		ProjY = sum(scatterhist[Low:High,:])
		PeakEstimate   = np.argmax(ProjY)
	if Neutrons:
		ProjY = sum(scatterhist[NLimit:High,:])
		PeakEstimate   = np.argmax(ProjY)
	return ProjY, PeakEstimate

# -----------------------------------------------------------

def gaus(x, *p0):
    amp, mean, sigma = p0
    return amp*np.exp(-(x-mean)**2/(2*sigma**2))
# -----------------------------------------------------------

def FitGaussian(ProjY, PeakEst):
#	ProjY = smooth(ProjY,3)
	Sigma = 1
	Amp = ProjY[PeakEst]
	PeakWidth = int(round(0.1*PeakEst))


	#print 'Peaks: ', GammaProjPeak, NeutronProjPeak
	#print 'Amplitudes: ', Amp1,  Amp2
	#print 'Widths: ' ,PeakGammaWidth, PeakNeutronWidth

	p0=[Amp, PeakEst, Sigma]
	fitrange = range(PeakEst-PeakWidth, PeakEst+PeakWidth)
	fitdata = ProjY[PeakEst-PeakWidth: PeakEst+PeakWidth]
	try:
		popt, pcov = curve_fit(gaus, fitrange, fitdata,  p0=p0)
		mean  = popt[1]
		sigma = popt[2]

	except RuntimeError:
		popt = [0,0,1]
		return popt

	except ValueError:
		print "Unable to fit peak. "
		popt = [0,0,1]
		return popt
	except TypeError:
		print "Unable to fit"
		popt = [0,0,1]
		return popt

	return popt

# -----------------------------------------------------------

def GetPeakRow(scatterhist, nYBins, YMin, YMax):
	print 'Select peak of gamma data points.'
	GammaRow = int(round(plt.ginput(1, timeout=0)[0][0]) * (nXBins/(XMax-XMin)) )
	print 'Select peak of neutron data points.'
	NeutronRow = int(round(plt.ginput(1, timeout=0)[0][0]) * (nXBins/(XMax-XMin)) )
	print 'Gamma Row:  ', GammaRow
	print 'Neutron Row: ', NeutronRow
	return GammaRow, NeutronRow

# -----------------------------------------------------------

def PlotHistograms(figScat, GammaProjY, GammaPeakEstimate, GammaFit, NeutronProjY, NeutronPeakEstimate, NeutronFit):
	plt.figure(figScat.number)
	plt.subplot(122)
	ax = plt.gca()
	ax.cla()
	plt.plot(NeutronProjY,'b', linewidth=1)
	plt.plot(GammaProjY,  'r', linewidth=1)
	data = range(GammaPeakEstimate-150, NeutronPeakEstimate+150)
	plt.plot(data, gaus(data, *NeutronFit), 'cyan', linewidth=2, )
	plt.plot(data, gaus(data, *GammaFit),   'green', linewidth=2, )
	plt.xlim(GammaPeakEstimate-100, NeutronPeakEstimate+100)
	plt.draw()
	plt.pause(0.001)

# -----------------------------------------------------------

def PlotFOM(NeutronMean, GammaMean, NeutronSigma, GammaSigma):
		FOM = float('{0:.3f}'.format(abs(NeutronMean-GammaMean)/(2.35*( abs(GammaSigma) + abs(NeutronSigma) )    ) ) )
		print 'FOM : ' ,  FOM
           #Plot FOM Map, Setting XAxis and YAxis Labels correctly
		FOMMap[P, Q] = FOM
		Underscore = NeutronDataFiles[FileNo].index('_')
		np.savetxt(NeutronDataFiles[FileNo][0:Underscore]+'_FOM_Interim.txt', FOMMap, fmt='%1.3f', delimiter=',')
		FOMFig = plt.figure(FOMPlot.number)
		ax = plt.gca()
		#ax.set_aspect('auto')

		# Add Tick Marks
		majorFormatter   = FormatStrFormatter('%1.2f')
		majorLocator     = MultipleLocator(1)
		ax.xaxis.set_major_locator(majorLocator)

		ax.set_xticklabels(map(str,  XAxis[:]) )
		ax.set_yticklabels(map(str,  YAxis[:]) )

		# Add Tick Marks
		majorFormatter = FormatStrFormatter('%1.2f')
	#	ax.xaxis.set_major_locator(majorXLocator)
	#	ax.yaxis.set_major_locator(majorYLocator)
		plt.imshow(FOMMap, cmap='hot', interpolation='nearest', extent=[XMin, XMax, YMin, YMax], aspect='auto')
		plt.title(NeutronDataFiles[FileNo][0:Underscore]+'_FOM_Map')
	#	plt.xlabel('Window 2 Integral (Starting at %1.2f us)' %  (W2Start*SampleRate*1e6))
	#	plt.ylabel('Window 1 Integral (Starting at %1.2f us)' %  (W1Start*SampleRate*1e6 ))


		plt.draw()
		plt.pause(0.0001)
		return FOMMap
