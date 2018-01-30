#Make FOM Map of Neutron and Gamma Data.
# Usage:python
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


def OpenFile(inputfile):
	filename = inputfile
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
	WaveformA      = []
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
		MaxSignal.append(max(Matrix[M,:]))
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

# =======================================================================	
# ======================================================================	
# 								Start program here	
# =======================================================================
# =======================================================================	

NeutronDataFiles  = ['CsI_22C.bin'] # ['CLLBC_20_Neutrons.bin'] #['CLLB-10_Neutrons.bin','CLLB00_Neutrons.bin','CLLB-20_Neutrons.bin','CLLB40_Neutrons.bin','CLLB50_Neutrons.bin']
GammaDataFiles    = ['CsI_0C.bin']  # ['CLLBC_20_Gammas.bin'] #['CLLB-10_Gammas.bin','CLLB00_Gammas.bin','CLLB-20_Gammas.bin','CLLB40_Gammas.bin','CLLB50_Gammas.bin']


W1Axis = []
W2Axis = []
WaveformMatrix = []

my_cmap = plt.cm.jet
my_cmap.set_under('w', 1)

CLYC  = 0
CLLBC = 1
CLLB  = CLLBC
CLLB_Delayed = 0
Test  = 0
SampleRate = 8e-9
#++++++++++++++++++++++++++++++++++
NWaveforms  = 1000
#+++++++++++++++++++++++++++++++++
PlotSignals = False
PlotHistos  = True
FitFOM 	 = True
PlotScatter = True
# Required to get colorbar only once on the live FOM Map
FirstFOMMap = True
FirstScatter = True

BaselineSamples = 50
Threshold = 0.005
XMin = 0
XMax = 50
YMin = 0
YMax = 20
# Bins for scatter plot
nXBins = 5*(XMax-XMin)
nYBins = 200*(YMax-YMin)
# FOR CLYC
# W1 and W2 start and end sample values are relative to the trigger
# NEEDS SETTING CORRECTLY
if (CLYC):
	W1StepSize = 2
	W2StepSize = 20

	W1Start = -5
	W1End   = range(W1Start+10,W1Start+70, W1StepSize)
	#print W1End
	W2Start = 40
	W2End   = range(W2Start+20, W2Start+260, W2StepSize)
     
	#print W2End
# FOR CLLB/CLLBC    20,12,220,20, 100, 200, 660, 40)
if (CLLB) or (CLLBC):
     W1StepSize = 5
     W2StepSize = 10
     W1Start = 0
     W1End = range(W1Start+100, W1Start+130, W1StepSize)
     W2Start = 130
     W2End   = range(W2Start+200, W2Start+600, W2StepSize)
     
if (CLLB_Delayed):
    W2StartStepSize = 10
    W2EndStepSize  = 40
    W1Start = 20
    W1End   = 40
    W2Start = range(100, 201,  W2StartStepSize)
    W2End   = range(200, 401,  W2EndStepSize)
   
#if (Test):
#    W1StepSize = 5
#    W2StepSize = 10
#    W1Start = 10
#    W1End = range(W1Start+40, W1Start+50, W1StepSize)
#    W2Start = 60
#    W2End   = range(W2Start+400, W2Start+500, W2StepSize)

if (CLLB or CLLBC or CLYC):
    FOMMap = np.zeros( (len(W1End), len(W2End)) )
elif (CLLB_Delayed):
    FOMMap = np.zeros( (len(W2Start), len(W2End)) )

# Integrating Sample Numbers
if (CLLB or CLLBC or CLYC):
    XAxis = np.arange(W2Start+W2End[0], W2Start+W2End[-1], W2StepSize)
    YAxis = np.arange(W1Start+W1End[0], W1Start+W1End[-1], W1StepSize)	
elif (CLLB_Delayed):
    XAxis = np.arange(W2End[0],   W2End[-1],   W2EndStepSize) 
    YAxis = np.arange(W2Start[0], W2Start[-1], W2StartStepSize)	

# Create Figures
figScat = plt.figure(figsize=(14,6), facecolor='white')
FOMPlot = plt.figure(figsize=(8,8),  facecolor='white')
plt.figure(figScat.number)
plt.subplot(121)


for FileNo, File in enumerate(NeutronDataFiles):
    del WaveformMatrix[:]
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Read In Neutron Data in to array
    N, cols, rows = OpenFile(NeutronDataFiles[FileNo])
    WaveformA = ReadBinaryFile(N, NWaveforms, cols, rows)
    # Convert in to numpy array for easier manipulation.
    NeutronMatrix = np.array(WaveformMatrix)
    
    del WaveformMatrix[:]
    # Read In Gamma Data in to array 
    G, cols, rows = OpenFile(GammaDataFiles[FileNo])
    WaveformA = ReadBinaryFile(G, NWaveforms, cols, rows)
    # Convert in to numpy array for easier manipulation.
    GammaMatrix = np.array(WaveformMatrix)

    #Subtract Baseline offset from Neutron Data
    NeutronMatrix = BaselineSubtraction(NeutronMatrix, BaselineSamples)
    #Subtract Baseline offset from Gamma Data
    GammaMatrix = BaselineSubtraction(GammaMatrix, BaselineSamples)
    #Get list of Neutron trigger sample indexes where data matrix row crosses threshold
    NeutronTriggerIndexList = GetTriggerPosition(NeutronMatrix, Threshold) 
    print 'Neutrons: Zero Trigger Samples ' , len(NeutronTriggerIndexList) - np.count_nonzero(NeutronTriggerIndexList)
    #Get list of Gamma trigger sample indexes where data matrix row crosses threshold
    GammaTriggerIndexList = GetTriggerPosition(GammaMatrix, Threshold) 
    print 'Gammas: Zero Trigger Samples ' , len(GammaTriggerIndexList) - np.count_nonzero(GammaTriggerIndexList)

    # Remove signals where trigger index = 0
    NeutronMatrix = NeutronMatrix[np.nonzero(NeutronTriggerIndexList)]
    NeutronTriggerIndexList = NeutronTriggerIndexList[np.nonzero(NeutronTriggerIndexList)]
    print 'Neutrons: Zero Trigger Samples ' , len(NeutronTriggerIndexList) - np.count_nonzero(NeutronTriggerIndexList)

    GammaMatrix = GammaMatrix[np.nonzero(GammaTriggerIndexList)]
    GammaTriggerIndexList = GammaTriggerIndexList[np.nonzero(GammaTriggerIndexList)]
    print 'Gammas: Zero Trigger Samples ' , len(GammaTriggerIndexList) - np.count_nonzero(GammaTriggerIndexList)

    # Plot Signals if enabled
    if (PlotSignals and (CLLB or CLLBC or CLYC)):
        for P in range(len(W1End)):
            for Q in range(len(W2End)):
                # Plot the signals with the initial Window Integral Boundaries.
                PlotWaveforms(NeutronMatrix, GammaMatrix, NeutronTriggerIndexList, GammaTriggerIndexList, W1Start, W1End[P], W2Start, W2End[Q])
    # Plot Signals if enabled
    elif (PlotSignals and (CLLB_Delayed)):
        for P in range(len(W2Start)):
            for Q in range(len(W2End)):
                # Plot the signals with the initial Window Integral Boundaries.
                PlotWaveforms(NeutronMatrix, GammaMatrix, NeutronTriggerIndexList, GammaTriggerIndexList, W1Start, W1End, W2Start[P], W2End[Q])



#==================================================================
#                     Start Looping over W1 and W2
#==================================================================


    print "Looping over Window Integrals... "
# Comment out the relevant loop.
    for P in range(len(W2Start)):
        for Q in range(len(W2End)):
            # for P in range(len(W1End)):	
            #     for Q in range(len(W2End)):
            gc.collect()
            if (CLLB or CLLBC or CLYC):
                print 'W1 End Point : %d.   W2 End Point : %d' %(W1End[P], W2End[Q])
            elif (CLLB_Delayed):
                print 'W2 Start Point : %d.  W2 End Point: %d' %(W2Start[P], W2End[Q])
            # Get Signal Peaks, W1Integrals, W2Integrals and FullIntegrals from neutron dat
            NeutronSignalPeak   = map(float, GetSignalPeak( NeutronMatrix))
            if (CLLB or CLLBC or CLYC):
                NeutronW1Integral   = map(float, GetW1Integral( NeutronMatrix, NeutronTriggerIndexList, W1Start, W1End[P]))
                NeutronW2Integral   = map(float, GetW2Integral( NeutronMatrix, NeutronTriggerIndexList, W2Start, W2End[Q]))
            elif (CLLB_Delayed):    
                NeutronW1Integral   = map(float, GetW1Integral( NeutronMatrix, NeutronTriggerIndexList, W1Start, W1End))
                NeutronW2Integral   = map(float, GetW2Integral( NeutronMatrix, NeutronTriggerIndexList, W2Start[P], W2End[Q]))
            NeutronFullIntegral = map(float, GetFullIntegral( NeutronMatrix )) 
            
            NeutronIntegrals = [NeutronSignalPeak, NeutronW1Integral, NeutronW2Integral, NeutronFullIntegral]
            # Get Signal Peaks, W1Integrals, W2Integrals and FullIntegrals from gamma data
            GammaSignalPeak   = map(float, GetSignalPeak( GammaMatrix ))
            if (CLLB or CLLBC or CLYC):
                GammaW1Integral   = map(float, GetW1Integral( GammaMatrix, GammaTriggerIndexList, W1Start, W1End[P]))
                GammaW2Integral   = map(float, GetW2Integral( GammaMatrix, GammaTriggerIndexList, W2Start, W2End[Q]))
            elif (CLLB_Delayed):    
                GammaW1Integral   = map(float, GetW1Integral( GammaMatrix, GammaTriggerIndexList, W1Start, W1End))
                GammaW2Integral   = map(float, GetW2Integral( GammaMatrix, GammaTriggerIndexList, W2Start[P], W2End[Q]))
            GammaFullIntegral = map(float, GetFullIntegral( GammaMatrix ))

            GammaIntegrals = [GammaSignalPeak, GammaW1Integral, GammaW2Integral, GammaFullIntegral]
            # Remove W2Integral Zero Data Entries (and corresponding entries in other arrays) in Neutron Data
            NeutronIntegrals = RemoveZeroValues( NeutronIntegrals )
            # Remove W2Integral Zero Data Entries (and corresponding entries in other arrays) in Gamma Data
            GammaIntegrals = RemoveZeroValues( GammaIntegrals ) 

            # Get Neutron PSDFactor List (PSDFactors for all waveforms) for Neutrons
            NeutronPSDFactor = GetPSDFactor(NeutronIntegrals)

            # Get Gamma PSDFactor List (PSDFactors for all waveforms) for Gammas
            GammaPSDFactor   = GetPSDFactor(GammaIntegrals)

            # ===== PLOT NEUTRON AND GAMMA SCATTER PLOTS SEPERATELY ====================
            # Make Neutron Scatter Plot from PSD Factors and FullIntegral
            #		print 'Size of Inputs:  Peak  W1Integral  W2Integral  FullIntegral'
            #		print 'Size of Inputs: ', len(NeutronIntegrals[0]), len(NeutronIntegrals[1]), len(NeutronIntegrals[2]), len(NeutronIntegrals[3])
            # Neutronscatterhist = MakeScatterPlot(figScat, NeutronPSDFactor, NeutronIntegrals[3], nXBins, nYBins, XMin,XMax, YMin, YMax, FirstScatter)
            # FirstScatter = False
            # Make Gamma Scatter Plot from PSD Factors and FullIntegral
            #Gammascatterhist   = MakeScatterPlot(figScat, GammaPSDFactor, GammaIntegrals[3], nXBins, nYBins, XMin,XMax, YMin, YMax, FirstScatter)

            # ======== Or Merge the PSDFactors and FullIntegral arrays and make scatterplot. =============
            CombinedPSDFactors = np.concatenate(( NeutronPSDFactor, GammaPSDFactor) )
            CombinedFullInts   = np.concatenate(( NeutronIntegrals[3], GammaIntegrals[3]) )

            CombinedScatterHist = MakeScatterPlot(figScat, CombinedPSDFactors, CombinedFullInts, nXBins, nYBins, XMin,XMax, YMin, YMax, FirstScatter)
            FirstScatter = False
            PlotScatter = False
            # On First Loop, Request Projection Limits
            if (P ==0 and Q == 0 and FileNo == 0):
                LowLimits, HighLimits, NeutronLimit = SetProjectionLimits()
            #    LowLimits = 9.84677419355 
            #    HighLimits = 33.2661290323 
            #    NeutronLimit = 16.7661290323
                print LowLimits, HighLimits, NeutronLimit
                
            # ==================================================================
            # ProjectY of Neutron Data
            # Fit Gaussian, Get Mean and Sigma
            # NeutronIntegrals[3] is the list of Full Integrals for the Neutron Data
            Neutron2DHist = Make2DHist(NeutronPSDFactor, NeutronIntegrals[3], nXBins, nYBins, XMin,XMax, YMin, YMax)
            NeutronProjY, NeutronPeakEstimate = ProjectY(Neutron2DHist, LowLimits, HighLimits, NeutronLimit, nXBins, XMax, XMin, True)
            NeutronFit = FitGaussian(NeutronProjY, NeutronPeakEstimate) 
            print 'Neutron Peak : ' , NeutronPeakEstimate
            if NeutronFit[0] == 0:
                print 'Skipping this one...'
                continue
            # ProjectY of Gamma Data
            # Fit Gaussian, Get Mean and Sigma
            # GammaIntegrals[3] is the list of Full Integrals for the Gamma Data
            Gamma2DHist = Make2DHist(GammaPSDFactor, GammaIntegrals[3], nXBins, nYBins, XMin,XMax, YMin, YMax)
            GammaProjY, GammaPeakEstimate = ProjectY(Gamma2DHist, LowLimits, HighLimits, NeutronLimit, nXBins, XMax, XMin, False)
            GammaFit = FitGaussian(GammaProjY, GammaPeakEstimate) 
            print 'Gamma Peak : ', GammaPeakEstimate
           
            #Overlay ProjectionY data and Gaussians
            if (PlotHistos):
                PlotHistograms(figScat, GammaProjY, GammaPeakEstimate, GammaFit, NeutronProjY, NeutronPeakEstimate, NeutronFit)
                # Calculate FOM and put in to 2d array
                # Amp  = popt[0]
                # Mean = popt[1]
                # Sigma = popt[2]
                GammaMean    = GammaFit[1]
                GammaSigma   = GammaFit[2]
                NeutronMean  = NeutronFit[1]
                NeutronSigma = NeutronFit[2]
            if NeutronFit[0] == 0:
                print 'Skipping this one...'
                continue
            if GammaFit[0] == 0:
                print 'Skipping this one...'
                continue
            print 'Neutron Mean %.2f.  Gamma Mean %.2f.  Neutron Sigma %.2f.  Gamma Sigma %.2f ' %(NeutronMean, GammaMean, NeutronSigma, GammaSigma)
            FOMMap = PlotFOM(NeutronMean, GammaMean, NeutronSigma, GammaSigma)

	#End of Loop

	# SaveFOMMap(FOMMap)
	# Save FOM Map as pdf and txt
    plt.imshow(FOMMap)
    if (FirstFOMMap):
#		plt.colorbar()		
     #   plt.xlabel('Window 2 Integral (Starting at %1.2f us)' %  (W2Start*SampleRate*1e6))
     #   plt.ylabel('Window 1 Integral (Starting at %1.2f us)' %  (W1Start*SampleRate*1e6 ))
        FirstFOMMap = False
    Underscore = NeutronDataFiles[FileNo].index('_')
    plt.title(NeutronDataFiles[FileNo][0:Underscore]+'_FOM_Map')
    np.savetxt(NeutronDataFiles[FileNo][0:Underscore]+'.txt', FOMMap, fmt='%1.3f',delimiter=',')
    plt.savefig(NeutronDataFiles[FileNo][0:Underscore]+'.pdf', bbox_inches='tight')
    # Delete FOM Map and other variables. 
    del FOMMap
    if (CLLB or CLLBC or CLYC):
        FOMMap = np.zeros( (len(W1End), len(W2End)) )
    elif (CLLB_Delayed):
        FOMMap = np.zeros( (len(W2Start), len(W2End)) )

    del WaveformMatrix
    WaveformMatrix = []
    

    FOMPlot.clf()
    figScat.clf()
    plt.close('all')
    gc.collect()

   # Create Figures
    figScat = plt.figure(figsize=(14,6), facecolor='white')
    FOMPlot = plt.figure(figsize=(8,8),  facecolor='white')
    plt.figure(figScat.number)
    plt.subplot(121)

print 'THE END'






















