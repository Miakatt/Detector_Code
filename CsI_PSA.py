# Scans through the .bin data for CsI located in /Users/ajbell/Work/Crystal_Studies/CsI/
# Automatically finds the bin files, extracts the temperature from the file name and
# reads the data.
# All signals are read in to a matrix and Baseline correction, time adjustment and saturation removal is done
# on all simultaneously (clever, eh?)
# Carries out a charge comparison PSA on either a random N selection (Set by NoOfWaveforms), or choses
# the largest N signals


import PSD
import glob
import os
import re
import random
import numpy as np
import matplotlib.pyplot as plt
from collections import deque


def PlotWaveforms(Waveforms, delay, F):

#	for wf in range(Waveforms.shape[0]):
	plt.plot( Waveforms , label=Temperature[F]+'C')
#	plt.ylim(-0.01,20)
	plt.legend(loc='best')
	plt.draw()
	plt.pause(delay)
	#plt.cla()


def TimeAlign(Waveforms, TriggerList):
	# Use to determine if bad data gets plotted
	Status = True

	TimeAlignedWaveforms = np.zeros(1000)
	for i, n in enumerate(TriggerList):
		try:
	#		print 'n  ' ,  n, '  i ', i
			TimeAlignedWaveforms += Waveforms[i, n-100:n+900]
			#plt.plot(Waveforms[i, n-100:n+900])
			#plt.draw()
			#plt.pause(0.1)
		except ValueError:
			continue
			#print 'Bad data'
			#Status = False

	return TimeAlignedWaveforms, Status


def RemoveSaturation(AllWaveforms, SatLimit):
	SatIndex = []
	Rows = np.shape(AllWaveforms)[0]
	print 'Rows ', Rows
	for i in range(Rows-1):
#		print (AllWaveforms[i,:])
		if max(AllWaveforms[i,:]) >= SatLimit:
			SatIndex.append(i)

	print 'Saturation Index ' , SatIndex
	npSatIndex  = np.array(SatIndex)
	AllWaveforms  = np.delete(AllWaveforms, npSatIndex, 0)
	print ' Length ', len(AllWaveforms)
	return AllWaveforms

def Initialise():
	DataList = []
	Temperature = []
	PSAFactor = []

	return DataList, Temperature, PSAFactor
# ----------------------------------------------------
#                  START PROGRAM
# ----------------------------------------------------

fig1 = plt.figure(figsize=(8,6), facecolor='white')
ax1 = fig1.add_subplot(111)
colormap = plt.cm.gist_ncar
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, 18)])

fig2 = plt.figure(figsize=(8,6), facecolor='white')

PlotSignals = True
RandomSample = True
# Chose whether or not to normalize the time aligned signals.
Normalise  = True
# Threshold at which to align the timing of the signals.
Threshold = 0.4
BaselineSamples = 20
ReadWaveforms = 1000
NoOfWaveforms = 500
SatLimit = 0.195 #Volts
# Delay = time (sec) between displaying waveforms, if PlotWaveforms is True
delay = 0.1
# Charge integration windows
Q1Min = 90
Q1Max = 290
Q2Min = 291
Q2Max = 491
# Get list of .bin files in directory
dirs = ["/Users/ajbell/Work/Crystal_Studies/CsI/Crystal_1", "/Users/ajbell/Work/Crystal_Studies/CsI/Crystal_2"]


DataList, Temperature, PSAFactor = Initialise()

for d, ThisDir in enumerate(dirs):
	os.chdir(ThisDir)
	del DataList[:]
	for file in glob.glob("*.bin"):
		DataList.append(file)
		T = re.search('_(.*)C.bin', file)
		print T.group(1)
		Temperature.append(T.group(1))
	print "Files to read: " , DataList
	print "Temperatures:  " , Temperature

#	NormalizeTimeAlignedWaveforms = np.empty([len(DataList), NoOfWaveforms])
	for F, file in enumerate(DataList):
		print "Reading " ,file
		## Read in N waveforms from bin file.
		## For each waveform with a peak above PeakThreshold,
		## do a baseline correction and time alignment of the sample at which
		## the signal crosses 'Threshold'.
		## Normalise each signal
		## Add signals together and divide by NoOfWaveforms to get the average signal.
		FID, cols, rows   = PSD.OpenFile(file)
		#print cols, rows
			## ReadBinary(..) returns a NxM matrix of waveforms. FID = fileID, NoOfWaveforms.
		AllWaveforms  = np.array(PSD.ReadBinaryFile(FID, ReadWaveforms, cols, rows ))

		# Remove waveforms that show saturation. i.e. Ones in which the maximum value is 100mV
		AllWaveforms = RemoveSaturation(AllWaveforms, SatLimit)
		if (RandomSample):
			# Randomly select NoOfWaveforms from all the waveforms.
			Selection = random.sample(AllWaveforms, np.shape(AllWaveforms)[0])
		if (not RandomSample):
			## Or... find the NoOfWaveforms that have the greatest sum (Largest Signals)
			SumAllWaveforms = PSD.GetFullIntegral(AllWaveforms)
			## Get NoOfWaveforms indices of the AllWaveforms
			N_Max = SumAllWaveforms.argsort()[-NoOfWaveforms:]
			Selection = AllWaveforms[N_Max]

		# Convert to Numpy Array
		Waveforms  = np.array(Selection)
		#print "Length of Waveforms : ", len(Waveforms)
		Waveforms = PSD.BaselineSubtraction(Waveforms, BaselineSamples)

		# Get peak value of each waveform
		PeakList =  PSD.GetSignalPeak(Waveforms)

		# Normalize each signal.
		for wf , val in enumerate(Waveforms):
			Waveforms[wf] = val/PeakList[wf]
		del PeakList
		# Get trigger sample (the sample at which the waveform crosses threshold)
		TriggerList = PSD.GetTriggerPosition(Waveforms, Threshold)
		#print TriggerList
		TimeAlignedWaveforms, Status = TimeAlign(Waveforms, TriggerList)


		# Normalize each timealigned averaged signal.
		if (Normalise):
			Peak = float('{0:0.8f}'.format(max(TimeAlignedWaveforms[:])))
		 	TimeAlignedWaveforms /= Peak

		if (PlotSignals and Status is True):
			plt.figure(fig1.number)
			ax1 = plt.gca()
			ax1.ticklabel_format(style='plain')
			plt.axvline(Q1Min, color='k', linestyle='--')
			plt.axvline(Q1Max, color='k', linestyle='--')
			plt.axvline(Q2Min, color='k', linestyle='--')
			plt.axvline(Q2Max, color='k', linestyle='--')

		PlotWaveforms(TimeAlignedWaveforms, delay, F)

		## Get the integrals of Q2 and Q1 for each individual waveform
		Q1 = float('{0:.6f}'.format(sum(TimeAlignedWaveforms[Q1Min:Q1Max])))
		Q2 = float('{0:.6f}'.format(sum(TimeAlignedWaveforms[Q1Min:Q2Max])))
		PSAFactor.append(Q1+Q2/Q1)

		# Tidy up
		del Waveforms
		del TimeAlignedWaveforms

	## Plot in a graph of Counts vs Q2/Q1
	plt.figure(fig2.number)
	if d==0:
		plt.plot(Temperature, PSAFactor,'ro')
	elif d==1:
		plt.plot(Temperature, PSAFactor,'bo')

	#Convert temperatures from string to floats
	Float_Temperatures = []
	for T, V in enumerate(Temperature):
		Float_Temperatures.append(float(V))
	print Float_Temperatures
	fit = np.polyfit(Float_Temperatures, PSAFactor, 1)
	p = np.poly1d(fit)
	print p
	plt.figure(fig2.number)
	if d==0:
		plt.plot(Float_Temperatures, p(Float_Temperatures), 'r--', label ='Crystal #1')
	if d==1:
		plt.plot(Float_Temperatures, p(Float_Temperatures), 'b--', label ='Crystal #2')
	plt.legend(loc='best')
#	plt.xlabel('Temperature ($^{\circ}$C)')
#	plt.ylabel('PSA Factor (Q1+Q2)/Q2')
	plt.draw()
	plt.xlim(min(Float_Temperatures)-5, max(Float_Temperatures)+5)
	print ' RESETTING EVERYTHING '
	DataList, Temperature, PSAFactor = Initialise()

#******* End loop through directories ******
plt.figure(fig2.number)
plt.xlabel('Temperature $^{\circ}$C')
plt.ylabel('PSA Factor (Q1+Q2)/Q1')
plt.ticklabel_format(style='plain', useOffset=False)
if RandomSample:
	plt.title('Temperature determination by PSA on CsI:(Th) - Random Sampling')
elif not RandomSample:
	plt.title('Temperature determination by PSA on CsI:(Th) - Largest Signals Sample')
plt.draw()
## Get Integrals of Q2 and Q1 for the averaged signal.
## Plot value as Q2/Q1 vs temperature.


## Calculate PSA value
##


plt.show()
