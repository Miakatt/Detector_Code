# Talk to Tektronix PSU. Increment voltage from 0 to 30V.
# Talk to Agilent Scope and read the trigger rate.
# Produce accurate turn on curve.



import visa
import sys
import os
import time
import numpy as np
import re
import matplotlib.pyplot as plt
from termcolor import colored
# ======================== TEKTRONIX COMMANDS ======================================
def InitializeTektronix():
	ClearInstr()
	time.sleep(0.2)
	Tek.timeout = 1000
	time.sleep(0.2)
	SetMaxVoltage('29.6V')
	time.sleep(0.2)
	PowerOn(True)
	Model = Tek.query("*IDN?")
	print("******************************************************************************")
	print("                 "+Model)
	print("******************************************************************************")
	time.sleep(0.2)
	V = (str(VOLTS[0])+'V')
	SetVoltage(V)

# -----------------------------------------------------------------------------------

def PowerOn(ON):
	if (ON is True):
		Tek.write("OUTP ON")
	elif (ON is False):
		Tek.write("OUTP OFF")

# -----------------------------------------------------------------------------------
def ClearInstr():
	Tek.write("*CLS")

# -----------------------------------------------------------------------------------
def GetVoltage():
	print "Getting Voltage"
	Tek.query("MEAS:VOLT?")
	voltage = Tek.query("FETC:VOLT?")
	
# -----------------------------------------------------------------------------------

def SetVoltage(volts):
	# Check request isn't larger than maximum allowed voltage
	OVP = GetMaxVoltage()
	
	if float(volts[0:4]) > float(OVP):
		print "Requested voltage is too high."
		return  
	print "Setting Voltage to %s" % volts
	Tek.write("SOURce:VOLTage %s" % volts)

# -----------------------------------------------------------------------------------

def SetMaxVoltage(MaxVoltage):
	print "Setting Maximum Voltage to %s" % MaxVoltage
	Tek.write("VOLT:PROT %s" % MaxVoltage)
# -----------------------------------------------------------------------------------

def GetMaxVoltage():
	OVP = Tek.query("VOLT:PROT?")
	return OVP
# -----------------------------------------------------------------------------------


# ======================== AGILENT COMMANDS ======================================

def InitializeAgilent():
	Agilent.write("*CLS")
	time.sleep(0.2)
	Agilent.timeout = 1000
	time.sleep(0.2)
	Agilent.term_chars = ';'
	Model = Agilent.query("*IDN?")
	print("******************************************************************************")
	print("          "+Model )
	print("******************************************************************************")
# -----------------------------------------------------------------------------------


def DisplayChannel(channel, state):
	print("Channel = %s" % channel)
	Agilent.write(":CHANnel%s:DISPlay %s" % (channel, state) )
	print("DISPLAY CHANNEL: %s " % Agilent.query(":CHANnel%s:DISPlay?" % channel))

# -----------------------------------------------------------------------------------

def SetTriggerSweep(Mode):
	Agilent.write("TRIGger:SWEep %s" % Mode)

# -----------------------------------------------------------------------------------

def SetTriggerChannel(channel):
	if (channel != 'EXTernal'):
		Agilent.write(":TRIGger:SOURce CHANnel %s" % channel )
		
# -----------------------------------------------------------------------------------

def SetTriggerLevel(level):
	Agilent.write(":TRIGger:LEVel%s" % level)	


# -----------------------------------------------------------------------------------

def SetTimeBaseScale(timeperdiv):
	Agilent.write(":TIMebase:SCALe %s" % timeperdiv)

# -----------------------------------------------------------------------------------

def SetChannelScale(channel, scale):
	Agilent.write(":CHANnel%s:SCALe %s" % (channel, scale))

# -----------------------------------------------------------------------------------

def GetFrequency(channel):
	Agilent.write(":MEASure:FREQuency%s" % channel)
	Freq = Agilent.query(":MEASure:FREQuency?")
	return float(Freq)
# -----------------------------------------------------------------------------------

def GetAmplitude(channel):
	Agilent.write(":MEASure:VAMPlitude%s" % channel)
	Amp = Agilent.query(":MEASure:VAMPlitude?")
	if float(Amp)>100:
		Amp = '0'
		print 'Amp ', Amp		
		return Amp
	print 'Amp ', Amp	
	return float(Amp)
# -----------------------------------------------------------------------------------

# Necessary variables

Amp = []
FreqSample = []
SubIndex    = []
SubCurrents = []
SubVOLTS    = []
Thr         = 3


# Set Voltages to scan
VOLTS = np.arange(24, 29.5, 0.01)

# CREATE VISA INSTRUMENT FOR TEKTRONIX PSU
rm = visa.ResourceManager()
ResourceList = rm.list_resources()
Tek = rm.open_resource('USB0::0x0699::0x0391::C011169::INSTR')


# CREATE VISA INSTRUMENT FOR AGILENT OSCILLOSCOPE
rm = visa.ResourceManager()
ResourceList = rm.list_resources()
Agilent = rm.open_resource('USB0::0x0957::0x17A2::MY52491398::INSTR')




# Initialize Tektronix power supply.
InitializeTektronix()
InitializeAgilent()
SetTriggerChannel(4)
SetTriggerLevel(20e-3)
DisplayChannel(1,1)
DisplayChannel(2,0)
DisplayChannel(3,0)
DisplayChannel(4,1)
SetChannelScale(1, 50e-3)
SetChannelScale(4, 2)
SetTimeBaseScale(50e-6)
SetTriggerSweep("AUTO")


fig = plt.figure(figsize=(8,6), facecolor='white')
ax = plt.gca()
#try:
for i, val in enumerate(VOLTS):
	volts = str(val)+'V'
	SetVoltage(volts)
	time.sleep(1)
	Amp.append(GetAmplitude(1))
		
plt.plot(VOLTS, Amp, 'ro')
ax = plt.gca()	
plt.draw()
plt.pause(0.01)	

print Amp
# Create and fill Report File
Temp = raw_input('Enter Temperature:  ')
ReportFile = open('BreakDownVoltage.txt', "a")
ReportFile.write(str(Temp)+', ')
for item in VOLTS:
	ReportFile.write(" %s," % item)
for item in Amp:
	ReportFile.write(" %s," % item)

ReportFile.write('\n')	
ReportFile.close()
plt.show()
#except:
	#Turn off SiPM Power
PowerOn(False)
Agilent.close()
Tek.close()
plt.show()


