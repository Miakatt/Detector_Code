# Talk to Tektronix PSU. Increment voltage from 0 to 30V.
# Talk to Keithley Multimeter and read the current.
# Produce accurate IV curve.
# Fit line to data over a given current threshold (Thr) and calculate x for y = 0 = mx + c
# i.e. Turn On Point = -c/m
# Save pdf of plot and create Report file containing the temperature, turn on point, voltage and currents for recreating the plot



import visa
import sys
import os
import time
import numpy as np
import re
import matplotlib.pyplot as plt

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


# ======================== KEITHLEY COMMANDS ======================================

def InitializeKeithley():
	Keithley.write("*CLS")
	time.sleep(0.2)
	Keithley.timeout = 1000
	time.sleep(0.2)
	SetToCurrentMeasurement()
	SetCurrentRange()
	Model = Keithley.query("*IDN?")
	print("******************************************************************************")
	print("          "+Model )
	print("******************************************************************************")
# -----------------------------------------------------------------------------------

def SetToCurrentMeasurement():
	print "Setting Keithley To Current DC Measurement"
	Keithley.write(":CONF:CURR:DC")


def SetCurrentRange():
	print "Setting Current Range to mA"
	Keithley.write(":CURR:RANGe 0.01")

def MeasureCurrent(Samples):
	Q = []
#	Keithley.write(":SAMPle:COUNt %d" % int(Samples))
#	print (Keithley.query(":SAMPle:COUNt?"))
	for i in range(Samples):
		ThisCurrent = (Keithley.query(":READ?"))
		Q.append(float(ThisCurrent))

	return '{0:.5f}'.format(1e6*(sum(Q)/Samples))
# ======================== END OF KEITHLEY COMMANDS ================================


# =================================================================
# ====================START PROGRAM HERE ==========================
# =================================================================

# Necessary variables
SubIndex    = []
SubCurrents = []
SubVOLTS    = []
Thr         = 3


# Set Voltages to scan
VOLTS = np.arange(22, 29.5, 0.2)



# CREATE VISA INSTRUMENT FOR TEKTRONIX PSU
rm = visa.ResourceManager()
ResourceList = rm.list_resources()
Tek = rm.open_resource('USB0::0x0699::0x0391::C011169::INSTR')


# CREATE VISA INSTRUMENT FOR KEITHLEY MULTIMETER
rm = visa.ResourceManager()
ResourceList = rm.list_resources()
Keithley = rm.open_resource('ASRL2::INSTR')




# Initialize Tektronix power supply.
InitializeTektronix()
InitializeKeithley()

time.sleep(1)
Currents = []
for i, val in enumerate(VOLTS):
	volts = str(val)+'V'
	SetVoltage(volts)
	time.sleep(2)
	current = MeasureCurrent(10)
	print 'Current: %f uA' % float(current)
	Currents.append(float(current))
	time.sleep(0.2)




# Plot the data

plt.figure(figsize=(8,6), facecolor='white')
plt.plot(VOLTS, Currents, 'ro')
plt.xlim(VOLTS[0]-1, VOLTS[-1]+1)
plt.xlabel('Bias Voltage (V)')
plt.ylabel('Current (uA)')
plt.draw()
plt.pause(0.01)

# Get subset of current and voltage data from Current > Thr
SubIndex[:]    = [x for x,v in enumerate(Currents) if v>Thr]
SubCurrents[:] = [v for x,v in enumerate(Currents) if v>Thr]
for i, v in enumerate(SubIndex):
	SubVOLTS.append(VOLTS[v])
	
# Fit line to the subset data
fit = np.polyfit(SubVOLTS, SubCurrents, 1)
p = np.poly1d(fit)
plt.plot(VOLTS, p(VOLTS), 'r--')
plt.axhline(linewidth=1, linestyle='--', color='k')
plt.ylim(-1, max(Currents)*1.05)
plt.draw()
plt.pause(0.001)

# Get equation of line fit, p. Calculate turn on point
c = p[0]
m = p[1]
TurnOnPoint = '{0:.3f}'.format(-c/m)
print 'Turn On Point: %s V' % TurnOnPoint

# Save the plot with title and sensible filename
Title = raw_input('Enter Title for Figure: ')
plt.title(Title)
plt.draw()
plt.pause(0.001)
SaveTitle = Title.replace(" ", "_") # I f&@king love Python!
plt.savefig(SaveTitle+'.pdf', bbox_inches='tight')

# Create and fill Report File
Temp = raw_input('Enter Temperature:  ')
ReportFile = open('IV_Curves'+Temp+'C.txt', "a")
ReportFile.write(str(Temp)+', '+str(TurnOnPoint)+', ')
for i in range(len(VOLTS)):
	ReportFile.write(str(VOLTS[i])+', ')
for i in range(len(Currents)):
	ReportFile.write(str(Currents[i])+',')	
ReportFile.write('\n')	
ReportFile.close()


#Turn off SiPM Power
PowerOn(False)

plt.show()


