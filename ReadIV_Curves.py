# READS THE IV CURVES FROM IV_CURVES_<Temp>C.txt and calculates
# the turn on point by one of 3 methods. 
# order = 1: Linear fit. Set ThrL and ThrH to the linear part of the data (curent thresholds)
# order = 2: Quadratic fit. Set ThrL and ThrH to the lower part of the data (curent thresholds)
# order = 3. Linear fit to Cubic curve. Set ThrL and ThrH to the lower part of the curve (Current). 
# Set x1 and x2 to high voltages (29.0 and 31.0, for example). This fits a linear function to the linear part
# of the cubic fit. 

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
#------------------------------------------------------------
def ReadCSV(arg):
	filename = sys.argv[arg]
	print filename

	with open(filename, 'r') as csvfile:
		read = csv.reader(csvfile, delimiter=',')
		for row in read:
			ClearLists()
			Temperature.append(float(row[0]))
			TOP.append(float(row[1]))
			for v in range(38):
				#print float(row[2+v]), "  ",  float(row[40+v]) 
				Volts.append(float(row[2+v]))
				Current.append(float(row[40+v]))
			print "-----"
			MakeIVPlots(Temperature[-1], Volts, Current)		
	return Temperature, TOP, Volts, Current		

#------------------------------------------------------------
def MakeIVPlots(Temp, Volts, Current):
	plt.figure(fig1.number)	
	plt.plot(Volts, Current, 'o', label=Temp)
	plt.hold(True)
	plt.legend(loc='center left')
	plt.ylabel('Current (uA)')
	plt.xlabel('Bias Voltage (V)')
	plt.title('IV Curve at %.1f$^{\circ}$C. Turn On estimated by %d-Order fit' % (Temp, order))
	#plt.xlim(21,30)
	plt.draw()
	plt.pause(0.01)
	CalTOP = FitLine(Volts, Current)
	return CalTOP

#------------------------------------------------------------

def GetDerivative(Current):
	Grad = np.gradient(Current)
	return Grad
#------------------------------------------------------------
def MakeTOPPlot():
	plt.figure(fig2.number)	
	plt.plot(Temperature, TOP, 'd-.')
	ax2 = plt.gca()
	ax2.ticklabel_format(useOffset=False, style='plain')
	plt.ylabel('Turn On Voltage (V)')
	plt.xlabel('Temperature $^{\circ}$C')
	plt.xlim(-21, 30)
	plt.hold(True)
	plt.draw()
	plt.pause(0.01)
#------------------------------------------------------------

def expofit(x, a,b,c):
	return a*np.exp(b*x) + c

#------------------------------------------------------------

def ClearLists():
	del Volts[:]
	del Current[:]
	del SubIndex[:]
	del SubCurrents[:]
	del SubVOLTS[:]
#------------------------------------------------------------
def FitLine(VOLTS, Currents):
	TurnOnPoint = 0
	print 'Thresholds: ' , ThrL, ThrH
	# Get subset of current and voltage data from Current > Thr
	SubIndex[:]    = [x for x,v in enumerate(Currents) if v>ThrL and v<ThrH]
	SubCurrents[:] = [v for x,v in enumerate(Currents) if v>ThrL and v<ThrH]
	for i, v in enumerate(SubIndex):
		SubVOLTS.append(VOLTS[v])
	
	# Fit line to the subset data
	if (order is not 'e'):
		fit = np.polyfit(SubVOLTS, SubCurrents, order)
		p = np.poly1d(fit)
	#	plt.plot(VOLTS, p(VOLTS), 'r--', linewidth=1)
	elif (order == 'e'):
		popt, pcov = curve_fit(expofit, SubVOLTS, SubCurrents)
	#	plt.plot(SubVOLTS, expofit(SubVOLTS, *popt), 'r--', linewidth=1)

	plt.axhline(linewidth=1, linestyle='--', color='b')
#	plt.ylim(-1, max(Currents)*1.05)
#	plt.draw()
#	plt.pause(0.001)

	# Get equation of line fit, p. Calculate turn on point
	if (order == 1):
		c = fit[1]
		m = fit[0]
		TurnOnPoint = '{0:.3f}'.format(-c/m)
	elif (order ==2):
		c = fit[2]
		b = fit[1]
		a = fit[0]
		Nom = 	-b+math.sqrt(b**2 - (4*a*c))
		Denom = (2*a)
		TurnOnPoint = '{0:.6}'.format(Nom/Denom)
		#print Nom, Denom, TurnOnPoint
		#TurnOnPoint = '{0:.3f}'.format( (-b+math.sqrt(b**2 - (4*a*c))/(2*a)) )
		print a, b, c
	elif (order==3):
		# Boundary parameters for the linear fit to the 3rd order polynomial if "order == '3' "
		x1 = 27.0
		x2 = 28.0
		ThirdOrderSubVOLTS = []
		ThirdOrderSubIndex = []
		ThirdOrderSubCurrents = []
		
		d = fit[3]	
		c = fit[2]
		b = fit[1]
		a = fit[0]
		# print a, b, c, d
		#Get equation of linear part of 3rd order polynomial. 
		# By evaluating at two points > 29V.
		# e.g y1 = a(x1)^3 + b(x1)^2 + c(x1) + d
		#     y2 = a(x2)^3 + b(x2)^2 + c(x2) + d
		#      m = (y2-y1)/(x2-x1)
		#      c = y1 - mx1 or c = y2 - mx2
		# Extraploate to y = 0.

		y1 = a*(x1)**3 + b*(x1)**2 + c*x1 + d
		y2 = a*(x2)**3 + b*(x2)**2 + c*x2 + d
		m = (y2-y1)/(x2-x1)
		c = y1 - m*x1

		ThirdOrderSubIndex[:] = [x for x,v in enumerate(VOLTS) if v>x1 and v<x2]
		for i, v in enumerate(ThirdOrderSubIndex):
			ThirdOrderSubVOLTS.append(VOLTS[v])
			ThirdOrderSubCurrents.append(m*ThirdOrderSubVOLTS[i] + c)
		fit = np.polyfit(ThirdOrderSubVOLTS, ThirdOrderSubCurrents, 1)
		p = np.poly1d(fit)
		#plt.plot(VOLTS, p(VOLTS), 'k--', linewidth=1)
		print m, c
		TurnOnPoint = '{0:.3f}'.format(-c/m)
	elif (order == 'e'):
		a = popt[0]
		b = popt[1]
		c = popt[2]
		print a, b, c
		TurnOnPoint = '{0:.6}'.format(math.log(-c/a)/b)

	print 'Turn On Point: %s V' % TurnOnPoint
	CalculatedTOP.append(TurnOnPoint)
	return CalculatedTOP
#------------------------------------------------------------
#------------------------------------------------------------

# Boundary parameters for the linear fit to the 3rd order polynomial if "order == '3' "
# x1 = 27.0
# x2 = 28
# Order = 1 for line, for 2nd order polynomial, 'e' for exponential.
order = 3
if order == 1:
	ThrL = 10
	ThrH = 20
	# for 20C data.
	#ThrL = 3
	#ThrH = 6
if order == 2:
	ThrL = 0.1
	ThrH = 20
if order == 3:
	ThrL = 0
	ThrH = 40

CalculatedTOP = []
Temperature = []
TOP         = []
Volts       = []
Current     = []
SubIndex    = []
SubVOLTS    = []
SubCurrents = []
fig1 = plt.figure(figsize=(12,6), facecolor='white')
#fig2 = plt.figure(figsize=(8,6), facecolor='white')
fig3 = plt.figure(figsize=(8,6), facecolor='white')
		

Temperature, TOP, Volts, Current = ReadCSV(1)			
plt.figure(fig3.number)
plt.plot(Temperature, CalculatedTOP, 'ro')
ax3 = plt.gca()
ax3.ticklabel_format(useOffset=False, style='plain')
# Get Spread of data
ymin, ymax = ax3.get_ylim()
plt.ylim(ymin-(ymin*0.0001), ymax*1.0001)
plt.ylabel('Turn On Voltage (V)')
plt.xlabel('Temperature $^{\circ}$C')
plt.title('Variability of SiPM Turn On Point extrapolation at %.1f $^{\circ}$C.' % Temperature[0] )
plt.draw()
#MakeTOPPlot()
FloatCalTop = map(float, CalculatedTOP)
TOPVariability =  max(FloatCalTop)-min(FloatCalTop)
print TOPVariability
plt.show()