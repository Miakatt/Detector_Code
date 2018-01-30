import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import csv
import sys
import os

def ReadCSV(arg):
	filename = sys.argv[arg]
	print filename
	with open(filename, 'r') as csvfile:
		read = csv.reader(csvfile, delimiter=',')
		for row in read:
			Temp.append(int(row[0]))
			Bias.append(float(row[1]))
			Mean122.append(float(row[2]))
			Peak122.append(float(row[3]))
			Reso122.append(float(row[4]))
			Peak344.append(float(row[5]))
			Mean344.append(float(row[6]))
			Reso344.append(float(row[7]))
			Peak662.append(float(row[8]))
			Mean622.append(float(row[9]))
			Reso662.append(float(row[10]))
			Peak778.append(float(row[11]))
			Mean778.append(float(row[12]))
			Reso778.append(float(row[13]))
			Peak964.append(float(row[14]))
			Mean964.append(float(row[15]))
			Reso964.append(float(row[16]))
			Peak1408.append(float(row[17]))
			Mean1408.append(float(row[18]))
			Reso1408.append(float(row[19]))

	return filename



def FitCurve(Temp, Bias, order):
	pol2 = np.polyfit(Temp, Bias, order)
	return pol2




 #==================================================================
 #==================================================================


# Font
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')



Temp     = []
Bias 	 = []
Peak122  = []
Peak245  = []
Peak344  = []
Peak662  = []
Peak778  = []
Peak964  = []
Peak1408 = []

Mean122  = []
Mean245  = []
Mean344  = []
Mean662  = []
Mean778  = []
Mean964  = []
Mean1408 = []

Reso122  = []
Reso245  = []
Reso344  = []
Reso662  = []
Reso778  = []
Reso964  = []
Reso1408 = []


filename = ReadCSV(1)

RTI = Temp.index(20)
Dev122  = [100*(x  - Peak122[RTI])/Peak122[RTI] for x in Peak122]
#Dev245  = [100*(x  - Peak245[RTI])/Peak245[RTI] for x in Peak245]
Dev344  = [100*(x  - Peak344[RTI])/Peak344[RTI] for x in Peak344]
Dev662  = [100*(x  - Peak662[RTI])/Peak662[RTI] for x in Peak662]
Dev778  = [100*(x  - Peak778[RTI])/Peak778[RTI] for x in Peak778]
Dev964  = [100*(x  - Peak964[RTI])/Peak964[RTI] for x in Peak964]
Dev1408 = [100*(x - Peak1408[RTI])/Peak1408[RTI] for x in Peak1408]



majorLocator = MultipleLocator(10)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(2)

majorYLocator = MultipleLocator(0.1)
minorYLocator = MultipleLocator(0.05)
#  PLOT BIAS VERSES TEMPERATURE

pol2 = FitCurve(Temp, Bias, 2)
p = np.poly1d(pol2)
xp = np.linspace(-20,50,70)
xerr = 1.0
fig1 = plt.figure(figsize=(10,8), facecolor='white')
ax1 = fig1.add_subplot(111)
plt.errorbar(Temp, Bias, 0, xerr, 'ro', label='Data')
plt.plot(xp, p(xp), 'b--', label='Polynomial Fit')
plt.xlim([min(Temp)-5, max(Temp)+5])
plt.ylim([28, 29.0])
plt.xlabel('Temperature $^{\circ}$C', fontsize=14)
plt.ylabel('Bias Voltage (V)', fontsize=14)
# X axis tick markers
ax1.xaxis.set_major_locator(majorLocator)
ax1.xaxis.set_major_formatter(majorFormatter)
ax1.xaxis.set_minor_locator(minorLocator)
# Y axis tick markers
ax1.yaxis.set_major_locator(majorYLocator)
ax1.yaxis.set_minor_locator(minorYLocator)
ax1.annotate('$V(T) = %.1E x^{2}  %.1E x + %.2f$' %(pol2[0], pol2[1], pol2[2]), xy=(-10,28.75), fontsize=16)


file, ext = os.path.splitext(filename)
plt.title(file+" Bias vs Temperature Curve")
plt.savefig(file+'_Bias_Curve.pdf', bbox_inches='tight')


# PLOT PHOTOPEAK POSITIONS BEFORE LINEARISATION
majorYLocator = MultipleLocator(500)
minorYLocator = MultipleLocator(100)

fig2 = plt.figure(figsize=(10,8), facecolor='white')
ax2 = fig2.add_subplot(111)
plt.plot(Temp, Peak122, 'go--', label='122keV')
#plt.plot(Temp, Peak245, 'mo--', label='245keV')
plt.plot(Temp, Peak344, 'bo--', label='344keV')
plt.plot(Temp, Peak662, 'ro--', label='662keV')
plt.plot(Temp, Peak778, 'ko--', label='778keV')
plt.plot(Temp, Peak964, 'yo--', label='964keV')
plt.plot(Temp, Peak1408, 'co--', label='1408keV')
plt.xlim([min(Temp)-5, max(Temp)+5])
plt.xlabel('Temperature $^{\circ}$C', fontsize=14)
plt.ylabel('Peak Position (keV)', fontsize=14)
# X axis tick markers
ax2.xaxis.set_major_locator(majorLocator)
ax2.xaxis.set_major_formatter(majorFormatter)
ax2.xaxis.set_minor_locator(minorLocator)
# Y axis tick markers
ax2.yaxis.set_major_locator(majorYLocator)
ax2.yaxis.set_minor_locator(minorYLocator)

plt.legend(loc='best')
plt.draw()
plt.title(file+" Photopeak Positions Pre-Linearisation")
plt.savefig(file+'_Peak_Positions.pdf', bbox_inches='tight')


# PLOT PHOTOPEAK CHANNEL

fig3 = plt.figure(figsize=(10,8), facecolor='white')
ax3 = fig3.add_subplot(111)
plt.plot(Temp, Mean122, 'go--', label='122keV')
plt.plot(Temp, Mean245, 'mo--', label='245keV')
plt.plot(Temp, Mean344, 'bo--', label='344keV')
plt.plot(Temp, Mean662, 'ro--', label='662keV')
plt.plot(Temp, Mean778, 'ko--', label='778keV')
plt.plot(Temp, Mean964, 'yo--', label='964keV')
plt.plot(Temp, Mean1408, 'co--', label='1408keV')
plt.ylim([-2, 2])
plt.xlim([min(Temp)-5, max(Temp)+5])
plt.xlabel('Temperature $^{\circ}$C', fontsize=14)
plt.ylabel('Mean Channel (Normalised to 20$^{\circ}$C)', fontsize=14)
plt.legend(loc='best')

ax3.xaxis.set_major_locator(majorLocator)
ax3.xaxis.set_major_formatter(majorFormatter)
ax3.xaxis.set_minor_locator(minorLocator)

plt.title(file+" Photopeak Channel Post-Linearisation")
plt.savefig(file+'_Peak_Channels.pdf', bbox_inches='tight')


# PLOT RESOLUTIONS AGAINST TEMPERATURE
majorYLocator = MultipleLocator(1)
minorYLocator = MultipleLocator(0.1)

fig4 = plt.figure(figsize=(10,8), facecolor='white')
ax4 = fig4.add_subplot(111)
plt.plot(Temp, Reso122, 'go--', label='122keV')
#plt.plot(Temp, Reso245, 'mo--', label='245keV')
plt.plot(Temp, Reso344, 'bo--', label='344keV')
plt.plot(Temp, Reso662, 'ro--', label='662keV')
plt.plot(Temp, Reso778, 'ko--', label='778keV')
plt.plot(Temp, Reso964, 'yo--', label='964keV')
plt.plot(Temp, Reso1408, 'co--', label='1408keV')
plt.xlim([min(Temp)-5, max(Temp)+5])
plt.xlabel('Temperature $^{\circ}$C', fontsize=14)
plt.ylabel('Resolution %', fontsize=14)
# X axis tick markers
ax4.xaxis.set_major_locator(majorLocator)
ax4.xaxis.set_major_formatter(majorFormatter)
ax4.xaxis.set_minor_locator(minorLocator)
# Y axis tick markers
ax4.yaxis.set_major_locator(majorYLocator)
ax4.yaxis.set_minor_locator(minorYLocator)
plt.legend(loc='best')
plt.title(file+" Photopeak Resolutions Post-Linearisation")
plt.savefig(file+'_Resolutions_v_Temperature.pdf', bbox_inches='tight')

plt.show()
