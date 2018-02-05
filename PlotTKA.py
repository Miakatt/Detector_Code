#/usr/bin/env python


import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def OpenFile(arg):
    fullfilename = sys.argv[arg]
    filename, file_ext = os.path.splitext(sys.argv[arg])
    print(filename, file_ext)
    return fullfilename, file_ext

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
def ReadSPE(fullfilename):

    print('Reading SPE')
    x = 0
    data = []
    with open(fullfilename,'r') as f:
        for line in f:

            if (x>1) and (x<=8190):

                datum = line

                data.append(float(datum))
            x = x+1

    return data

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def FitGaussian(XClick, channels, data):
    mean = int(XClick[0][0])
    amp  = XClick[0][1]
    sigma = np.sqrt(mean)
    slope = -5
    intercept = 200
    p0 = [amp, mean, sigma, slope, intercept]
    R = 300
    fitRange = channels[mean-R:mean+R]
    fitdata  = data[mean-R:mean+R]

    popt, pcov = curve_fit(gausWithBackground, fitRange, fitdata, p0=p0)
    print ('Mean: ', popt[1])
    print ('FWHM: ', 2.35*popt[2])
    Resolution = 235 * popt[2] / popt[1]
    print (Resolution)
    ymin, ymax = plt.ylim()
    xmin, xmax = plt.xlim()
    plt.ylim(-10, ymax)
    plt.xlim(-10, 5000)
    plt.xlabel('Energy (Ch)', fontsize=10)
    plt.ylabel('Counts', fontsize=10)
    plt.plot(fitRange, gausWithBackground(fitRange, *popt), 'black', linewidth=2, )
    plt.text(1000, ymax - 100, r'Resolution: %2.1f%% +/- 0.1%%' % abs(Resolution), color='black', fontsize=10)
#   plt.hold(True)
    plt.draw()

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def gausWithBackground(x, *p0):
    amp, mean, sigma, slope, intercept = p0

    return amp*np.exp(-(x-mean)**2/(2*sigma**2)) + slope*x + intercept

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def gaus(x, *p0):
    amp, mean, sigma = p0
    return amp*np.exp(-(x-mean)**2/(2*sigma**2))

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo




filename, ext = OpenFile(1)
data= ReadSPE(filename)
channels = range(len(data))
print (channels)
plt.plot(channels, data, 'r')
plt.draw()
plt.pause(0.01)
xclick = plt.ginput(1)

FitGaussian(xclick, channels, data)
#plt.figure(num=1, facecolor='white', figsize=(8,6))
#plt.plot(channels, data, 'r', label = 'Label')

plt.legend(loc='best')
plt.savefig('CS137_2.png')
plt.show()