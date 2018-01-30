#!/usr/bin/env python


#  Alan Bell. April 2017
# Usage Linearity.py <Ref spe or TKA file> <Cs137 spe or TKA file>
# Linearity correction of Cs137 spectrum using either NG3 or Eu152 spectrum as a reference.


import sys
from future_builtins import map
from pylab import *
from scipy.optimize import curve_fit
import numpy.polynomial.polynomial as poly
from operator import truediv, sub
import os
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import time
import bisect
import pickle

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

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
            x = x+1

            datum = line
            if x==11:
            #    print datum
                # Get number of data points from file
                dataLength = int(datum[2:])
            # Read those data points from file
            if (x>11 and x<dataLength+11):

                data.append(float(datum))

    return data
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def ReadTKA(fullfilename):
    print('Reading TKA File')
    print(fullfilename)
    x = 0
    data = []
    with open(fullfilename,'r') as f:
        # Set number of data points from file
        dataLength = 8190

        for line in f:
            x = x+1
            datum = line

            # Read those data points from file
            if (x>2):

                data.append(float(datum))

    return data

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def PlotData(data):
    channels = range(len(data))

    plt.plot(channels, data, color='blue')
    plt.xlabel('Channel Number',fontsize=10)
    plt.ylabel('Counts', fontsize=10)
    ymin, ymax = ylim()
    plt.ylim(1,ymax)
    plt.semilogy()
    plt.xlim(0, len(data))
    # ginput(1) autofits the point when clicked
    # ginput(2) places a marker and waits for user to hit <Enter>
    x = plt.ginput(1)
    plt.hold(True)

    return channels, data, x

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def FitGaussian(XClick, channels, data):
    mean = XClick[0][0]
    amp  = XClick[0][1]
    popt = []
    sigma = np.sqrt(mean)

    try:
        for i in range(10):
            # fit reference peaks multiple times. Record mean position each time.
            # Calculate uncertainties by +/- std dev of means.
            RandomNoise = np.random.uniform(-mean*0.01,mean*0.01)
            LowerRange = int(mean+RandomNoise-sigma)
            UpperRange = int(mean+RandomNoise+sigma)
            LowerRange = int(mean+RandomNoise-sigma)
            UpperRange = int(mean+RandomNoise+sigma)
            fitRange = range(LowerRange, UpperRange)
            fitdata  = data[LowerRange:UpperRange]
            p0=[amp, mean, sigma]
            popt, pcov = curve_fit(gaus, fitRange, fitdata, p0=p0)
            errList.append(popt[1])
    except RuntimeError:
        popt = [0,0,0]
        # print(popt)
        return popt
    plt.plot(channels, data, 'blue')
    plt.plot(channels, gaus(channels, *popt), 'black', linewidth=1, )
#   plt.hold(True)
    plt.draw()

    return popt, pcov, errList

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def gaus(x, *p0):
    amp, mean, sigma = p0
    return amp*np.exp(-(x-mean)**2/(2*sigma**2))

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def FitQuadratic(energy, means, w):
	xq = np.poly1d(np.polyfit(means[:], energy[:], 2, w = w))
	return xq

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def FitLinear(energy,means, w):
    xlin = np.poly1d(np.polyfit(means[:], energy[:], 1, w = w))
    # Constrain fit to >=0 at y-intercept
    return xlin

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def gausWithBackground(x, *p0):
    amp, mean, sigma, slope, intercept = p0

    return amp*np.exp(-(x-mean)**2/(2*sigma**2)) + slope*x + intercept

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def PlotLinearity(energy, means, w):

    # Get a 2nd order polynomial and linear fit (of the first 4 points) for
    # the Eu152 mean peak values.
    xquad = FitQuadratic(energy, means,w)
    xlin  = FitLinear(energy, means, w)

    coefs = xquad.c

   # print 'Coefs:  ', coefs
    x = range(dataLength)
    if PlotNoise is False:
        # Plot both fits to the mean points
        #fig1 = figure(facecolor='white')
        if Source=='Eu152':
            plt.plot(means, energy, 'rx',  label='Eu152 Peaks')
        elif Source=='NG3':
            plt.plot(means, energy, 'rx',  label='NG3 Peaks')

        plt.plot(x, xquad(x), 'r-' ,    label='Quadratic Fit')
        plt.plot(x, xlin(x),  'b-.',    label='Linear Fit')
        plt.ylabel('Energy (keV)', fontsize=10)
        plt.xlabel('Channel Number', fontsize=10)

        plt.legend(loc='upper left', fontsize=10)
        plt.xlim(0, dataLength)

        ax = gca()
        if coefs[2] >=0:
        	ax.text(2000, 0, r'y = %.2f${\times}$10$^{-6}$x$^{2}$ + %.2fx  + %.2f ' %( coefs[0]*1e6, coefs[1], coefs[2] ) , fontsize=12 )
        elif coefs[2] < 0:
    		ax.text(2000, 0, r'y = %.2f${\times}$10$^{-6}$x$^{2}$ + %.2fx  - %.2f ' %( coefs[0]*1e6, coefs[1], abs(coefs[2]) ) , fontsize=12 )
        plt.draw()


    #fig1.savefig('Linearity.pdf', bbox_inches='tight')

    return xquad, xlin

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

###### Not currently implemented #######
def PlotDeviation(xquad, xlin):
    ch = range(dataLength)
    print 'PlotDeviation'
    QuadFit = np.array(xquad(ch))
    LinFit  = np.array(xlin(ch))
    Deviation = (QuadFit-LinFit)
    plt.plot(ch, Deviation, 'b--')
    plt.xlabel('Channel')
    plt.ylabel('Deviation from Linear (ch)')
    plt.draw()
    plt.pause(0.0001)
    return Deviation, ch

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def PlotResiduals(xlin, xquad, means, energy, mean_err):
    # Evaluate the quadratic at each mean. Calculate the difference
    # between the energy  and the quadratic
    residual = []
    pct = []
    for i, mean in enumerate(means):
        residual.append( energy[i] - xquad(mean) )
        pct.append(100*residual[-1]/energy[i])

    plt.figure(num=3, figsize=(9, 7), dpi=100, facecolor='w', edgecolor='k')
    plt.subplot(211)
    plt.errorbar(energy, residual, yerr=mean_err , fmt='rd--')
    plt.axhline(0,  linewidth=1, linestyle='--', color='k')
    for i, val in enumerate(energy):
        if residual[i] >= 0:
            plt.text(energy[i], residual[i]+0.7, str(round(residual[i],1))+'keV\n'+str(val)+'keV', fontsize=9, multialignment='center')
        elif residual[i] < 0:
            plt.text(energy[i], residual[i]-0.7, str(round(residual[i],1))+'keV\n'+str(val)+'keV', fontsize=9, multialignment='center')
    plt.xlabel('Isotope Energy (keV)')
    plt.ylabel('Deviation from $f(E)$ (keV)')
    ymin, ymax = plt.ylim()
    plt.ylim(ymin-2, ymax+2)
    plt.draw()
    plt.pause(0.001)

    plt.subplot(212)
    plt.plot(energy, pct, 'rd--')
    plt.axhline(0,  linewidth=1, linestyle='--', color='k')
    for i, val in enumerate(energy):
        if residual[i] >= 0:
            plt.text(energy[i], pct[i]+0.2, str(round(pct[i],1))+'%\n'+str(val)+'keV', fontsize=9, multialignment='center')
        elif residual[i] < 0:
            plt.text(energy[i], pct[i]-0.2, str(round(pct[i],2))+'%\n'+str(val)+'keV', fontsize=9, multialignment='center')
    plt.xlabel('Isotope Energy (keV)')
    plt.ylabel('Deviation from $f(E)$ (%)')
    ymin, ymax = plt.ylim()
    plt.ylim(ymin-0.5, ymax+0.5)
    plt.draw()
    plt.pause(0.001)

    return residual, pct

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


def Linearize(Cs137data, xquad, xlin):

    chan_num = range(len(Cs137data))                #[0, 1, 2 .... dataLength]
    chan_center = [bin+0.5 for bin in chan_num]     #[0.5, 1.5, 2.5 .... dataLength + 0.5]
    chan_center2= [bin+1.5 for bin in chan_num]     #[1.5, 2.5, 3.5 .... dataLength + 1.5]



    # Python equivalent of a*x[i]*x[i] + b*x[i] * c  within a for loop.
    # Quadratic fit normalisation
    chan_cal  = xquad(chan_center)
    chan_cal2 = xquad(chan_center2)
    # Linear fit normalisation
    chan_lin   = xlin(chan_center)
    chan_lin2  = xlin(chan_center2)

    chan_quad_width = np.array(chan_cal2) -  np.array(chan_cal) # quadratic fit bin width
    chan_lin_width  = np.array(chan_lin2) -  np.array(chan_lin) # linear fit bin width

    # plt.figure(num=99, figsize = (8,6), facecolor='white')
    # c = np.arange(0,len(chan_quad_width))
    # plt.plot(c, chan_quad_width, 'g', c, chan_lin_width,'b')
    # plt.draw()
    # plt.pause(1)
    # plt.figure(1)

    NormCs137data = np.array(Cs137data) / np.array(chan_quad_width)
    LinCs137data  = np.array(Cs137data) / np.array(chan_lin_width)

    return chan_num, chan_lin, LinCs137data, chan_cal, NormCs137data

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def CalculateResolution(chan, Cs137data, color, placement):
    # Find the channel (>Ch 500) that hold the maximum data value. [0][0] takes it out of the list structure.

    meanCh = np.where(Cs137data==max(Cs137data[CUT:]))[0][0]

    # Get the mean of the photopeak (each channel does not equal one increment on 'x-axis'
    PhotoPeakMean = chan[meanCh]
    PhotoPeakAmp  = Cs137data[meanCh]
    sig = sqrt(PhotoPeakMean)
    slope = -0.001
    intercept = 5
    p1=[PhotoPeakAmp, PhotoPeakMean, sig, slope, intercept]
    # Fit to a range of data around the mean
    fitRange = chan[meanCh-LOWRANGE:meanCh+HIGHRANGE]
    fitData  = Cs137data[meanCh-LOWRANGE:meanCh+HIGHRANGE]

    # Fit the photo-peak with a linear background
    popt, pcov = curve_fit(gausWithBackground, fitRange, fitData, p0=p1)

    #popt, pcov = curve_fit(gaus, fitRange, fitData, p0=[amp, mean, sig])
    plt.plot(fitRange, gausWithBackground(fitRange, *popt), color, linewidth=2)

    #plt.plot(fitRange, gaus(fitRange, *popt), color, linewidth=2)
    plt.draw()
    plt.pause(0.01)
    figure(1)

    #  plt.plot(chan_lin, LinCs137data, 'blue', label='Corrected')
    plt.plot(fitRange, gausWithBackground(fitRange, *popt), color, linewidth=2, label='Gaussian Fit')
    #plt.plot(fitRange, gaus(fitRange, *popt), color, linewidth=2, label='Gaussian Fit')
    plt.legend(loc='upper right', fontsize=8)
    ymin, ymax = plt.ylim()
    xmin, xmax = plt.xlim()
    ylim(-10, ymax)
    xlim(-10, 1000)
    xlabel('Energy (keV)', fontsize=10)
    ylabel('Counts', fontsize=10)

    # Calculate resolution and error
    perr = np.sqrt(np.diag(pcov))
    Resolution = 2.35*popt[2]/popt[1]
    # Get fractional error of mean and FWHM
    MeanErr = perr[1]/popt[1]
    FWHMErr = perr[2]/popt[2]

    ResolutionErr = Resolution * sqrt(FWHMErr**2 + MeanErr**2)
    # annotate plot with resolution.
    plt.text(150, ymax-placement, r'Resolution: %2.1f%% +/- 0.1%%' % abs(100*Resolution) , color=color, fontsize=10)
    plt.draw()
    plt.pause(0.001)
    return abs(Resolution), ResolutionErr, popt, perr



#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#--------------------- Start Here ------------------------
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


#  PlotNoise will change the middle plot to show residuals and bottom plot to show a close up of the noise region.
PlotNoise = False
# Channel above which the code looks for the peak.
CUT = 1000
# Range of data to fit Gaussian to
LOWRANGE  = 250
HIGHRANGE = 350

errList = []
amps = []
means  = []
sigmas = []

fullfilename1, ext1 = OpenFile(1)
print fullfilename1
if 'Eu152' or 'eu152' or 'EU152' in fullfilename1:
    Source = 'Eu152'
    # Eu152 Source
    # **** Delete Entry In Energy if Peak Cannot Be Easily Fitted. ****
    # **** The Annotation will take care of itself                  ****
    Annotation = {32.0: '32.0keV', 121.78:'121.78keV', 244.6:'244.6keV', 344.3:'344.3keV', 661.7:'662keV', 778.9:'778.9keV', 964.08:'964.1keV', 1408:'1408keV'}
    #energy = [121.78, 244.6, 344.3, 661.7, 778.9, 964.08, 1408.0] FOR REFERENCE : DO NOT EDIT
    energy  = [ 121.78, 244.6, 344.3, 661.7,  778.9, 964.08, 1408.0]
elif 'NG3' or 'ng3' in fullfilename1:
    # NG3 Source
    Source = 'NG3'
    Annotation = { 59.5:'Am241', 88.0:'Cd109', 122.1:'Co57', 165.9:'Ce139', 391.7:'Hg203', 661.7:'Cs137', 898.0:'Y88(898keV)',1173.2: 'Co60 (1170keV)', 1332.5:'Co60 (1330keV)',1836.1: 'Y88(1836keV)'}
 #   energy = [ 59.5, 88.0, 122.1, 165.9, 391.7, 661.7, 898.0, 1173.2, 1332.5, 1836.1] # FOR REFERENCE: DO NOT EDIT
    energy = [  59.5, 88.0, 122.1,  661.7,  1173.2, 1332.5, 1836.1]

# Create weighting array.
# Lists the penalties for missing the data point.
w = ones(len(energy))

if ext1 == '.spe' or ext1 == '.Spe':
       data = ReadSPE(fullfilename1)
       print 'Length of input data:' , len(data)
       dataLength = len(data)

if (ext1 =='.TKA' or ext1 == '.dat'):
        data = ReadTKA(fullfilename1)
        print 'Length of input data:' , len(data)
        dataLength = len(data)

data = np.array(data)


# Create figure window with 3 subplots
# with correct spacing
fig1 = plt.figure(num=1, figsize=(9, 7), dpi=100, facecolor='w', edgecolor='k')
plt.subplot(311)

fig2 = plt.figure(num=2, figsize=(9, 5), dpi=100, facecolor='w', edgecolor='k')

figure(1)

left  = 0.125  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.05   # the bottom of the subplots of the figure
top = 0.92      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.3   # the amount of height reserved for white space between subplots

subplots_adjust(left, bottom, right, top, wspace, hspace)

# Get user to click Eu152 photopeaks.
# For each, fit a gaussian and pass back the mean value
# Save each mean value to 'means' list
ii = 0
mean_err = []


while ii  < len(energy):
    print("Please click %3.2f " % energy[ii])

    channels, data, XClick = PlotData(data)
    popt, pcov, errList = FitGaussian(XClick, channels, data)
    this_amp, this_mean, this_sigma = popt
    # Write isotope on plot next to corresponding photo-peak
    plt.text(this_mean+50, this_amp*1.1, Annotation[energy[ii]], fontsize=6)
 #  plt.text(this_mean+50, this_amp*1.1, 'Am-241 (59.5keV)', fontsize = 10)
    plt.draw()
    plt.pause(0.001)
    stddev = np.sqrt(np.diag(pcov))
    mean_err.append(np.std(errList, ddof=1))
#    mean_err.append(stddev[1])
    del errList[:]
    if (this_mean==0):
        print '\a'
    else:
        means.append(this_mean)
        ii += 1

#means.insert(0,0)
#energy.insert(0,0)
if PlotNoise is False:
    plt.subplot(312)
xquad, xlin = PlotLinearity(energy, means, w)

# New as of 9/3/17. Plot seperately, the difference between xquad and xlim

#Deviation, ch, = PlotDeviation(xquad, xlin)
residual, pct = PlotResiduals(xlin, xquad, means, energy, mean_err)

for i, val in enumerate(residual):
    print '&', float('{0:.2f}'.format(val) ) ,
print '\n'

# Return current figure to figure 1.
plt.figure(1)

# Open the 2nd file (Cs137)
fullfilename2, ext2 = OpenFile(2)
if ext2 == '.spe' or ext2 == '.Spe':
    Cs137data = ReadSPE(fullfilename2)

if (ext2 =='.TKA' or ext2 == '.dat'):
    Cs137data = ReadTKA(fullfilename2)

chan_num, chan_lin, LinCs137data, chan_cal, NormCs137data = Linearize(Cs137data, xquad, xlin)


# Plot both the reference and the target spectra to check for overlaps

fig10 = plt.figure(num=10, figsize=(8,6), facecolor='white')
plt.plot(channels, data, 'r-', channels, Cs137data, 'b')
plt.xlabel("Channels")
plt.ylabel("Counts")
plt.title("Ensure Cs137 peaks overlap.")
plt.draw()
plt.pause(0.01)

# Return current figure to figure(1)
plt.figure(1)



plt.subplot(2)
plt.errorbar(energy, residual, yerr=mean_err , fmt='rd--')
plt.axhline(0,  linewidth=1, linestyle='--', color='k')
plt.xlabel('Isotope Energy (keV)')
plt.ylabel('Deviation from $f(E)$ (keV)')

plt.subplot(313)
plt.plot(chan_lin, LinCs137data,  'b-', label='Linear Fit')
plt.plot(chan_cal, NormCs137data, 'r-', label='Quadratic Fit')
plt.xlabel('energy (keV)')
plt.xlim(0,1000)

plt.draw()
plt.pause(0.1)



Resolution, ResolutionErr, popt, perr = CalculateResolution(chan_lin, LinCs137data,'blue', 0.2*ymax)
LinResolution = 100*Resolution
LinResolutionErr = ResolutionErr
LinMean     = abs(popt[1])
LinMeanErr  = abs(perr[1])
LinSigma    = abs(popt[2])
LinSigmaErr = abs(perr[2])

Resolution, ResolutionErr, popt, pcov = CalculateResolution(chan_cal, NormCs137data,'red', 0.3*ymax)
QuadResolution = 100*Resolution
QuadResolutionErr = ResolutionErr
QuadMean     = abs(popt[1])
QuadMeanErr  = abs(perr[1])
QuadSigma    = abs(popt[2])
QuadSigmaErr = abs(perr[2])
# Write Report file.
DateStamp = time.strftime('%H:%M   %d-%m-%Y')

Title = raw_input('Enter Title for Figure: ')
plt.subplot(311)
plt.title(Title)
SaveTitle = Title.replace(" ", "_") # I f&@king love Python!
Deviations = open('Deviations.txt', 'a')
ReportFile = open(SaveTitle+'.txt', "w")
ReportFile.write("Date: %s\n " % DateStamp)
ReportFile.write("Reference Spectrum: \t %s \n" % fullfilename1)
ReportFile.write("Cs137 Spectrum: \t %s \n" % fullfilename2)
ReportFile.write("----------------------------------\n")
ReportFile.write("Linear Fit \n" )
ReportFile.write("----------------------------------\n")
ReportFile.write("Resolution:\t %3.2f%% +/- %1.2f%%   \n" % (LinResolution, LinResolutionErr) )
ReportFile.write("Mean Channel:\t %3.1f  +/-%3.2f \n" % (LinMean, LinMeanErr))
ReportFile.write("Sigma:\t \t %3.1f  +/-%3.2f \n\n" % (LinSigma, LinSigmaErr))
ReportFile.write("----------------------------------\n")
ReportFile.write("Quadratic Fit \n" )
ReportFile.write("----------------------------------\n")
ReportFile.write("Resolution:\t %3.2f%% +/- %1.2f%%   \n" % (QuadResolution, QuadResolutionErr) )
ReportFile.write("Mean Channel:\t %3.1f  +/-%3.2f \n" % (QuadMean, QuadMeanErr))
ReportFile.write("Sigma:\t \t %3.1f  +/-%3.2f \n" % (QuadSigma, QuadSigmaErr))
ReportFile.write("Residuals (keV)  \t   Residuals(%) \n")

for i , val in enumerate(residual):
    ReportFile.write(str(round(val, 2))+"\t \t"+str(round(pct[i],2))+"\n")

ReportFile.close()

Deviations = open('Deviations_Table.csv', 'a')
#for i, val in enumerate(energy):
#    Deviations.write(','+ str(val) + 'keV ,')
#Deviations.write('\n')
for i, val in enumerate(residual):
    Deviations.write(','+ str(round(val, 2)) +'keV  ')
Deviations.write('\n')
for i, val in enumerate(pct):
    Deviations.write(','+str(round(val, 2)) +'% ')
Deviations.write('\n')

if 'SaveTitle' not in locals():
    Title = raw_input('Enter Title for Figure: ')
    plt.subplot(311)
    plt.title(Title)
    SaveTitle = Title.replace(" ", "_") # I f&@king love Python!

plt.savefig(SaveTitle+'.pdf', bbox_inches='tight')
plt.savefig(SaveTitle+'.png', bbox_inches='tight')
figure(1)
plt.draw()
plt.pause(0.01)
ax1 = plt.gca()
pickle.dump(ax1 , file('Output.pickle', 'w'))

#plt.figure(2)
#plt.title(Title+' - Noise Region.')
#plt.savefig(SaveTitle+'_Noise_Region.pdf', bbox_inches='tight')
#plt.savefig(SaveTitle+'_Noise_Region.png', bbox_inches='tight')


plt.figure(3)
plt.subplot(211)
plt.title(Title+' - Deviation from Polynomial.')
plt.savefig(SaveTitle+'_Deviation.pdf', bbox_inches='tight')
plt.savefig(SaveTitle+'_Deviation.png', bbox_inches='tight')


plt.figure(1)
plt.draw()
plt.pause(0.001)

plt.figure(2)
plt.draw()
plt.pause(0.001)

print('Done');
