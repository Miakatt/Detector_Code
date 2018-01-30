#!/usr/bin/env python

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
                # Get number of data points from file
                dataLength = int(datum[2:])
            # Read those data points from file
            if (x>11 and x<dataLength+12):
                
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
        for i in range(50):
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

    # Temporary XQuad
#    xquad = np.poly1d([2.95e-5, 0.55, -74.1575])    

    xlin  = FitLinear(energy, means, w)
    
    coefs = xquad.c
    print coefs
    x = range(dataLength)
    
    # Plot both fits to the mean points
    #fig1 = figure(facecolor='white')
    if Source=='Eu152':
        plt.plot(means, energy, 'ro',  label='Eu152 Peaks')
    elif Source=='NG3':
        plt.plot(means, energy, 'ro',  label='NG3 Peaks')
    
    plt.plot(x, xquad(x), 'r-',   label='Quadratic Fit')
    plt.plot(x, xlin(x),'b-.',     label='Linear Fit')

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
    Deviation = np.zeros(len(ch))    
    np.divide(( QuadFit - LinFit ), LinFit, Deviation ) 
    
    return Deviation, ch

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def Linearize(Cs137data, xquad, xlin):
    
    chan_num = range(len(Cs137data))                #[0, 1, 2 .... dataLength]
    chan_center = [bin+0.5 for bin in chan_num]     #[0.5, 1.5, 2.5 .... dataLength + 0.5]
    chan_center2= [bin+1.0 for bin in chan_center]  #[1.5, 2.5, 3.5 .... dataLength + 1.5]



    # Python equivalent of a*x[i]*x[i] + b*x[i] * c  within a for loop.
    # Quadratic fit normalisation
    chan_cal  = xquad(chan_center)
    chan_cal2 = xquad(chan_center2)
    # Linear fit normalisation
    chan_lin   = xlin(chan_center)
    chan_lin2  = xlin(chan_center2)

    chan_quad_width = np.array(chan_cal2) -  np.array(chan_cal) # quadratic fit bin width
    chan_lin_width  = np.array(chan_lin2) -  np.array(chan_lin) # linear fit bin width
    
    NormCs137data = np.array(Cs137data) / np.array(chan_quad_width)
    LinCs137data  = np.array(Cs137data) / np.array(chan_lin_width)

    return chan_num, chan_lin, LinCs137data, chan_cal, NormCs137data

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def CalculateResolution(chan, Cs137data, color, placement, means):
    # Find the channel (>Ch 500) that hold the maximum data value. [0][0] takes it out of the list structure.
    #meanCh = np.where(Cs137data==max(Cs137data[1200:]))[0][0]
    # Get the mean of the photopeak (each channel does not equal one increment on 'x-axis'
    for i in range(len(means)):
        PhotoPeakMean = chan[means[i]]
        PhotoPeakAmp  = Cs137data[means[i]]
        print 'Mean: ', PhotoPeakMean 
        print 'Amp: ', PhotoPeakAmp

        sig = sqrt(PhotoPeakMean)
        slope = -0.001
        intercept = 5
        p0=[PhotoPeakAmp, PhotoPeakMean, sig, slope, intercept]
        # Fit to a range of data around the mean
        #if means[i] <= 665:
        #    halfRange = means[i]*0.12
        #elif means[i] > 665:
        #    halfRange = means[i]*0.075
        halfRange = 2*np.sqrt(means[i])
        fitRange = chan[means[i]-halfRange:means[i]+halfRange]
        fitData  = Cs137data[means[i]-halfRange:means[i]+halfRange]

        # Fit the photopeak with a linear background
        popt, pcov = curve_fit(gausWithBackground, fitRange, fitData, p0=p0)

        #popt, pcov = curve_fit(gaus, fitRange, fitData, p0=[amp, mean, sig])
        plt.plot(fitRange, gausWithBackground(fitRange, *popt), color, linewidth=2)
    
        #plt.plot(fitRange, gaus(fitRange, *popt), color, linewidth=2)
        plt.draw()
        plt.pause(0.01)
        figure(1)


    
         # Calculate resolution and error
        perr = np.sqrt(np.diag(pcov))
      #  print perr
      #  print popt
        Resolution.append( abs(2.35*popt[2]/popt[1]) )
        # Get fractional error of mean and FWHM
        MeanErr = perr[1]/popt[1]
        FWHMErr = perr[2]/popt[2]

        ResolutionErr.append(Resolution * sqrt(FWHMErr**2 + MeanErr**2))
        print ResolutionErr[i]

        #  plt.plot(chan_lin, LinCs137data, 'blue', label='Corrected')
        if color=='blue':
            plt.plot(fitRange, gausWithBackground(fitRange, *popt), color, linewidth=2)

        if color=='black':
           plt.plot(fitRange, gausWithBackground(fitRange, *popt), color, linewidth=1, label='{0:.1f}'.format(abs(Resolution[i])*100)+'% @{0:.1f}'.format(energy[i])+'keV' )   
        #plt.plot(fitRange, gaus(fitRange, *popt), color, linewidth=2, label='Gaussian Fit')
        

        # annotate plot with resolution.
     #   if color=='red':
     #       plt.text(PhotoPeakMean, PhotoPeakAmp-placement, r'ER: %2.2f%% +/- %1.2f%%' % (abs(100*Resolution) , abs(100*ResolutionErr)), color=color, fontsize=5)

    plt.semilogy()    
  #  plt.legend(loc='upper right', fontsize=8)    
    ymin, ymax = plt.ylim()
    xmin, xmax = plt.xlim()        
    ylim(-10, ymax)
    xlim(-10, 2000)
    xlabel('Energy (keV)', fontsize=10)
    ylabel('Counts', fontsize=10)
    
    plt.plot()
    plt.draw()

    return Resolution, ResolutionErr, popt, perr

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

def PlotResiduals(xlin, xquad, means, energy, mean_err):
    # Evaluate the quadratic at each mean. Calculate the difference
    # between the energy  and the quadratic
    residual = []
    for i, mean in enumerate(means):
        residual.append( energy[i] - xquad(mean) )

    print energy
    print residual
    print mean_err   
    print float(max(residual)), float(min(residual) )
    plt.figure(num=2, figsize=(9, 7), dpi=100, facecolor='w', edgecolor='k')
    plt.errorbar(energy, residual, yerr=mean_err , fmt='rd--')
    plt.axhline(0,  linewidth=1, linestyle='--', color='k')
    plt.axhline(float(max(residual)), linewidth=1, linestyle='.', color='k')
    plt.xlabel('Isotope Energy (keV)')
    plt.ylabel('Deviation from $f(E)$ (keV)' ) 

    plt.draw()
    plt.pause(0.001)

    return residual

#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#--------------------- Start Here ------------------------
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# Create weighting array to force fit through 0, 0
# Lists the penalties for missing the data point.
# Channel above which the code looks for the peak.
CUT = 1000
# Range of data to fit Gaussian to
LOWRANGE  = 250
HIGHRANGE = 350

errList = []
amps = []
means  = []
#means.insert(0,0)
sigmas = []

#Open first file (Eu152 spectrum)
fullfilename1, ext1 = OpenFile(1)
if ext1 == '.spe' or ext1 == '.Spe':
       data = ReadSPE(fullfilename1)
       print 'Length of input data:' , len(data)
       dataLength = len(data)
    
if (ext1 =='.TKA' or ext1 == '.dat'):
        data = ReadTKA(fullfilename1)
        print 'Length of input data:' , len(data)
        dataLength = len(data)

data = np.array(data)


if 'Eu152' in fullfilename1:
    Source = 'Eu152'
    # Eu152 Source
    # **** Delete Entry In Energy if Peak Cannot Be Easily Fitted. ****
    # **** The Anotation will take care of itself                  ****
    Comps = {121.78:'121.78keV', 244.6:'244.6keV', 344.3:'344.3keV', 661.7:'662keV', 778.9:'778.9keV', 964.08:'964.1keV', 1408:'1408keV'}
    #energy = [121.78, 244.7, 344.3, 661.7, 778.9, 964.08, 1408.0] FOR REFERENCE : DO NOT EDIT
    energy = [121.78, 244.6, 344.3 , 778.9, 964.08, 1408.0]
elif 'NG3' in fullfilename1:
    # NG3 Source
    Source = 'NG3'
    Comps = { 59.5:'Am241', 88.0:'Cd109', 122.1:'Co57', 165.9:'Ce139', 391.7:'Hg203', 661.7:'Cs137', 898.0:'Y88(898keV)',1173.2: 'Co60 (1170keV)', 1332.5:'Co60 (1330keV)',1836.1: 'Y88(1836keV)'}
#    energy = [ 59.5, 88.0, 122.1, 165.9, 391.7, 661.7, 898.0, 1173.2, 1332.5, 1836.1] FOR REFERENCE: DO NOT EDIT
    energy = [  59.5, 88.0, 122.1, 165.9, 391.7, 661.7, 898.0, 1173.2, 1332.5, 1836.1]    
# Create weighting array to force fit through 0, 0
# Lists the penalties for missing the data point.
w = ones(len(energy))
w[0] = 0

def FitGaussian(XClick, channels, data):
    mean = XClick[0][0]
    amp  = XClick[0][1]
    popt = []
    sigma = np.sqrt(mean)

    try:
        for i in range(50):
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

#dataLength = 4096

# Create figure window with 3 subplots
# with correct spacing
plt.figure(num=1, figsize=(9, 7), dpi=100, facecolor='w', edgecolor='k')
plt.subplot(311)

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
    plt.text(this_mean+50, this_amp*1.1, Comps[energy[ii]], fontsize=6)
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

plt.subplot(312)
xquad, xlin = PlotLinearity(energy, means, w)


# New as of 6/4/17. Plot seperately, the difference between xquad and xlim
#Deviation, ch, = PlotDeviation(xquad, xlin)
residual = PlotResiduals(xlin, xquad, means, energy, mean_err)

# Return current figure to figure 1. 
plt.figure(1)
plt.subplot(311)

# Open the 2nd file (Cs137)
fullfilename2, ext2 = OpenFile(1)
if ext2 == '.spe' or ext2 == '.Spe':
    Cs137data = ReadSPE(fullfilename2)

if (ext2 =='.TKA' or ext2 == '.dat'):
    Cs137data = ReadTKA(fullfilename2)

chan_num, chan_lin, LinCs137data, chan_cal, NormCs137data = Linearize(Cs137data, xquad, xlin)

plt.subplot(313)
plt.plot(chan_lin, LinCs137data,'b--', label='Linear Fit')
plt.plot(chan_cal, NormCs137data, 'r-.', label='Quadratic Fit')
#plt.xlim(0,1000)
plt.pause(0.001)
ymin, ymax = plt.ylim()


Resolution = []
ResolutionErr = []

# CalculateResolution(mean, amplitude, x, y, line color, text y position below ymax )
Resolution, ResolutionErr, popt, perr = CalculateResolution(chan_lin, LinCs137data,'blue', 0.2*ymax, means)
LinResolution = (100*Resolution)
LinResolutionErr = (ResolutionErr)

LinMean     = abs(popt[1])
LinMeanErr  = abs(perr[1])
LinSigma    = abs(popt[2])
LinSigmaErr = abs(perr[2])


Resolution, ResolutionErr, popt, pcov = CalculateResolution(chan_cal, NormCs137data,'black', 0.3*ymax, means)
PeakMeans = []
for i, val in enumerate(means):
    PeakMeans.append(chan_cal[val])

QuadResolution = 100*Resolution
QuadResolutionErr = ResolutionErr

QuadMean     = abs(popt[1])
QuadMeanErr  = abs(perr[1])
QuadSigma    = abs(popt[2])
QuadSigmaErr = abs(perr[2])

Title = raw_input('Enter Title for Figure: ')
Crystal     = raw_input('Enter Crystal Type: ')
Temperature = raw_input('Enter Temperature: ')
#Voltage     = raw_input('Enter Bias Voltage: ')

ReportFile = open(Crystal+'.txt', "a")
DeviationFile = open(Crystal+'_Deviation.txt', "a")
ReportFile.write("%d,   " % ( float(Temperature) ) )
coefs = xquad.c
for r in range(len(means)):
    ReportFile.write(", %0.1f , %1.1f " % (means[r], 100*Resolution[r]) )
#    ReportFile.write(" %1.1f,   %2.1f," %  (chan_cal[means[r]], mean_err[r])  )
ReportFile.write(" %2.2e, %2.2f, %2.2f" % (coefs[0], coefs[1], coefs[2]))        
ReportFile.write("\n")
ReportFile.close()

DeviationFile.write("%d,   " % ( float(Temperature) ) )
for r in range(len(residual)):
    DeviationFile.write(" %1.1f, %2.2f," %  (residual[r], mean_err[r]) )    
DeviationFile.write("\n")



plt.subplot(311)
plt.title(Title)
SaveTitle = Title.replace(" ", "_") # I f&@king love Python!
plt.savefig(SaveTitle+'.pdf', bbox_inches='tight')
plt.draw()
plt.pause(0.01)

plt.figure(2)
plt.title(Title+' - Deviation from Polynomial.')
plt.savefig(SaveTitle+'_Deviation.pdf', bbox_inches='tight')


print('Done');

