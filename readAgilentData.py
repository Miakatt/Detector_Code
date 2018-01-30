#!/usr/bin/env python
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
import datetime as dt
import matplotlib.dates as md
import time
import collections
from matplotlib import animation
from numpy.random import normal
from mpl_toolkits.mplot3d import Axes3D
from subprocess import call
from scipy.optimize import curve_fit
from scipy import stats
from scipy.signal import argrelextrema


amp = []
BINS = 250

def ReadFile(inputfile):
    print ("Reading FileL %s" % inputfile)
    with open(inputfile,'r') as f:
        for line in f:
            data = line.split()
            amp.append(float(data[0]))
#    print(amp)
    return amp
        
def PlotIt(amp):
    fig = plt.figure(figsize=(10,8))
    
    h, thebins, patch = plt.hist(amp, bins=BINS, range=[0,1], normed=True)

    plt.title("Photon Peaks")
    plt.xlabel("Amplitude (V)", fontsize=18, family='serif')
    plt.ylabel("Frequency", fontsize=18, family='serif')
    plt.tick_params(axis='both',labelsize=15)
    fig.patch.set_facecolor('white')
    

    return h, thebins

def FitCurve(h, thebins):
    movingAvg = 5
    curve = np.cumsum(h, dtype=float)
    curve[movingAvg:] = (curve[movingAvg:] - curve[:-movingAvg])/movingAvg
    
    peaks = argrelextrema(curve, np.greater, 0, 3)[0]

    plt.plot(thebins[:BINS], curve, 'r', linewidth=3)
#    plt.show()

    return peaks

def PlotSlope(peaks):
    n = 4
    x = np.linspace(1,n,n)
    print(x)

    PeakSlope = []
    for ii in range(1,n+1):
        PeakSlope.append(peaks[ii])
   
    fig1 = plt.figure(figsize=(8,6))
    plt.axis([0,8,0,110])
    plt.plot(x, PeakSlope,'o')

    p = np.polyfit(x,PeakSlope,1)
    print(p)
    y = p[0]*x+p[1]
    print(PeakSlope)
    plt.plot(x, p[0]*x + p[1] , 'r-')
    fig1.patch.set_facecolor('white')
    plt.text(1,84,'y = %1.1fx + %1.1f'% (p[0], p[1]),fontsize=24, family='serif', style='italic')
    plt.text(1,94,'%sC'% Temperature,fontsize=24, family='serif', style='italic')

    plt.xlabel("Photon Peak", fontsize=18, family='serif')
    plt.ylabel("Bin Value", fontsize=18, family='serif')
    plt.tick_params(axis='both',labelsize=15)


inputfile = sys.argv[1]
Temperature = sys.argv[2]
ReadFile(inputfile)
h, thebins = PlotIt(amp)
peaks = FitCurve(h, thebins)
print(peaks)
PlotSlope(peaks)

plt.show()
