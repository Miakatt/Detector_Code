#!/usr/bin/env python
import sys
import os
import visa
import time

from pypeaks import Data, Intervals
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
from scipy.signal import argrelextrema

# Bin Size for histogram
BINS = 250
# Number of measurements
N = 1000
# Print out crap?
verbose = 0
#Slow it down a bit.
Wait = 0.01
# Amplitude array
amp = []

# Get it talking
rm = visa.ResourceManager()
inst = rm.open_resource("USB0::0x0957::0x17A2::MY52441236::INSTR") 

#Now some functions
def Initialize():
    ClearInstr()
    time.sleep(Wait)
    inst.timeout = 1000
    time.sleep(Wait)
    inst.term_chars = ";"
    time.sleep(Wait)
    Model = inst.query("*IDN?")
    time.sleep(5*Wait)
    print("******************************************")
    print(Model)
    print("******************************************")


def ClearInstr():
    inst.write("*CLS")

def DisplayChannel(channel, state):
    print("Channel = %s" % channel)
    inst.write(":CHANnel%s:DISplay %s" % (channel, state) )
    print("DISPLAY CHANNEL: %s" % inst.query(":CHANnel%s:DISPlay?" % channel))

def SetTriggerMode(Mode):
    inst.write("TRIGger:MODE %s" % Mode)
    print("TRIGGER MODE: %s" % inst.query(":TRIGger:MODE?"))

def SetTriggerChannel(channel):
    if (channel != "EXTernal"):
        inst.write(":TRIGger:SOURce CHANnel%s" % channel)
    if (channel == "EXTernal"):
        inst.write(":TRIGger:SOURce %s" % channel)
    source = inst.query(":TRIGger:SOURce?")
    time.sleep(Wait)
    print("TRIGGER CHANNEL: %s" % source)

def SetTimeBaseScale(timeperdiv):
    inst.write(":TIMebase:SCALe %s" % timeperdiv)
    timebase = float( inst.query(":TIMebase:SCALE?") )
    print("TIME BASE: %s s/div" % timebase)

def SetChannelOffset(channel, offset):
    inst.write(":CHANnel%s:OFFSet %s" % (channel, offset))

def SetChannelScale(channel, scale):
    inst.write(":CHANnel%s:SCALe %s" % (channel, scale))

def SetAcqType(type):
    inst.write(":ACQuire:TYPE %s" % type)

def SetAcqCount(count):
    print(":ACQuire:COUNT %s" % count)
    inst.write(":ACQuire:COUNT %s" % count)

def SetTimeBaseMode(mode):
    inst.write(":TIMebase:MODE %s" % mode)

def SetSource(channel):
    inst.write(":MEASure:SOURce CHANnel%s" % channel)
    source = inst.query(":MEASure:SOURce?")
    print("Measurement Source: %s" % source)

def Digitize(channel):
    if verbose:
        print(":DIGitize CHANnel%s" % channel)
    inst.write(":DIGitize CHANnel%s" % channel)

def GetChannelAmplitude(channel):
    inst.write(":MEASure:VAMPlitude CHANnel%d" % channel)
    amplitude = inst.query(":MEASure:VAMPlitude?")
    if (verbose):
        print("Amplitude : %s" % amplitude)

    return amplitude

def GetAmplitude():
    inst.write(":MEASure:VAMPlitude")
    amplitude = inst.query(":MEASure:VAMPlitude?")

    if (verbose):
        print("Channel 1 : %s" % amplitude)
    return amplitude

def GetChannelTop(channel):
    inst.write(":MEASure:VTOP CHANnel%s" % channel)
    top = inst.query(":MEASure:VAMPlitude?")
    if (verbose):
        print("Amplitude : %s" % top)
    
    return top

def SetMeasureWindow(window):
    inst.write(":MEASure:WINDow %s" % window)
    if (verbose):
        print("Measure Window : %s" % window)




def FitCurve(h, thebins):
    print(len(h))
    print(len(thebins))
    movingAvg = 5
    curve = np.cumsum(h, dtype=float)
    curve[movingAvg:] = (curve[movingAvg:] - curve[:-movingAvg])/movingAvg
    
    peaks = argrelextrema(curve, np.greater, 0, 3)[0]
    print(len(curve))
    print(len(thebins[:BINS]) )
    plt.plot(thebins[:BINS], curve, 'r', linewidth=3)
    
    return peaks


def PlotSlope(peaks):
    n = 4
    x = np.linspace(1,n,n)
    print(x)
    print(peaks)
    
    PeakSlope = []
    for ii in range(0,n):
        PeakSlope.append(peaks[ii])
    
    fig1 = plt.figure(figsize=(8,6))
    print(len(x))
    print(len(PeakSlope))

    plt.axis([0,8,0,110])
    plt.plot(x, PeakSlope,'o')

    p = np.polyfit(x,PeakSlope,1)
    print(p)
    y = p[0]*x+p[1]
    print(PeakSlope)
    plt.plot(x, p[0]*x + p[1] , 'r-')
    fig1.patch.set_facecolor('white')
    plt.text(1,84,'y = %1.1fx + %1.1f'% (p[0], p[1]),fontsize=24, family='serif', style='italic')
    plt.xlabel("Photon Peak", fontsize=18, family='serif')
    plt.ylabel("Bin Value", fontsize=18, family='serif')
    plt.tick_params(axis='both',labelsize=15)


def Hist(data):
    normed_value = 1
    # Make a canvas of 10x8cm
    fig = plt.figure(figsize=(8,6))
    h, thebins, patch = plt.hist(data, bins=BINS, range=[0,1], label='SiPM', normed=True)
    plt.title("SiPM Amplitudes")
    plt.xlabel("Volts", fontsize=18, family='serif')
    plt.ylabel("Frequency", fontsize=18, family='serif')
    plt.tick_params(axis='both',labelsize=15)

    fig.patch.set_facecolor('white')
    return h, thebins
def WriteToFile(data):
    Output = open('Agilent_36C.txt', 'w')
    print(Output)
    for item in data:
        Output.write("%s\n" % item)
    Output.close()

############################################################################
############################################################################

Initialize()
DisplayChannel(1, 1)
time.sleep(Wait)
DisplayChannel(4, 1)
time.sleep(Wait)
DisplayChannel(3, 0)
time.sleep(Wait)

SetTriggerMode("EDGE")
SetTriggerChannel("4")

SetTimeBaseMode("WINDow")
SetTimeBaseScale(20E-6)

SetChannelScale(1, 100E-3)
SetChannelScale(4, 1.0)
#SetChannelScale(3, 100E-3)

SetChannelOffset(1, 120E-3)
SetChannelOffset(4, 2.9)
#SetChannelOffset(3, 325E-3)
SetAcqType("NORMal")
#SetAcqCount(10)

SetMeasureWindow("ZOOM")
#SetSource(1)

############################################################################
############################################################################
print("STARTING LOOP")
#plt.axis([0,0.5,0,100])
#plt.show(block=False)
#time.sleep(2)

#oooOOOOOOooo Program Loop oooOOOOOOoooo
for ii in range(1,N):
    if (ii%100==0):
        print(ii)
    Event = int( inst.query(":TER?") )
    if (verbose):
       print("Triggered Event: %d" % Event)
    if (Event):
        #  Digitize(3)
        #        time.sleep(Wait)
        result = GetChannelTop(1)
        #  print(result)
        amp.append(float(result))

    else:
        time.sleep(Wait)

WriteToFile(amp)
h, thebins = Hist(amp)
#s = FitCurve(h,thebins)
#PlotSlope(peaks)

#plt.show()
