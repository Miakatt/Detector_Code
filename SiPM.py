#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on 19-6-17

@author: Alan J Bell.
SiPM / Scintillator Model for D3S.
Useage: SiPM.py <config file>
Generates output text file for reading with LTSPICE Model 'SiPMs + pre-amp + peak detect - offset corrected'
Running LTSpice simulation, and probing the output produces an output text file of the peak hold waveforms.
This text file can be read with 'PeakPlotting.py' or 'AllPeakPlotting.py' to get the channel numbers.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import ConfigParser
import random
import math
import json
from termcolor import colored
from scipy.interpolate import interp1d
from scipy.stats import poisson
import csv

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def ReadDarkCurrentTable(DarkCurrentCSV):
    # Simply read the dark current from the csv file in to a list (array)
    DarkCurrentCSV = list(csv.reader(open(DarkCurrentCSV)))

    return DarkCurrentCSV

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def ReadEmissionSpectrum(EmissionSpectrumCSV):
    # Read CSV file of Emission Intensities in to a dictionary
    # required for calculating overall detection efficiency with PDE_WLEff.
    # PDE_WL values are used as the key in the EmissionWL dictionary to extract the emission intensity
    EmissionDict = {}
    with open(EmissionSpectrumCSV, 'r') as csvfile:
        read = csv.reader(csvfile, delimiter=',')
        # Skip Header Info
        for i in range(4):
            next(read)
        try:
            for row in read:
                EmissionDict.update({int(row[0]): float(row[1])})
        except(IndexError):
            print 'End of file...'

    return EmissionDict

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def DictToList(EmissionDict):
    # Extract dictionary keys and values in to 2 lists for plotting.c
    Wavelength = []
    Intensity = []

    for key, value in EmissionDict.iteritems():
        Wavelength.append(key)
        Intensity.append(value)

    return Wavelength, Intensity

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def ReadPDEvsWavelength(PDE_WavelengthCSV):
    PDE_WL = []
    PDE_WLEff = []
    with open(PDE_WavelengthCSV, 'r') as csvfile:

        read = csv.reader(csvfile, delimiter=',')
        # Skip Header Info
        for i in range(4):
            next(read)
        try:
            for row in read:
                PDE_WL.append(int(row[0]))
                PDE_WLEff.append(float(row[1]))
        except(IndexError):
            print 'End of file...'

    return PDE_WL, PDE_WLEff

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def CombineEmissionAndDetectionEff(PDE_WL, PDE_WLEff, EmissionDict):
    CombinedEff = []
    for w, val in enumerate(PDE_WL):
        ThisPDEEff = PDE_WLEff[w]
        try:
            # Get the value whose key is the PDE Wavelength value (if it exists)
            ThisEmissionEff = EmissionDict[val]

        except (KeyError):
            # if the key does not exist, set to value = zero
            ThisEmissionEff = 0

        CombinedEff.append((ThisEmissionEff * ThisPDEEff))
    return CombinedEff

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def EnergyToPhotons(E, PhotonsPerkeV):
    # Take energy input as keV. Convert to No of Photons with
    # the non-proportionality function for CsI. (This is for a 3us shaping time. 12us isn't so large):
    # f(E) = 1 +C*(1-1/(1+exp(-k*(E-m))))
    C = 1.74317
    K = 0.003895
    m = -626.578

    f = 1 + C*(1-(1/(1+np.exp(-K*(E-m)))))
    print 'Low-Energy Correction Factor: ' , f

    NPhotonsGenerated = int((E * PhotonsPerkeV * f)/6)


    return NPhotonsGenerated

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def GetDoubleHitProbability(k, NoOfPixels, NPhotons):
    mu = NPhotons/NoOfPixels
    DoubleHitProbability = 0
    for i, val in enumerate(k):
        DoubleHitProbability += np.exp(-mu) * mu**val / math.factorial(val)

    return DoubleHitProbability*NoOfPixels





# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def PhotonsCollected(PDE_WL, PDE_WLEff, PDE_OV, PDE_OVEff, CombinedEff, OverVoltage, NPhotonsAtSiPM):
    # CombinedEff is the combination of SiPM PDE(WL) and Xtal emission spectra.
    AverageEff = sum(CombinedEff)/sum(x > 0 for x in CombinedEff)
    fOV = interp1d(PDE_OV, PDE_OVEff)
    # PDE change due to over voltage for each SiPM, relative to reference PDE (2.5V Overvoltage @21C).
    PDE_OverVolts = []
    PDE_OverVolts[:] = [OV/fOV(2.5) for OV in fOV(OverVoltage)]
    print 'PDE ratio due to Voltage : ',
    print [round(P,3) for P in PDE_OverVolts[:]]
    NPhotonsCollected = 0
    # Get Photons collected with each PDE due to each SiPMs overvoltage. Assume
    # all PhotonsAtSiPM are equally dispersed over all NoOfSiPMs.
    NPhotonsCollected = [int((NPhotonsAtSiPM * P * AverageEff) /NoOfSiPMs) for P in PDE_OverVolts[:]]
    # Sum them to get total Photons Collected from N sipms.
    return NPhotonsCollected, PDE_OverVolts, AverageEff


# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def CrystalLosses(NPhotons, CrystalCrossSection, SiPMArea):
    # Account for the photons that his the SiPM facing surface but miss the SiPMs
    print SiPMArea, CrystalCrossSection

    if (SiPMArea/CrystalCrossSection)>1:
        EffectiveArea = 1
    else:
        EffectiveArea = SiPMArea/CrystalCrossSection
    CrystalLoss = REFLECTIVE_FACTOR * EffectiveArea
    print EffectiveArea
    print 'Geometric Acceptance: %.1f%%' % (100*CrystalLoss)

    NPhotonsAtSiPM = int(NPhotons*(CrystalLoss))
    return NPhotonsAtSiPM

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

#def RejectDoubleHits(NPhotonsAtSiPM, NoOfPixels):
#    NFire = NoOfPixels*(1-np.exp(NPhotonsAtSiPM/NoOfPixels))
#    print NFire
#    return NFire

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def GetCrossTalkPercentage(CrossTalkvsOV, OverVoltage):
    CrossTalk = []
    CrossTalk[:] = [CrossTalkvsOV * OV for OV in OverVoltage]
    return CrossTalk

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def GetDarkCurrent(NoOfSiPMs, DarkCurrentTable, Bias, Temperature):
    # Find index (Row) equal to Temperature
    List =  zip(*DarkCurrentTable)[0][1:]
    ClosestTemperature = min(List, key=lambda x:abs(float(x)-Temperature))
    Row = List.index(ClosestTemperature)+1

    List = DarkCurrentTable[0][1:]
    ClosestBias = min(List, key=lambda x:abs(float(x)-float(Bias)))
    Col = List.index(ClosestBias)+1


    DarkCurrent = NoOfSiPMs * float(DarkCurrentTable[Row][Col])
    return DarkCurrent

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def GetDarkRate(DarkCurrent, Gain):
    DCR = 0
    DCR = sum([(1e-9*DarkCurrent/(G * q)) for G in Gain])

    return DCR

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def GetFiringPixels(NPhotonsCollected, NoOfPixels, PDE_OverVolts, PDE_Wavelength, CrossTalk):
    NPC = NPhotonsCollected
    NP  = NoOfPixels
    #CT  = CrossTalk

    FiringPixelsIncCrossTalk = []
    FiringPixels = []
    for i, CT in enumerate(CrossTalk):
        FiringPixelsIncCrossTalk.append(int(NP*(1-np.exp((-NPC[i]*(1+(CT)))/NP))) )  # ??? Still required?

    for i, FPIC in enumerate(FiringPixelsIncCrossTalk):
        FiringPixels.append(int(FPIC/(1+(CrossTalk[i]))))

    print 'Pixels Firing (Inc XTalk) : ', FiringPixelsIncCrossTalk
    print 'Total Pixels Firing (Inc XTalk) : ', sum(FiringPixelsIncCrossTalk)
    print 'Total Pixels Firing             : ', sum(FiringPixels)

    return sum(FiringPixelsIncCrossTalk), sum(FiringPixels)
# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
def GetIntrinsicResolution(E, IntRes_vs_Energy, IntResolutionEnergy):
    fIR = interp1d(IntResolutionEnergy, IntRes_vs_Energy)
    IntrinsicResolution = fIR(E)

    if (PlotData):
        fig97 = plt.figure(figsize=(8,6), facecolor='white')
        IR = np.arange(min(IntResolutionEnergy), max(IntResolutionEnergy))
        plt.plot(IntResolutionEnergy, IntRes_vs_Energy, 'bs')
        plt.plot(IntResolutionEnergy, fIR(IntResolutionEnergy) )
        plt.semilogx()
        plt.xlabel("Energy(keV)")
        plt.ylabel("Intrisic Resolution (%) [A.Syntfeld-Kazuch et. al]")
        plt.grid(True)
        plt.draw()
        plt.pause(0.01)
    return IntrinsicResolution



def GetPhotonCountingResolution(GainVariance, FiringPixels):
    Res_Photons = (100*2.35*np.sqrt((1+GainVariance)/FiringPixels))
    return Res_Photons

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def GetNoiseResolution(DarkRate, IntegrationTime, FiringPixels):
    Res_Noise = 100*(2.35*np.sqrt(DarkRate*IntegrationTime))/(FiringPixels)
    return Res_Noise

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def GetResolution(PhotonResolution, NoiseResolution, IntrinsicResolution):
    #print colored(str(PhotonResolution)+'  '+str(NoiseResolution)+'  '+str(IntrinsicResolution),'red')
    Resolution = (np.sqrt(PhotonResolution**2 + NoiseResolution**2 + IntrinsicResolution**2))
    return Resolution

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def GetGainDueToTemperature(GainTemperature, GainvsTemp, Temperature):
    # Calculate the overall gain due to temperature
    fGTemp = interp1d(GainTemperature, GainvsTemp)
    GainAtThisTemp = fGTemp(Temperature)
    GainChange = int(GainAtThisTemp - fGTemp(20))

    if (PlotData):
        plt.figure(figsize=(8,6), facecolor='white')
        plt.plot(GainTemperature, fGTemp(GainTemperature), 'ro-')
        plt.title('Temperature Gain Dependence')
        plt.xlabel('Temperature (C)')
        plt.ylabel('Gain')
        ax = plt.gca()
        plt.semilogy()
        ax.yaxis.grid(True, which="both",ls="--")
        ax.grid(True)
        plt.draw()
        plt.pause(0.01)
    print 'Gain at this temp (2.5V OverVoltage) ', GainAtThisTemp
    return GainAtThisTemp, GainChange

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def GetGainDueToOvervoltage(GainOV, GainvsOV, OverVoltage):
    # Calculate the overall gain due to overvoltage
    # requires GainAtThisTemp first, which represents the gain at 2.5V overvoltage
    GainAtThisTemp, GainChange = GetGainDueToTemperature(GainTemperature, GainvsTemp, Temperature)
    # Modify Overvoltage/Gain curve to incorporate temperature gain change.
    ModGainvsOV = []
    for N in GainvsOV:
        ModGainvsOV.append(N+GainChange)

    # Interpolate the point-wise data
    fGOV = interp1d(GainOV, ModGainvsOV)
    GainAtThisBias = fGOV(OverVoltage)
    if (PlotData):
        plt.figure(figsize=(8,6), facecolor='white')
        plt.plot(GainOV, fGOV(GainOV), 'ro-')
        plt.title('Overvoltage Gain')
        plt.xlabel('Overvoltage (V)')
        plt.ylabel('Gain')
        ax = plt.gca()
        plt.semilogy()
        ax.yaxis.grid(True, which="both",ls="--")
        plt.pause(0.01)
    print 'Gain at this overbias ', GainAtThisBias
    return GainAtThisBias

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def PhotonsToCharge(FiringPixelsIncCrossTalk, Capacitance_nF, IntegrationTime, Gain):
    # Convert number of photons to charge.
    #print 'a ', FiringPixelsIncCrossTalk
    # print 'b ', DarkRate
    # print 'c ', Gain
    # print 'd ', IntegrationTime
    # print 'e ', OverVoltage
    Charge_nC = 0
    for G in Gain:
        Charge_nC += (1e9*((FiringPixelsIncCrossTalk/NoOfSiPMs) * G * q))

    Current_uA = 1e-3*Charge_nC/IntegrationTime  # 1e-3 = 1e-9 * 1e6 conversion from nanoCoulombs to microAmps

    return Charge_nC, Current_uA

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def ChargeToPulse(Q, Capacitance, Rq, Rs, tauList, SampleStart, SampleEnd, SampleStep):
    tau_scint1, tau1_Intensity, tau_scint2, tau2_Intensity, SiPMtauRise, SiPMtauFall = TauList

    ScintillationTau = (tau1_Intensity*tau_scint1 + tau2_Intensity*tau_scint2 - tau_rise )
    Imax =  Charge_nC/ScintillationTau
    print 'IMax ', Imax

    # Calculate the scintillation decay curve as 1ns steps
    I_scintArray = []
    I_scint = 0
    I_scintRise = 0
    # DEFINE HOW THE PIXELS FIRE AND THE CHARGE ARRIVES DUE TO SCINTILLATION
    T1 = range(SampleStart, SampleEnd+1, SampleStep )

    for x, val in enumerate(T1):
        I_scintRise  = Imax * (1-np.exp(-val/tau_rise))
        I_scintFall  = I_scintRise * ((tau1_Intensity * (np.exp(-val/tau_scint1)) + tau2_Intensity * (np.exp(-val/tau_scint2)))) # V = Qt/Capacitance
        I_scintArray.append(I_scintFall * 1000) # CONVERT TO mA

    if (0):
        plt.plot(T1, I_scintArray, 'r')
        plt.draw()

    # SiPMRise = 0
    # SiPMCurrent= []
    # SiPMFall = []
    # SiPMVoltage = []
    # T = range(SampleStart, SampleEnd-SampleStep+1, SampleStep)
    # for x, val in enumerate(T):
    #     # CURRENT FROM SIPM
    #     SiPMRise = (I_scintArray[val]) * (1-np.exp(-val/SiPMtauRise))
    #     SiPMCurrent.append(SiPMRise)
    # #    SiPMFall = SiPMRise *  (1/np.exp(val/SiPMtauFall))
    #     SiPMVoltage.append(SiPMRise*OutputLoad)

    return I_scintArray, T1



def GetNoise(NoiseOffset, NoiseLevel, NoOfSamples):

    noise = np.random.normal(NoiseOffset, NoiseLevel, NoOfSamples)
    return noise

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def PlotWaveforms(ScintSignal, T):
    if 'fig99' not in locals():
        fig99 = plt.figure(num=99, figsize=(9,6), facecolor='white')

    from matplotlib.ticker import MultipleLocator
    plt.hold(True)
    #plt.clf()
    # Convert ns time axis to microsecond
    t_us = [ns*1e-3 for ns in T]
    # Convert V to mV
    Sig  = [v for v in ScintSignal]
    ax = plt.gca()
    plt.axhline(0, linestyle='-.', linewidth=2, color='k')
    plt.plot(t_us, Sig, linewidth=2)
    mlx = MultipleLocator(1)
    mly = MultipleLocator(10)
    plt.axes().yaxis.set_minor_locator(mly)
    plt.axes().xaxis.set_minor_locator(mlx)
    plt.xlabel("Time (us)")
    plt.ylabel("Scintillator Current (mA)")
    plt.title("Scintillator Charge Pulse")
    plt.grid()
    plt.draw()
    plt.pause(1)

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

def WriteWaveformToTextFile(E, AcqTime, SampleStep,  ScintSignal, Period_us):
    Filename = sys.argv[1]
    #remove filename extention.
    Name = os.path.splitext(Filename)[0]
    WaveformSave = Name+'_WaveformData_'+str(E)+'keV_'+str(NoOfEvents)+'evts_'+str(IdenticalRepeats)+'Reps_'+str(Period_us/1000)+'us_'+str(SampleStep)+'ns.txt'

    if os.path.isfile(WaveformSave):
        # delete waveform output file if it already exists
        print 'Removing previous Waveform file...'
        os.remove(WaveformSave)

    Waveform = open(WaveformSave, "a")

    # Ensure last Vo data point is zero (so we're not feeding in constant, non-zero charge during the gap)
    ScintSignal[-1] = 0.0
    for i in range(0, IdenticalRepeats):
        for x in range(len(ScintSignal)):
            Waveform.write( ('%dn'+'\t'+'%.3fm\n') % (AcqTime, ScintSignal[x]) )
            AcqTime += SampleStep

        AcqTime = (i+1)*Period_us
    Waveform.close()

    return AcqTime

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#                             CODE STARTS HERE
# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
AcqTime = 0
#REFLECTIVE_FACTOR = 1.68
q = 1.602e-19 # electron charge
PulseLength = 500
config      = ConfigParser.RawConfigParser()
configFilePath = sys.argv[1]
print colored('Using '+configFilePath,'blue')
config.read(configFilePath)
PlotData = config.getboolean("Sim_Settings","PlotData")
PlotWaveform = config.getboolean("Sim_Settings","PlotWaveform")
Verbose  = config.getboolean("Sim_Settings","Verbose")
NoOfEvents = config.getint("Sim_Settings","NoOfEvents")
# ------------------------------------------------------------------------------
Temperature         = config.getfloat("Environment", "Temperature")
# ------------------------------------------------------------------------------
Energy              = json.loads(config.get("Source", "Energy"))
# ------------------------------------------------------------------------------
PhotonsPerkeV       = int(config.get("Xtal", "PhotonsPerkeV"))
CrystalCrossSection = config.getfloat("Xtal", "CrystalCrossSection")
IntRes_vs_Energy = json.loads(config.get("Xtal", "IntRes_vs_Energy"))
IntResolutionEnergy = json.loads(config.get("Xtal", "IntResolutionEnergy"))
EmissionSpectrumCSV = config.get("Xtal","EmissionSpectrum")
tauRiseIntensity    = config.getfloat("Xtal", "TauRiseIntensity")
tau_rise            = config.getfloat("Xtal", "Tau_Rise")
tau_scint1          = config.getfloat("Xtal","Tau_Scint1")
tau1_Intensity      = config.getfloat("Xtal","Tau1_Intensity")
tau_scint2          = config.getfloat("Xtal","Tau_Scint2")
tau2_Intensity      = config.getfloat("Xtal","Tau2_Intensity")
REFLECTIVE_FACTOR   = config.getfloat("Xtal","REFLECTIVE_FACTOR")

# ------------------------------------------------------------------------------
SiPMCapacitance_nF    = config.getfloat("SiPM", "SiPMCapacitance_nF")
Rs                = config.getfloat("SiPM", "Rs")
Rq                = config.getfloat("SiPM", "Rq")
TempDependence    = config.getfloat("SiPM", "TempDependence")
Breakdown21C      = json.loads(config.get("SiPM","Breakdown21C"))
NoOfSiPMs         = config.getint("SiPM", "NoOfSiPMs")
SiPMSize          = config.getint("SiPM", "SiPMSize")
NoOfPixels        = config.getint("SiPM", "NoOfPixels")
#FillFactor        = config.getfloat("SiPM", "FillFactor")
NoOfSiPMs         = config.getint("SiPM", "NoOfSiPMs")
DarkCurrentCSV    = config.get("SiPM", "DarkCurrentCSV")
CrossTalkvsOV     = config.getfloat("SiPM", "CrossTalkvsOV")
GainVariance      = config.getfloat("SiPM", "GainVariance")
RecoveryTimens    = config.getfloat("SiPM", "RecoveryTimens")
GainOV            = json.loads(config.get("SiPM", "GainOV"))
GainvsOV          = json.loads(config.get("SiPM", "GainvsOV"))
GainTemperature   = json.loads(config.get("SiPM", "GainTemperature"))
GainvsTemp        = json.loads(config.get("SiPM", "GainvsTemp"))
PDE_OV            = json.loads(config.get("SiPM","Overvoltage"))
PDE_OVEff         = json.loads(config.get("SiPM","PDE_Overvoltage"))
PDE_WavelengthCSV = config.get("SiPM","PDE_Wavelength")

# ------------------------------------------------------------------------------
Bias              = config.getfloat("Electronics","Bias")
#InputImp          = config.getfloat("Electronics","InputImp")
OutputLoad        = config.getfloat("Electronics", "OutputLoad")
IntegrationTime   = config.getfloat("Electronics", "IntegrationTime")
NoiseOffset       = config.getfloat("Electronics", "NoiseOffset")
NoiseLevel        = config.getfloat("Electronics", "NoiseLevel")
# ------------------------------------------------------------------------------
SampleStart       = config.getint("Samples","SampleStart")
SampleEnd         = config.getint("Samples","SampleEnd")
SampleStep        = config.getint("Samples","SampleStep")
IdenticalRepeats  = config.getint("Samples","IdenticalRepeats")
Period_us         = config.getint("Samples","Period_us")*1000
# ------------------------------------------------------------------------------
if len(Breakdown21C) != NoOfSiPMs:
    print colored('---------------------------------------------------------------------------------------------------------', 'red')
    print colored('Warning: You have provided more or less breakdown values than there are SiPMs. Please amend config file! ' ,'red')
    print colored('---------------------------------------------------------------------------------------------------------', 'red')
    sys.exit()
OverVoltage = []
OverVoltage[:] = [round(Bias-BD,2) for BD in Breakdown21C]
# THE FOLLOWING ARE NOT USED IN THE FUNCTIONS. INSTEAD, THESE PARAMETERS ARE MODELLED
# IN LTSPICE
Capacitance_nF = SiPMCapacitance_nF * NoOfSiPMs
SiPMtauRise  = (Rs)*Capacitance_nF
SiPMtauFall  = (Rq)*Capacitance_nF

NoOfSamples = (SampleEnd-SampleStart)/SampleStep


PDE_WL, PDE_WLEff = ReadPDEvsWavelength(PDE_WavelengthCSV)
EmissionDict = ReadEmissionSpectrum(EmissionSpectrumCSV)

CombinedEff = CombineEmissionAndDetectionEff(PDE_WL, PDE_WLEff, EmissionDict)
if (PlotData):
    Wavelength, Intensity = DictToList(EmissionDict)
    plt.figure(figsize=(8,6), facecolor='white')
    plt.plot(PDE_WL, PDE_WLEff, 'r', label='PDE $\lambda$')
    plt.plot(Wavelength, Intensity, 'b.', label='Scintillation $\lambda$')
    plt.plot(PDE_WL, CombinedEff, 'g', label='Combined Efficiency')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Normalised Intensity')
    plt.legend(loc='best')
#===============================================================================
#===============================================================================

for E in Energy:
    AcqTime = 0
    # Create new report file
    if (Verbose):
        Name = os.path.splitext(sys.argv[1])[0]
        Report = open(Name+'_Report_'+str(E)+'keV_'+str(NoOfEvents)+'evts_'+str(IdenticalRepeats)+'Reps_'+str(Period_us/1000)+'us_'+str(SampleStep)+'ns.txt', "w")


    if (Verbose):
        print colored('Using '+DarkCurrentCSV, 'blue')
        print colored('Using '+PDE_WavelengthCSV, 'blue')
        print colored('Bias : '+str(Bias)+'V', 'blue')
        print colored('Temperature : '+str(Temperature)+'C', 'blue')
        print colored('Input Energy : '+str(E)+'keV', 'blue')
        print '\n'
        Report.write( 'Using '+DarkCurrentCSV+'\n' )
        Report.write( 'Using '+PDE_WavelengthCSV+'\n')
        Report.write( 'Bias : '+str(Bias)+'V\n')
        Report.write( 'Temperature : '+str(Temperature)+'C\n')
        Report.write( 'Input Energy : '+str(E)+'keV\n')
        Report.write( '-----------------------------------------------------------------\n')

    DarkCurrentTable = ReadDarkCurrentTable(DarkCurrentCSV)

    if (Verbose):
        print 'Number of SiPMs : %d' %NoOfSiPMs
        Report.write('Number of SiPMs: %d \n' % NoOfSiPMs)

    SiPMArea = NoOfSiPMs * SiPMSize**2

    if (Verbose):
        print 'Total SiPM Area: %.1f mm2' % SiPMArea
        Report.write('Total SiPM Area: %.1f mm2 \n' % SiPMArea)
        Report.write('\n')
        print colored('------------------------ PHOTON TRANSPORT -----------------------', 'red')
        Report.write('------------------------ PHOTON TRANSPORT -----------------------\n')



    for Nevt in range(NoOfEvents):
        NPhotonsGenerated = EnergyToPhotons(E, PhotonsPerkeV)
        if (Verbose):
            print 'Photon Yield: ', NPhotonsGenerated
            Report.write('Photon Yield: %d\n' % NPhotonsGenerated)

        # Put in later
        #NPhotons = GetDoubleHitProbability(k, NoOfPixels, NPhotonsGenerated)
        # Delete this line...
        NPhotons = NPhotonsGenerated
        NPhotonsAtSiPM   = CrystalLosses(NPhotons, CrystalCrossSection, SiPMArea)

        if (Verbose):
            print 'NPhotons At SiPM Face: %d' % int(NPhotonsAtSiPM)
            Report.write('NPhotons At SiPM Face: %d\n' % int(NPhotonsAtSiPM))



        NPhotonsCollected, PDE_OverVolts, AverageEff = PhotonsCollected(PDE_WL, PDE_WLEff, PDE_OV, PDE_OVEff, CombinedEff, OverVoltage, NPhotonsAtSiPM)
        if (Verbose):
            print NPhotonsCollected
            print 'NPhotons Collected: ', int(sum(NPhotonsCollected))
            Report.write('NPhotons Collected %d \n' % int(sum(NPhotonsCollected)))
            print 'Average Efficiency: %.2f%%' % (100*AverageEff)
            Report.write('Average Efficiency: %.2f%% \n' % (100*AverageEff))
            print colored('--------------------- CROSSTALK & GAIN --------------------------', 'red')
            Report.write('-----------------------CROSSTALK & GAIN -------------------------\n')


        # ---------------------------------------------------------------------------------------------
        #                    CROSS TALK & GAIN CALCULATIONS
        # ---------------------------------------------------------------------------------------------
        CrossTalk    = GetCrossTalkPercentage(CrossTalkvsOV, OverVoltage)

        if (Verbose):
            print 'Cross Talk  ',
            Report.write('Cross Talk : \n')
            for i in range(NoOfSiPMs):
                print 'SiPM%d: %.1f%% '  % (i, 100*CrossTalk[i]),
                Report.write('SiPM%d: %.1f%%    '  % (i+1, 100*CrossTalk[i]))
            print ''
            Report.write('\n\n')

        Gain = GetGainDueToOvervoltage(GainOV, GainvsOV, OverVoltage)

        FiringPixelsIncCrossTalk, FiringPixels = GetFiringPixels(NPhotonsCollected, NoOfPixels, PDE_OverVolts, AverageEff, CrossTalk)
        if (Verbose):
            print 'Firing Pixels (Incl Crosstalk): %d \n' % FiringPixelsIncCrossTalk
            Report.write('Firing Pixels (Incl Crosstalk): %d \n' % FiringPixelsIncCrossTalk)
            print 'Firing Pixels (Excl Crosstalk): %d \n' % FiringPixels
            Report.write('Firing Pixels (Excl Crosstalk): %d \n' % FiringPixels)

        #NFire = RejectDoubleHits(NPhotonsAtSiPM, NoOfPixels)


        # ---------------------------------------------------------------------------------------------
        #                    DARK CURRENT CALCULATIONS
        # ---------------------------------------------------------------------------------------------
        if (Verbose):
            print colored('------------------ DARK CURRENT CALCULATONS --------------------', 'red')
            Report.write('------------------ DARK CURRENT CALCULATONS ---------------------\n')

        DarkCurrent  = GetDarkCurrent(NoOfSiPMs, DarkCurrentTable, Bias, Temperature)

        if (Verbose):
            print 'Dark Current (nA): ', DarkCurrent
            Report.write('Dark Current (nA): %.3f \n' % DarkCurrent)

        DarkRate = GetDarkRate(DarkCurrent, Gain)
        if (Verbose):
            print 'Dark Count Rate %.2f MHz: ' % (DarkRate/1e6)
            Report.write('Dark Count Rate %.2f MHz: \n' % (DarkRate/1e6))
            Report.write('\n\n')

            print colored('------------------------ OUTPUT CHARGE --------------------------', 'red')
            Report.write('------------------------ OUTPUT CHARGE --------------------------\n')


            print colored('--------------------------- RESOLUTION --------------------------', 'red')
            Report.write('--------------------------- RESOLUTION --------------------------\n')


        PhotonResolution = GetPhotonCountingResolution(GainVariance, FiringPixelsIncCrossTalk)


        if (Verbose):
            print 'Photon Resolution %.1f%%: ' %PhotonResolution
            Report.write('Photon Resolution %.1f%%: \n' %PhotonResolution)


        NoiseResolution  = GetNoiseResolution(DarkRate, IntegrationTime, FiringPixelsIncCrossTalk)


        if (Verbose):
            print 'Noise Resolution %.1f%%: ' %NoiseResolution
            Report.write('Noise Resolution %.1f%%: \n' %NoiseResolution)

        IntrinsicResolution = GetIntrinsicResolution(E, IntRes_vs_Energy, IntResolutionEnergy)
        if (Verbose):
            print 'Intrinsic Resolution %.1f%%: ' %IntrinsicResolution
            Report.write('Intrinsic Resolution %.1f%%: \n' %IntrinsicResolution)

        Resolution       = GetResolution(PhotonResolution, NoiseResolution, IntrinsicResolution)


        if (Verbose):
            print 'Resolution : %.1f %%' % Resolution
            Report.write('Resolution : %.1f %% \n' % Resolution)

        Charge_nC, Current_uA = PhotonsToCharge(FiringPixelsIncCrossTalk, Capacitance_nF, IntegrationTime, Gain)

        if (Verbose):
            print 'Charge : %.2f nC' %  Charge_nC
            print 'Current : %.2f uA' % Current_uA
            Report.write('Charge : %.2f nC  \n' %  Charge_nC)
            Report.write('Current : %.2f uA \n' % Current_uA)
            Report.write('\n\n')



        TauList = [tau_scint1, tau1_Intensity, tau_scint2, tau2_Intensity, SiPMtauRise, SiPMtauFall]

        ScintSignal, T = ChargeToPulse(Charge_nC, Capacitance_nF, Rq, Rs, TauList, SampleStart, SampleEnd, SampleStep)

        ScintSignal += GetNoise(NoiseOffset, NoiseLevel, len(T))
        if (PlotWaveform):
            PlotWaveforms(ScintSignal, T)


        AcqTime = WriteWaveformToTextFile(E, AcqTime, SampleStep, ScintSignal, Period_us)
        # If Identical repeated waveforms are being recorded,
        # skip non-interesting datapoints, not required by SPICE by putting a time gap in the output data


if (Verbose):
    Report.close()

plt.show()
