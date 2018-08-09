#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Started on 19-6-17

@author: Alan J Bell.
Extend the photon yield deviation beyond 1MeV based on extrapolation of existing data.

"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def main():
    E = range(10, 1100, 20)
    NPhotons = []
    NPhotonsLinear = []
    f_Factor = []
    PhotonsPerkeV = 54
    for i, e in enumerate(E):
        A, B, f = EnergyToPhotons(e, PhotonsPerkeV)
        NPhotons.append(A)
        NPhotonsLinear.append(B)
        f_Factor.append(f)
    PlotData(E, NPhotons, NPhotonsLinear)
    Plot_f(E, f_Factor)
    Extend_f(E, f_Factor)


def PlotData(E, NPhotons, NPhotonsLinear):
    fig1 = plt.figure(num=1, facecolor='white', figsize=(9,6))
    plt.plot(E, NPhotons, 'r-')
    plt.plot(E, NPhotonsLinear, 'k--')
    plt.semilogy()
    plt.draw()

def Plot_f(E, f_Factor):
    fig2 = plt.figure(num=2, facecolor='white', figsize=(9,6))
    plt.plot(E, f_Factor, 'r-', linewidth = 2, label='Fit to data')
    #plt.semilogx()
    plt.xlabel('Energy keV')
    plt.ylabel('f Factor')
    plt.draw()

def Expo(x, *p0):
    a, b, c = p0
    return a*np.exp(-b/x)


def Extend_f(E, f_Factor):
    sub_E = range(800, 1100, 10)
    C = 1.74317
    K = 0.003895
    m = -626.578
    sub_f = [1 + C*(1-(1/(1+np.exp(-K*(e-m))))) for e in sub_E]

    #line = np.polyfit(sub_E, sub_f, 1)
    #line = np.poly1d(line)
    p0 = [1.01, -12, 1]
    popt, pcov = curve_fit(Expo, sub_E, sub_f, p0=p0)
    print popt
    plt.plot(sub_E, Expo(sub_E, *popt),'b-')
    plt.semilogx()
    plt.draw()
    plt.pause(1)
    ext_E = range(1000, 3200, 100)
#    ext_f = line(ext_E)
#    plt.plot(ext_E, Expo(ext_E, *popt),'b.-')
#    plt.draw()
#    plt.pause(.01)
    plt.plot(ext_E, Expo(ext_E, *popt), 'b-.', linewidth=2, label='Extrapolation')
    plt.draw()
    LUT = open('f-Factor_LUT.csv','w')
    for i, e in enumerate(range(20, 880, 29)):
        f = 1 + C*(1-(1/(1+np.exp(-K*(e-m)))))
        print 'At  %.1f keV, f = %.4f ' % (float(e), float(f))
        LUT.write('%.1f , %.4f \n' % (float(e), float(f)))
    for i, e in enumerate(range(900, 3200, 50)):
        print 'At  %.1f keV, f = %.4f ' % (float(e), float(Expo(e, *popt)))
        LUT.write('%.1f , %.4f \n' % (float(e), float(Expo(e, *popt))))
    #plt.semilogx()

def EnergyToPhotons(E, PhotonsPerkeV):
    C = 1.74317
    K = 0.003895
    m = -626.578

    f = 1 + C*(1-(1/(1+np.exp(-K*(E-m)))))

    NPhotonsGenerated = int((E * PhotonsPerkeV * f)/6.0)
    NPhotonsLinear = int(E * PhotonsPerkeV)/6.0

    return NPhotonsGenerated, NPhotonsLinear, f


if __name__=='__main__':
    main()
    plt.show()
