#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 23:54:02 2024

@author: lguser
"""

''' 
    Input
-----
    w: Row vector containing frecuency grid in 1/s; w = f *2pi

    Output
-----
    epsilonm: Row vector containing electrical permittivity of Au
'''

import numpy as np

# Fix constants
hbar = 6.5821195 * 10**(-16); # Dirac constant, in eVs
c = 2.99792458 * 10**8

# # Convert lambda (nm) in freq (Hz)
# c = 2.99792458 * 10**8
# x_nm = np.array(x)
# x_m = x_nm*10**(-9)
# w = 2*np.pi*c/x_m

def epsilonAu(w):
    # Constants for Au, from Yu, García de Abajo, Liz-Marzán
    #epsb = 9.5;
    wp = 9.06 / hbar #in 1/s
    Tau = 0.071 / hbar #in 1/s
    A = 0.132
    B = -1.755
    C = 20.43
    w1 = 2.43 / hbar #in 1/s
    w2 = 1.52 / hbar #in 1/s
    Tau1 = 0.0716 / hbar #in 1/s
    
    # Initiate vectors
    nPts = np.size(w);
    epsilonm = epsb = [0]*nPts
    w = np.array(w)
        
    # Create metal permittivity epsilonm vector
    epsb = A + B * np.log((w1-w-1j*Tau1)/(w1+w+1j*Tau1)) + C*np.exp(-w/w2)
    epsilonm = epsb - (wp**2 / (w * (w + 1j*Tau)))
    
    return [epsilonm, epsb]

def epsilonAg(w):
    # Constants for Au, from Yu, García de Abajo, Liz-Marzán
    # epsb = 4.0;
    wp = 9.17 /hbar #in 1/s
    Tau = 0.021 /hbar #in 1/s
    A = -9.71
    B = -1.111
    C = 13.77
    w1 = 4.02 /hbar #in 1/s
    w2 = 18.5 /hbar #in 1/s
    Tau1 = 0.076 /hbar #in 1/s
        
    # Initiate vectors
    nPts = np.size(w);
    epsilonm = epsb = [0]*nPts
    w = np.array(w)
        
    # Create metal permittivity epsilonm vector
    epsb = A + B * np.log((w1-w-1j*Tau1)/(w1+w+1j*Tau1)) + C*np.exp(-w/w2)
    epsilonm = epsb - (wp**2 / (w * (w + 1j*Tau)))
    
    return [epsilonm, epsb]

def epsilonCu(w):
    # Constants for Au, from Yu, García de Abajo, Liz-Marzán
    #epsb = 8.0;
    wp = 8.88 /hbar  #in 1/s
    Tau = 0.103 /hbar; #in 1/s
    A = -4.36
    B = -1.655
    C = 12.3
    w1 = 2.12 /hbar #in 1/s
    w2 = 5.43 /hbar #in 1/s
    Tau1 = 0.0528 /hbar #in 1/s

    # Initiate vectors
    nPts = np.size(w);
    epsilonm = epsb = [0]*nPts
    w = np.array(w)
        
    # Create metal permittivity epsilonm vector
    epsb = A + B * np.log((w1-w-1j*Tau1)/(w1+w+1j*Tau1)) + C*np.exp(-w/w2)
    epsilonm = epsb - (wp**2 / (w * (w + 1j*Tau)))
    
    return [epsilonm, epsb]