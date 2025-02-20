#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Aug 11 21:23:27 2024
@author: Carlos RENERO LECUNA based on the MatLab version of Natalie Fehn
    
Easy and fast calculation of extinction cross-section, scattering cross-section, and quantum yield
of nanoparticles with given shape, dimensions, and material.

Tutorial review: Renwen Yu, Luis M. Liz-Marzán, and F. Javier García de Abajo (Chemical Society Reviews, 2017).
Matlab code: Natalie Fehn
Pythond code: Carlos Renero Lecuna
"""

# Librarie to import

import numpy as np
import shapes as sh
import transv_mode_NR as tr
import eps_function as eps

''' List of parameters to inittiate the simulation if not present in the mail
# x = list(range(300, 1201, 1)) # Spectral range in nm

# R = float(1) # Aspect Ratio

# l = float(20) # Length of the nanoparticle (nm)

# shape = 1    # Index of shape:
#              # 1: rod,   2: triangle,    3: ellipsoid,   4: disk,
#              # 5: ring,  6: bipyramid,   7: cage,        8: bicone

# mat = 1      # Material:
#              # 1: gold,  2: silver,      3: copper
'''

def giveResonance(x, R, l, shape, mode, mat, refind):
    
    # Defining constants
    epsh = refind**2            # relative permittivity of host medium, water epsh = n^2
    #refind = np.sqrt(epsh)     # refractive index of host medium
    c = 2.99792458 * 10**8      # speed of light in vacuum, in m/s

    # Create wavelength and frecuency vectors
    nPts = np.size(x)           # number of points of wavelength grid
    x = np.array(x)             # Transform it to manipulate it mathematically
    lambda_m = x*10**(-9)       # row vector wavelength, in m
    w = 2*np.pi*c/ lambda_m     # row vector frecuency, in 1/s
    
    # Create shape dependent constants
    L = l * 10**(-9)            # characteristic length in m

    if (shape == 0):
        [V, eps1, a12, a14, V1] = sh.shapeRod(L, R)
    elif (shape == 1):
        [V, eps1, a12, a14, V1] = sh.shapeTriangle(L, R)
    elif (shape == 2):
        [V, eps1, a12, a14, V1] = sh.shapeEllipsoid(L, R)
    elif (shape == 3):
        [V, eps1, a12, a14, V1] = sh.shapeDisk(L, R)
    elif (shape == 4):
        [V, eps1, a12, a14, V1] = sh.shapeRing(L, R)
    elif (shape == 5):
        [V, eps1, a12, a14, V1] = sh.shapeBipyramid(L, R)
    elif (shape == 6):
        [V, eps1, a12, a14, V1] = sh.shapeCage(L, R)
    elif (shape == 7):
        [V, eps1, a12, a14, V1] = sh.shapeBicone(L, R)
    elif (shape == 8): 
        print('No more shapes')
        return [V, eps1, a12, a14, V1]
            
    if (mode == 1):
        [V, eps1, a12, a14, V1] = tr.mode1(L, R)
    elif (mode == 2):
        [V, eps1, a12, a14, V1] = tr.mode2(L, R)
    elif (mode == 3):
        [V, eps1, a12, a14, V1] = tr.mode3(L, R)
    elif (shape == 8): 
        print('No more shapes')
        return [V, eps1, a12, a14, V1]
        
    # Factor for calculation of perturbation, 3rd order of s
    a13 = 4 * np.pi**2 * 1j * V1 / (3*L**3)
    
    # Create metal permittivity
    if (mat == 0):
      [epsilonm,epsb] = eps.epsilonAu(w)
    elif (mat == 1):
      [epsilonm,epsb] = eps.epsilonAg(w)
    elif (mat == 2):
      [epsilonm,epsb] = eps.epsilonCu(w)
    else:
        print("There is an error")
        
        return [epsilonm, epsb]
  
    # Initiate vectors
    sigma_ext = np.array([0]*nPts)  # extinction cross-section, in m^2
    a = np.array([0]*nPts)          # polarisability, in m^3
    s = np.array([0]*nPts)          # size factor
    A1 = np.array([0]*nPts)         # retardation effects
    sigma_sca = np.array([0]*nPts)  # scattering cross-section, in m^2
    sigma_abs = np.array([0]*nPts)  # absorption cross-section, in m^2
    Y = np.array([0]*nPts)          # quantum yield  
  
    # Create vectors for size, retardation and polarisability
    s = refind * L / lambda_m
    A1 = a12 * s**2 + a13 * s**3 + a14 * s**4
    a = epsh / (np.pi * 4) * V1/(1/(epsilonm/epsh - 1) - 1 /(eps1 - 1) - A1)
    
    # Create vector for extinction cross-section
    sigma_ext = 8 * np.pi**2 / (refind * lambda_m) * a.imag
    
    # Create vectors for scattering cross-section 
    sigma_sca = 128 * np.pi**5 / (3 * lambda_m**4) * np.abs(a)**2
    
    #Create a vector for absorption cross-section
    sigma_abs = sigma_ext - sigma_sca
    
    # Create a vector for quantum yield (defined by Abs/Ext)
    Y = sigma_abs / sigma_ext # Quanutm Yield
    
    # Display resonance wavelength and quantum yield at resonance wavelength
    index_max = np.argmax(sigma_ext)
    
    # Maximum of extinction cross-section and its index
    lamdamax = x[index_max]
    Ymax = Y[index_max]
    
    # Return all the resultst
    return [sigma_ext, sigma_sca, sigma_abs, 
            Y, epsilonm, epsb, lamdamax, Ymax]
