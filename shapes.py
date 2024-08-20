#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 23:10:54 2024

@author: lguser
"""

''' Input
-----
    L:  Length of particle, in m
    R:  Aspect ratio (length over width)

    Output
------
Row vector containing 
 - resonant permittivity eps1, 
 - mode volume V1,
 - size factors a12, a14,
 - particle volume V
'''

import numpy as np

# Definition of different shapes

# Index of shape:
             # 0: rod,   1: triangle,    2: ellipsoid,   3: disk,
             # 4: ring,  5: bipyramid,   6: cage,        7: bicone


def shapeRod(L, R):

    eps1 = -1.73 * R**(1.45) - 0.296
    a12 = 6.92 / (1-eps1)
    a14 = -6.69 / (1-eps1)
    V = np.pi * (3*R - 1) / (12*R**3) * L**3
    V1 = 0.896 * V;

    return [V, eps1, a12, a14, V1]

def shapeTriangle(L,R):
    eps1 = -0.87 * R**1.12 - 4.33;
    a12 = 5.57 / (1-eps1);
    a14 = -6.83 / (1-eps1);
    V = ( -0.00544 / R**2 + 0.433 / R )* L**3;
    V1 = V * (- 0.645 * R**(-1.24) + 0.678);
    
    return [V, eps1, a12, a14, V1]

def shapeEllipsoid(L,R):
    eps1 = -0.871 - 1.35 * R**1.54;
    a12 = 5.52 / (1-eps1);
    a14 = -9.75 / R**2.53;
    V = (np.pi /(6 * R**2) )* L**3;
    V1 = V * 0.994;
    
    return [V, eps1, a12, a14, V1]

def shapeDisk(L,R):
    eps1 = -0.479 - 1.36 * R**0.872;
    a12 = 7.05 / (1-eps1);
    a14 = -10.9 / R**0.98;
    V = ( 4 + 3 * (R-1) * (2*R + np.pi-2) ) * np.pi / (24*R**3) * L**3;
    V1 = V * 0.944;
        
    return [V, eps1, a12, a14, V1]

def shapeRing(L,R):
    eps1 = 1.39 - 1.31 * R**1.73; 
    a12 = 7.24 / (1-eps1);
    a14 = -19.1 / (1-eps1);
    V = np.pi**2 * (R - 1) / (4*R**3) * L**3;
    V1 = (0.514 + 2.07 / R**2.67) * V;
    
    return [V, eps1, a12, a14, V1]

def shapeBipyramid(L,R):
    eps1 = 1.43 - 4.52 * R*1.12; 
    a12 = 2.89 / (1-eps1);
    a14 = -1.79 / (1-eps1);
    V = 0.219 / R**2 * L**3;
    V1 = (1.96 - 1.73/R**0.207) * V;

    return [V, eps1, a12, a14, V1]

def shapeCage(L,R):
    eps1 = -0.0678 * R**2.02 - 3.42;
    a12 = -0.00405 * R**2.59 + 2.21;
    a14 = -13.9;
    V = (8.04 / R**3 - 12/ R**2 + 6 / R - 0.00138) * L**3;
    V1 = V * (-0.008 * R**2 + 0.103 * R + 0.316);

    return [V, eps1, a12, a14, V1]

def shapeBicone(L,R):
    eps1 = -0.687 - 2.54 * R**1.5;
    a12 = 1.34 / (1-eps1);
    a14 = -1.04 / (1-eps1);
    V = (0.262/R**2)* L**3;
    V1 = V * (0.648 - 0.441 / R**0.687);
    
    return [V, eps1, a12, a14, V1]