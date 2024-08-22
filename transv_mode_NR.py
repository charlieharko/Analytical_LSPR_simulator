#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 13:34:56 2024

@author: Carlos Renero-Lecuna

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

"""



import numpy as np

# Definition of different shapes

# Index of shape:
             # 0: mode 1,   1: mode 2,    2: mode 3


def mode1(L, R):

    eps1 = -1.75 + 3.19 / R**(6.14)
    a12 = 0.0148 + 3.69  / R**(2.85)
    a14 = 0.0142 - 16.9 / R**(3.58)
    V = np.pi * (3*R - 1) / (12*R**3) * L**3
    V1 = (0.0679 + 1.83 / R**(2.1)) * V

    return [V, eps1, a12, a14, V1]

def mode2(L, R):

    eps1 = -0.978 - 0.661 / R**(1.1)
    a12 = -21.7 + 22.7  / R**(0.0232)
    a14 = 1.48 - 3.67 / R**(0.458)
    V = np.pi * (3*R - 1) / (12*R**3) * L**3
    V1 = (0.891 - 2.28 / R**(2.53)) * V

    return [V, eps1, a12, a14, V1]

def mode3(L, R):

    eps1 = -1.57 + 0.0446 * R
    a12 = -0.0117 + 0.773  / R**(1.46)
    a14 = -0.256 + 0.0554 * R**(0.758)
    V = np.pi * (3*R - 1) / (12*R**3) * L**3
    V1 = (-0.0346 + 0.0222 * R) * V

    return [V, eps1, a12, a14, V1]