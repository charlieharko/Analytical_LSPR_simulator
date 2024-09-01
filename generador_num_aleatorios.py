#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Size and AR random generator following a normal distribution. As input:
    
    Size = size of the nanorod
    S_stdev = Standard deviation of the size of the nanorod
    
    AR= aspect ration of the nanorod
    AR_stdev = Standard deviation of the AR 
    
    scale = number of nanoparticles to be generated (N)

@author: Carlos Renero Lecuna
"""

import numpy as np


def sizeDistribution (size, s_stdev, AR, AR_stdev, scale):
     
    size_dist = np.random.normal(size,s_stdev, scale)
    AR_dist = np.random.normal(AR,AR_stdev, scale)
        
    # Return sizes and AR normal distribution
    
    return [size_dist, AR_dist]