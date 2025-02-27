#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 22:24:38 2024

@author: lguser
"""

import generador_num_aleatorios as gdist
import giveResonance as gr
import numpy as np

def size_modified_LSPR (x, R, AR_stdev, l, l_stdev,
                        shape, mat, refind, mode, nPts,seed):
    
    
    rows = np.size(x) # Rows

    # Normal distribution for values of the size of the NPs
    sizes = gdist.sizeDistribution(l, l_stdev, R, AR_stdev, nPts,seed)

    # Simulate the LSPR of all the NPs 

    LSPR = 0
    # Array containing the results for Ext, Sca and Abs
    Ext = np.zeros([rows, nPts])
    Sca = np.zeros([rows, nPts])
    Abs = np.zeros([rows, nPts])
    QY = np.zeros([rows, nPts]) 

    for i in range(nPts):
        LSPR = gr.giveResonance(x, sizes[1][i], sizes[0][i], shape, mode, 
                                    mat, refind)
        
        Ext[:,i] = LSPR[0]
        Sca[:,i] = LSPR[1]
        Abs[:,i] = LSPR[2]
        QY[:,i] = LSPR[3]   

    # Removing possible NaN from the data
    Ext = Ext[:, ~np.isnan(Ext).any(axis=0)]
    Sca = Sca[:, ~np.isnan(Sca).any(axis=0)]
    Abs = Abs[:, ~np.isnan(Abs).any(axis=0)]
    QY = QY[:, ~np.isnan(QY).any(axis=0)]

    # Axis = 1 across columns & Axis = 0 across rows
    Ext_av = np.mean(Ext,axis=1)   
    Sca_av = np.mean(Sca,axis=1)   
    Abs_av = np.mean(Abs,axis=1)
    Y_av = Abs_av/Ext_av
    # Y_av = np.mean(QY,axis=1)
    
    #The rest of results, that do not averate
    
    epsilonm = LSPR[4]
    epsb = LSPR[5]
    lamdamax = LSPR[6]
    Ymax = LSPR[7]
    diam = sizes[0]
    AR = sizes[1]
    
    return [Ext_av, Sca_av, Abs_av, Y_av, epsilonm, epsb, lamdamax, Ymax, diam, AR]
