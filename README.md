Created on Sun Aug 11 21:23:27 2024
@author: Carlos RENERO LECUNA based on the MatLab version of Natalie Fehn
    
Easy and fast calculation of extinction, scattering, absorbance cross-sections of plasmonic nanoparticles together with their
Quantum Yield and also the dielectric function of Au, Ag and Cu with given shape, dimensions, and material.

The fast calculation is due to the analytical nature of the calculations. This work is part of the project of Natalie Fehn based on the

Tutorial review by Renwen Yu, Luis M. Liz-Marzán, and F. Javier García de Abajo (Chemical Society Reviews, 2017).

The original Matlab code was made by Natalie Fehn and the translation into Python was done by Carlos Renero Lecuna

UPDATED - V.1.2   2024.08.22
+++++++++++++++++++++++++++++
The new version allows to simulate the LSPR with the polarization of the electromagnetic field along 
the transversal axis of the nanorods (only) to excite the transversal modes of the LSPR.

UPDATED - V.2     2024.09.01
+++++++++++++++++++++++++++++
This new version takes into account the size distribution of the NPs to generate a more accurate LSPR simulation,
where the SD of the size and AR can be introduced and a normal size distribution is generated

UPDATED - V.2.1    2024.10.23
+++++++++++++++++++++++++++++
This new release allows the user to use of different seeds in the generation of the size distribution of length and AR of the nanoparticles, so opens up
the possibility of comparing results of reproducibility. Also, includes in the lower part the size distribution histograms and violin plots used for the simulations.

UPDATED - V.2.1.1   2024.10.30
+++++++++++++++++++++++++++++
This version includes a new normalization that takes into account the 1/3 and 2/3 factors for Transversal and Longitudinal LSPR respectively
