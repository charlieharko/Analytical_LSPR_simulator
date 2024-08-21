#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:42:57 2024

@author: Carlos RENERO LECUNA
"""

import streamlit as st
import numpy as np
import pandas as pd
import giveResonance as gr
import plotly.express as px


# The plots will fill all the webpage space
st.set_page_config(layout="wide")

# Title of the webpage
st.write("""
         # LSPR simulation under longitudinal E excitation
         based on: Renwen Yu, Luis M. Liz-Marzán, and F. Javier García de Abajo (Chemical Society Reviews, 2017)
Different materials (Au, Ag, Cu) and shapes can be selected

The simulation ONLY takes into account the electric field in the longitudinal dimension of the NP
       """)
         
# Simulation parameter: x, R, l, shape, mat, refind
st.sidebar.write("""
         # Choose the Simulation Parameters
         """)
spectral_range = st.sidebar.slider(
    "Spectral range to be simulated (nm):", 190, 2500, value=(350,1200), 
    key="run_every", help = 'Select the maximun and minimun values of the ' 
    'wavelengths to to run the simulation. The step is always 1 nm.'
)

length = st.sidebar.text_input("Introduce the length of the nanoparticle (nm):",
                               value='20', key="size")

AR = st.sidebar.text_input("Introduce the aspect ratio",
                               value='1', key="AR", help = "The AR is the division"
                               " of the length over the size of the NP.")

shapes = ['Nanorod', 'Triangle', 'Ellipsoid', 'Disk', 'Ring','Bipyramid',
            'Cage','Bicone']

shp = st.sidebar.selectbox("Select the shape of the nanoparticle:", 
                           shapes, 
                            help = 'Select the shape of the nanoparticle '
                            'to be simulated.', index = 0)

n = st.sidebar.text_input("Introduce refractive index of the media:",
                               value='1', key="refind", help = "The n of water "
                               "is 1.33 and the vacuum is 1.")

materials = ['Gold', 'Silver', 'Copper']
mtr = st.sidebar.selectbox("Select the material of the nanoparticle:", 
                           materials, 
                            help = 'Select the material of the nanoparticle '
                            'to be simulated.', index = 0)

# Data to feed the giveResonance method

x = list(range(spectral_range[0], spectral_range[1], 1)) # Spectral range in nm
R = float(AR)                   # Aspect Ration (taken from the text)
l = float(length)               # Length of the nanoparticle (nm)
shape = shapes.index(shp)       # Shapes
mat = materials.index(mtr)      # Materials
refind = float(n)               # Medium refractive index

# Running the method
Results = gr.giveResonance(x, R, l, shape, mat, refind)

# Chart organization in columns
col1, col2 = st.columns(2)

# Plotting the LSPR
col1.write("""
        ## LSPR simulation
        """)
        
# Preparing the data for the LSPR
column_name = ['Ext', 'Sca', 'Abs']
A = pd.DataFrame(list(zip(Results[0],Results[1], Results[2])), x, 
                 columns = column_name)

# Creating the Figures
fig1 = px.line(A)
# Set x-axis title
fig1.update_xaxes(title_text="Wavelength (nm)")
# Set y-axes titles
fig1.update_yaxes(title_text="Extinction Cross section")
# Plot the figure
col1.plotly_chart(fig1, use_container_width=True)

# Plotting the QY (%)
col2.write("""
        ## Quantum Yield
        """)
# Preparing the data for the LSPR
QY = ['QY']
B = pd.DataFrame(Results[3]*100, x, columns = QY)
# Creating the Figures
fig2 = px.line(B)
# Set x-axis title
fig2.update_xaxes(title_text="Wavelength (nm)")
# Set y-axes titles
fig2.update_yaxes(title_text="Quantum Yield (%)")
# Plot the figure
col2.plotly_chart(fig2, use_container_width=True)


# We plot the dielectric functions Eps and Epb imag and real
col1.write("""
          ## Imaginary part of Dielectric functions
          """)

# Data for the dielectric functions of the material and the media
# We separte real and imaginary parts
eps_imag = np.imag(Results[4])
eps_real = np.real(Results[4])

epb_imag = np.imag(Results[5])
epb_real = np.real(Results[5])


column_name2 = ['Eps_imag', 'Epb_imag']
C = pd.DataFrame(list(zip(eps_imag, epb_imag)), x, 
                  columns = column_name2)         
# Creating the Figures
fig3 = px.line(C)
# Set x-axis title
fig3.update_xaxes(title_text="Wavelength (nm)")
# Set y-axes titles
fig3.update_yaxes(title_text="Imaginary part Eps / Epb")
# Plot the figure
col1.plotly_chart(fig3, use_container_width=True)


col2.write("""
          ## Real part of Dielectric functions
          """)
         
column_name3 = ['Eps_real', 'Epb_real']
D = pd.DataFrame(list(zip(eps_real, epb_real)), x, 
                  columns = column_name3)
# Creating the Figures
fig4 = px.line(D)
# Set x-axis title
fig4.update_xaxes(title_text="Wavelength (nm)")
# Set y-axes titles
fig4.update_yaxes(title_text="Real part Eps / Epb")
# Plot the figure
col2.plotly_chart(fig4, use_container_width=True)

# Download the generated data
dfs = [B, C, D]
data = A.join(dfs)
data = data.rename_axis("Wavelength (nm)")
csv = data.to_csv().encode("utf-8")

st.sidebar.download_button("Download data as CSV.csv",
                            data = csv,
                            file_name = "Simulated data.txt",
                            mime="txt/csv")
