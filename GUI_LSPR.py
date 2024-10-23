#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sunday Sept 1 22:39:57 2024

Modification of the program GUI_LSPR to see the influence of size distribution
of the nanoparticles produced, applying a normal size distribution of the sample

@author: Carlos RENERO LECUNA
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import size_mod_LSPR as smLSPR


# The plots will fill all the webpage space
st.set_page_config(layout="wide")

# Title of the webpage
st.write("""
         # LSPR simulation under longitudinal and/or transversal excitation
         based on: Renwen Yu, Luis M. Liz-Marzán, and F. Javier García de Abajo (Chemical Society Reviews, 2017)

Different materials (Au, Ag, Cu), size distributions and shapes can be simulated
       """)
         
# Simulation parameter: x, R, l, shape, mat, refind
st.sidebar.write("""
         # Simulation Parameters
         """)
spectral_range = st.sidebar.slider(
    "Spectral range to be simulated (nm):", 190, 2500, value=(350,1200), 
    key="run_every", help = 'Select the maximun and minimun values of the ' 
    'wavelengths to to run the simulation. The step is always 1 nm.'
)

semilla = st.sidebar.text_input("Random Generator Seed (when using size distribution):",
                               value='1', key="seed", help = "When using size distribution, you can introduce"
                               " a seed number in the number random generator for reproducibility tests. If you "
                               "leave it blank, the simulation will change when changing parameters.")

materials = ['Gold', 'Silver', 'Copper']
mtr = st.sidebar.selectbox("Select the material of the nanoparticle:", 
                           materials, 
                            help = 'Select the material of the nanoparticle '
                            'to be simulated.', index = 0)

shapes = ['Nanorod', 'Triangle', 'Ellipsoid', 'Disk', 'Ring','Bipyramid',
                        'Cage','Bicone']

shp = st.sidebar.selectbox("Select the shape of the nanoparticle:", 
                           shapes, 
                            help = 'Select the shape of the nanoparticle '
                            'to be simulated.', index = 0)

# Disabled option to calculate transversal modes when other option rather
# than Nanorod is selected
def transmode (shp):
    if (shp == 'Nanorod'):
        dis = False
    else:
        dis = True
    
    return dis

dis = transmode(shp)

# Checkbox to be activated only for Nanorod shape in the case that the
# transversal modes want to be calculated

agree = st.sidebar.checkbox("If you select Nanorod, would you like to calculate the transversal modes?",
                            disabled = dis)

modes = ['Mode 1', 'Mode 2', 'Mode 3', 'Average']

# Select the transversal mode that you want to calculate
modes_option = st.sidebar.selectbox("Select the transversal mode:", 
                           modes, 
                            help = 'Select the transveral mode of the nanoparticle '
                            'to be simulated.', index = 0, disabled = not agree)


length = st.sidebar.text_input("Introduce the length of the nanoparticle (nm):",
                               value='20', key="size")

SD_length = st.sidebar.text_input("Introduce the SD of the length of the nanopartile (nm):",
                                     value = '0', key='SD_size', help='If you do not want to consider it, write 0')

AR = st.sidebar.text_input("Introduce the aspect ratio",
                               value='1', key="AR", help = "The AR is the division"
                               " of the length over the size of the NP.")

AR_stedv = st.sidebar.text_input("Introduce the SD of the AR of the nanopartile:",
                                     value = '0', key='SD_AR', help='If you do not want to consider it, write 0')

nPts = st.sidebar.text_input("Introduce the number of nanopartiles to make the size distribution:",
                                     value = '1', key='nPts', help='If you do not want to consider it, write 1')

# Introduce the refractive index of the media 
n = st.sidebar.text_input("Introduce refractive index of the media:",
                               value='1.33', key="refind", help = "The n of water "
                               "is 1.33 and the vacuum is 1.")

# Data to feed the giveResonance method

x = list(range(spectral_range[0], spectral_range[1], 1)) # Spectral range in nm

# When entering blank that cannot be converted into integer - no seed
try:
   val = int(semilla)
except ValueError:
   # foo is something that cannot be converted to
   # a number. It could be an empty string, or a
   # string like 'hello'
   # provide a default value
   val = None
seed = val

R = float(AR)                   # Aspect Ration (taken from the text)
AR_stdev = float(AR_stedv)       # SD for AR            
l = float(length)               # Length of the nanoparticle (nm)
l_stdev = float(SD_length)      # SD for NPs length
nPts = int(nPts)              # Numer of nanoparticles in the size distribution
shape = shapes.index(shp)       # Shapes
mode = agree                    # Transversal mode calculation ONLY for rods
mat = materials.index(mtr)      # Materials
refind = float(n)               # Medium refractive index

# Running the method
# Results = gr.giveResonance(x, R, l, shape, mode, mat, refind)

Results = smLSPR.size_modified_LSPR (x, R, AR_stdev, l, l_stdev,
                        shape, mat, refind, mode, nPts,seed)

# All the options
if (mode == 0):
    # Chart organization in columns
    col1, col2 = st.columns(2)
    
    # Plotting the LSPR
            
    # Preparing the data for the LSPR
    column_name = ['Ext', 'Sca', 'Abs']
    A = pd.DataFrame(list(zip(Results[0],Results[1], Results[2])), x, 
                     columns = column_name)
    
    # Creating the Figures
    fig1 = px.line(A, title='Longitudinal LSPR analytical simulation')
    # Set x-axis title
    fig1.update_xaxes(title_text="Wavelength (nm)")
    # Set y-axes titles
    fig1.update_yaxes(title_text="Extinction Cross section")
    # Plot the figure
    col1.plotly_chart(fig1, use_container_width=True)
    
    # Plotting the QY (%)
    
    # Preparing the data for the LSPR
    QY = ['QY']
    B = pd.DataFrame(Results[3]*100, x, columns = QY)
    # Creating the Figures
    fig2 = px.line(B, title='Quantum Yield (%) - Scattering / Extinction cross section')
    # Set x-axis title
    fig2.update_xaxes(title_text="Wavelength (nm)")
    # Set y-axes titles
    fig2.update_yaxes(title_text="Quantum Yield (%)")
    # Plot the figure
    col2.plotly_chart(fig2, use_container_width=True)
    
    
    # We plot the dielectric functions Eps and Epb imag and real

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
    fig3 = px.line(C, title='Imaginary part of dielectric function')
    # Set x-axis title
    fig3.update_xaxes(title_text="Wavelength (nm)")
    # Set y-axes titles
    fig3.update_yaxes(title_text="Imaginary part Eps / Epb")
    # Plot the figure
    col1.plotly_chart(fig3, use_container_width=True)
    
             
    column_name3 = ['Eps_real', 'Epb_real']
    D = pd.DataFrame(list(zip(eps_real, epb_real)), x, 
                      columns = column_name3)
    # Creating the Figures
    fig4 = px.line(D, title='Real part of dielectric function')
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

else:
    
    # Chart organization in columns
    col1, col2 = st.columns(2)

    # Plotting the LSPR
            
    if (modes_option == 'Mode 1'):
        
        Results0 = smLSPR.size_modified_LSPR(x, R, AR_stdev, l, l_stdev, shape, mat, refind, 0, nPts, seed) # Long
        Results1 = smLSPR.size_modified_LSPR(x, R, AR_stdev, l, l_stdev, shape, mat, refind, 1, nPts, seed) # Trans
        
        # Only longitudinal mode
        Ext_mode0 = Results0[0]
        Sca_mode0 = Results0[1]
        Abs_mode0 = Results0[2]
        # Only transversal mode
        Ext_mode1 = Results1[0]
        Sca_mode1 = Results1[1]
        Abs_mode1 = Results1[2]
        
        # Both modes combined
        Ext_average = (Results0[0]+Results1[0])/2
        Sca_average = (Results0[1]+Results1[1])/2
        Abs_average = (Results0[2]+Results1[2])/2
        
        QY_average = (Results0[3]+Results1[3])/2
        
        # Preparing the data for the LSPR
        column_name = ['Ext', 'Sca', 'Abs']
        A = pd.DataFrame(list(zip(Ext_average,Sca_average, Abs_average)), x, 
                         columns = column_name)

        # Transversal
        column_name_trans = ['Ext_trans', 'Sca_trans', 'Abs_trans']
        A_trans = pd.DataFrame(list(zip(Ext_mode1,Sca_mode1, Abs_mode1)), x, 
                         columns = column_name_trans)
        # Longitudinal
        column_name_long = ['Ext_long', 'Sca_long', 'Abs_long']
        A_long = pd.DataFrame(list(zip(Ext_mode0,Sca_mode0, Abs_mode0)), x, 
                 columns = column_name_long)

        # Creating the Figures
        fig1 = px.line(A, title='Complete LSPR analytical simulation')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col1.plotly_chart(fig1, use_container_width=True)
        
        # Plotting the QY (%)
        QY = ['QY']
        B = pd.DataFrame(QY_average*100, x, columns = QY)
        # Creating the Figures
        fig2 = px.line(B, title='Quantum Yield (%) - Scattering / Extinction cross section')
        # Set x-axis title
        fig2.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig2.update_yaxes(title_text="Quantum Yield (%)")
        # Plot the figure
        col2.plotly_chart(fig2, use_container_width=True)
        
        # Creating the Figures (Transversal part)
        fig1 = px.line(A_trans, title='Transversal (mode 1) contribution to LSPR')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col1.plotly_chart(fig1, use_container_width=True)
        
        # Creating the Figures (Long. mode)
        fig1 = px.line(A_long, title='Longitudinal contribution to LSPR')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col2.plotly_chart(fig1, use_container_width=True)
        
        # We plot the dielectric functions Eps and Epb imag and real

        # Data for the dielectric functions of the material and the media
        # We separte real and imaginary parts
        eps_imag = np.imag(Results0[4])
        eps_real = np.real(Results0[4])

        epb_imag = np.imag(Results0[5])
        epb_real = np.real(Results0[5])

        column_name2 = ['Eps_imag', 'Epb_imag']
        C = pd.DataFrame(list(zip(eps_imag, epb_imag)), x, 
                          columns = column_name2)         
        # Creating the Figures
        fig3 = px.line(C, title='Imaginary part of dielectric function')
        # Set x-axis title
        fig3.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig3.update_yaxes(title_text="Imaginary part Eps / Epb")
        # Plot the figure
        col1.plotly_chart(fig3, use_container_width=True)

        column_name3 = ['Eps_real', 'Epb_real']
        D = pd.DataFrame(list(zip(eps_real, epb_real)), x, 
                          columns = column_name3)
        # Creating the Figures
        fig4 = px.line(D, title='Real part of dielectric function')
        # Set x-axis title
        fig4.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig4.update_yaxes(title_text="Real part Eps / Epb")
        # Plot the figure
        col2.plotly_chart(fig4, use_container_width=True)
        
        # Download the generated data
        dfs = [A_trans, A_long, B, C, D]
        data = A.join(dfs)
        data = data.rename_axis("Wavelength (nm)")
        csv = data.to_csv().encode("utf-8")
    
        st.sidebar.download_button("Download data as CSV.csv",
                                    data = csv,
                                    file_name = "Simulated data.txt",
                                    mime="txt/csv")

    elif (modes_option == 'Mode 2'):
        
        Results0 = smLSPR.size_modified_LSPR(x, R, AR_stdev, l, l_stdev, shape, mat, refind, 0, nPts, seed)
        Results2 = smLSPR.size_modified_LSPR(x, R, AR_stdev, l, l_stdev, shape, mat, refind, 2, nPts, seed)
        
        # Only longitudinal mode
        Ext_mode0 = Results0[0]
        Sca_mode0 = Results0[1]
        Abs_mode0 = Results0[2]
        # Only transversal mode
        Ext_mode2 = Results2[0]
        Sca_mode2 = Results2[1]
        Abs_mode2 = Results2[2]
        
        # Both modes combined
        Ext_average = (Results0[0]+Results2[0])/2
        Sca_average = (Results0[1]+Results2[1])/2
        Abs_average = (Results0[2]+Results2[2])/2
        
        QY_average = (Results0[3]+Results2[3])/2
        
        # Preparing the data for the LSPR
        column_name = ['Ext', 'Sca', 'Abs']
        A = pd.DataFrame(list(zip(Ext_average,Sca_average, Abs_average)), x, 
                         columns = column_name)

        # Transversal
        column_name_trans = ['Ext_trans', 'Sca_trans', 'Abs_trans']
        A_trans = pd.DataFrame(list(zip(Ext_mode2,Sca_mode2, Abs_mode2)), x, 
                         columns = column_name_trans)
        # Longitudinal
        column_name_long = ['Ext_long', 'Sca_long', 'Abs_long']        
        A_long = pd.DataFrame(list(zip(Ext_mode0,Sca_mode0, Abs_mode0)), x, 
                 columns = column_name_long)
        

        # Creating the Figures
        fig1 = px.line(A, title='Complete LSPR analytical simulation')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col1.plotly_chart(fig1, use_container_width=True)
        
        # Plotting the QY (%)

        # Preparing the data for the LSPR
        QY = ['QY']
        B = pd.DataFrame(QY_average*100, x, columns = QY)
        # Creating the Figures
        fig2 = px.line(B, title='Quantum Yield (%) - Scattering / Extinction cross section')
        # Set x-axis title
        fig2.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig2.update_yaxes(title_text="Quantum Yield (%)")
        # Plot the figure
        col2.plotly_chart(fig2, use_container_width=True)
        
        # Creating the Figures (Transversal part)
        fig1 = px.line(A_trans, title='Transversal (mode 2) contribution to LSPR')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col1.plotly_chart(fig1, use_container_width=True)
        
        # Creating the Figures (Long. mode)
        fig1 = px.line(A_long, title='Longitudinal contribution to LSPR')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col2.plotly_chart(fig1, use_container_width=True)
        
        # We plot the dielectric functions Eps and Epb imag and real

        # Data for the dielectric functions of the material and the media
        # We separte real and imaginary parts
        eps_imag = np.imag(Results0[4])
        eps_real = np.real(Results0[4])

        epb_imag = np.imag(Results0[5])
        epb_real = np.real(Results0[5])

        column_name2 = ['Eps_imag', 'Epb_imag']
        C = pd.DataFrame(list(zip(eps_imag, epb_imag)), x, 
                          columns = column_name2)         
        # Creating the Figures
        fig3 = px.line(C, title='Imaginary part of dielectric function')
        # Set x-axis title
        fig3.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig3.update_yaxes(title_text="Imaginary part Eps / Epb")
        # Plot the figure
        col1.plotly_chart(fig3, use_container_width=True)

        column_name3 = ['Eps_real', 'Epb_real']
        D = pd.DataFrame(list(zip(eps_real, epb_real)), x, 
                          columns = column_name3)
        # Creating the Figures
        fig4 = px.line(D, title='Real part of dielectric function')
        # Set x-axis title
        fig4.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig4.update_yaxes(title_text="Real part Eps / Epb")
        # Plot the figure
        col2.plotly_chart(fig4, use_container_width=True)
        
        # Download the generated data
        dfs = [A_trans, A_long, B, C, D]
        data = A.join(dfs)
        data = data.rename_axis("Wavelength (nm)")
        csv = data.to_csv().encode("utf-8")
    
        st.sidebar.download_button("Download data as CSV.csv",
                                    data = csv,
                                    file_name = "Simulated data.txt",
                                    mime="txt/csv")
    
    elif (modes_option == 'Mode 3'):     
        
        Results0 = smLSPR.size_modified_LSPR(x, R, AR_stdev, l, l_stdev, shape, mat, refind, 0, nPts, seed)
        Results3 = smLSPR.size_modified_LSPR(x, R, AR_stdev, l, l_stdev, shape, mat, refind, 3, nPts, seed)
        
        # Only longitudinal mode
        Ext_mode0 = Results0[0]
        Sca_mode0 = Results0[1]
        Abs_mode0 = Results0[2]
        # Only transversal mode
        Ext_mode3 = Results3[0]
        Sca_mode3 = Results3[1]
        Abs_mode3 = Results3[2]
        
        # Both modes combined
        Ext_average = (Results0[0]+Results3[0])/2
        Sca_average = (Results0[1]+Results3[1])/2
        Abs_average = (Results0[2]+Results3[2])/2
        
        QY_average = (Results0[3]+Results3[3])/2
        
        # Preparing the data for the LSPR
        column_name = ['Ext', 'Sca', 'Abs']
        A = pd.DataFrame(list(zip(Ext_average,Sca_average, Abs_average)), x, 
                         columns = column_name)

        # Transversal
        column_name_trans = ['Ext_trans', 'Sca_trans', 'Abs_trans']
        A_trans = pd.DataFrame(list(zip(Ext_mode3,Sca_mode3, Abs_mode3)), x, 
                         columns = column_name_trans)
        # Longitudinal
        column_name_long = ['Ext_long', 'Sca_long', 'Abs_long']        
        A_long = pd.DataFrame(list(zip(Ext_mode0,Sca_mode0, Abs_mode0)), x, 
                 columns = column_name_long)
        

        # Creating the Figures
        fig1 = px.line(A, title='Complete LSPR analytical simulation')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col1.plotly_chart(fig1, use_container_width=True)
        
        # Plotting the QY (%)

        # Preparing the data for the LSPR
        QY = ['QY']
        B = pd.DataFrame(QY_average*100, x, columns = QY)
        # Creating the Figures
        fig2 = px.line(B, title='Quantum Yield (%) - Scattering / Extinction cross section')
        # Set x-axis title
        fig2.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig2.update_yaxes(title_text="Quantum Yield (%)")
        # Plot the figure
        col2.plotly_chart(fig2, use_container_width=True)
        
        # Creating the Figures (Transversal part)
        fig1 = px.line(A_trans, title='Transversal (mode 3) contribution to LSPR')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col1.plotly_chart(fig1, use_container_width=True)
        
        # Creating the Figures (Long. mode)
        fig1 = px.line(A_long, title='Longitudinal contribution to LSPR')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col2.plotly_chart(fig1, use_container_width=True)
        
        # We plot the dielectric functions Eps and Epb imag and real

        # Data for the dielectric functions of the material and the media
        # We separte real and imaginary parts
        eps_imag = np.imag(Results0[4])
        eps_real = np.real(Results0[4])

        epb_imag = np.imag(Results0[5])
        epb_real = np.real(Results0[5])

        column_name2 = ['Eps_imag', 'Epb_imag']
        C = pd.DataFrame(list(zip(eps_imag, epb_imag)), x, 
                          columns = column_name2)         
        # Creating the Figures
        fig3 = px.line(C, title='Imaginary part of dielectric function')
        # Set x-axis title
        fig3.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig3.update_yaxes(title_text="Imaginary part Eps / Epb")
        # Plot the figure
        col1.plotly_chart(fig3, use_container_width=True)

        column_name3 = ['Eps_real', 'Epb_real']
        D = pd.DataFrame(list(zip(eps_real, epb_real)), x, 
                          columns = column_name3)
        # Creating the Figures
        fig4 = px.line(D, title='Real part of dielectric function')
        # Set x-axis title
        fig4.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig4.update_yaxes(title_text="Real part Eps / Epb")
        # Plot the figure
        col2.plotly_chart(fig4, use_container_width=True)
        
        # Download the generated data
        dfs = [A_trans, A_long, B, C, D]
        data = A.join(dfs)
        data = data.rename_axis("Wavelength (nm)")
        csv = data.to_csv().encode("utf-8")
    
        st.sidebar.download_button("Download data as CSV.csv",
                                    data = csv,
                                    file_name = "Simulated data.txt",
                                    mime="txt/csv")
    
    elif (modes_option == 'Average'):

        Results0 = smLSPR.size_modified_LSPR(x, R, AR_stdev, l, l_stdev, shape, mat, refind, 0, nPts, seed)
        Results1 = smLSPR.size_modified_LSPR(x, R, AR_stdev, l, l_stdev, shape, mat, refind, 1, nPts, seed)
        Results2 = smLSPR.size_modified_LSPR(x, R, AR_stdev, l, l_stdev, shape, mat, refind, 2, nPts, seed)
        Results3 = smLSPR.size_modified_LSPR(x, R, AR_stdev, l, l_stdev, shape, mat, refind, 3, nPts, seed)
        
        # Only longitudinal mode
        Ext_mode0 = Results0[0]
        Sca_mode0 = Results0[1]
        Abs_mode0 = Results0[2]
       
        # Only transversal mode
        Ext_3trans = (Results1[0]+Results2[0]+Results3[0])/3
        Sca_3trans = (Results1[1]+Results2[1]+Results3[1])/3
        Abs_3trans = (Results1[2]+Results2[2]+Results3[2])/3
        
        # Both combined
        Ext_average = (Results0[0]+Results1[0]+Results2[0]+Results3[0])/4
        Sca_average = (Results0[1]+Results1[1]+Results2[1]+Results3[1])/4
        Abs_average = (Results0[2]+Results1[2]+Results2[2]+Results3[2])/4
        
        QY_average = (Results0[3]+Results1[3]+Results2[3]+Results3[3])/4
        
        # Preparing the data for the LSPR
        column_name = ['Ext', 'Sca', 'Abs']
        A = pd.DataFrame(list(zip(Ext_average,Sca_average, Abs_average)), x, 
                         columns = column_name)
        
        # Transversal
        column_name_trans = ['Ext_trans', 'Sca_trans', 'Abs_trans']
        A_trans = pd.DataFrame(list(zip(Ext_3trans,Sca_3trans, Abs_3trans)), x, 
                         columns = column_name_trans)
        # Longitudinal 
        column_name_long = ['Ext_long', 'Sca_long', 'Abs_long']
        A_long = pd.DataFrame(list(zip(Ext_mode0,Sca_mode0, Abs_mode0)), x, 
                 columns = column_name_long)

        # Creating the Figures
        fig1 = px.line(A, title='Complete LSPR analytical simulation')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col1.plotly_chart(fig1, use_container_width=True)
        
        # Plotting the QY (%)
        
        # Preparing the data for the LSPR
        QY = ['QY']
        B = pd.DataFrame(QY_average*100, x, columns = QY)
        
        # Creating the Figures
        fig2 = px.line(B, title='Quantum Yield (%) - Scattering / Extinction cross section')
        # Set x-axis title
        fig2.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig2.update_yaxes(title_text="Quantum Yield (%)")
        # Plot the figure
        col2.plotly_chart(fig2, use_container_width=True)
        
        # Creating the Figures (Transversal part)
        fig1 = px.line(A_trans, title='Transversal (averaged) contribution to LSPR')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col1.plotly_chart(fig1, use_container_width=True)
        
        # Creating the Figures (Long. mode)
        fig1 = px.line(A_long, title='Longitudinal contribution to LSPR')
        # Set x-axis title
        fig1.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig1.update_yaxes(title_text="Extinction Cross section")
        # Plot the figure
        col2.plotly_chart(fig1, use_container_width=True)
        
        # We plot the dielectric functions Eps and Epb imag and real

        # Data for the dielectric functions of the material and the media
        # We separte real and imaginary parts
        eps_imag = np.imag(Results0[4])
        eps_real = np.real(Results0[4])

        epb_imag = np.imag(Results0[5])
        epb_real = np.real(Results0[5])

        column_name2 = ['Eps_imag', 'Epb_imag']
        C = pd.DataFrame(list(zip(eps_imag, epb_imag)), x, 
                          columns = column_name2)
         
        # Creating the Figures
        fig3 = px.line(C, title='Imaginary part of dielectric function')
        # Set x-axis title
        fig3.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig3.update_yaxes(title_text="Imaginary part Eps / Epb")
        # Plot the figure
        col1.plotly_chart(fig3, use_container_width=True)

        column_name3 = ['Eps_real', 'Epb_real']
        D = pd.DataFrame(list(zip(eps_real, epb_real)), x, 
                          columns = column_name3)
        # Creating the Figures
        fig4 = px.line(D, title='Real part of dielectric function')
        # Set x-axis title
        fig4.update_xaxes(title_text="Wavelength (nm)")
        # Set y-axes titles
        fig4.update_yaxes(title_text="Real part Eps / Epb")
        # Plot the figure
        col2.plotly_chart(fig4, use_container_width=True)
        
        # Download the generated data
        dfs = [A_trans, A_long, B, C, D]
        data = A.join(dfs)
        data = data.rename_axis("Wavelength (nm)")
        csv = data.to_csv().encode("utf-8")
    
        st.sidebar.download_button("Download data as CSV.csv",
                                    data = csv,
                                    file_name = "Simulated data.txt",
                                    mime="txt/csv")

# Creating the Figures
fig5 = px.histogram(Results[8], title='Size distribution of NP\'s length', marginal = 'violin')
# Set x-axis title
fig5.update_xaxes(title_text="Length (nm)", row = 1, col = 1)
# Set y-axes titles
fig5.update_yaxes(title_text="# of NPs with this length", row = 1, col = 1)
# Plot the figure
col1.plotly_chart(fig5, use_container_width=True)


# Creating the Figures
fig6 = px.histogram(Results[9], title='Size distribution of NP\'s AR', marginal = 'violin')
# Set x-axis title
fig6.update_xaxes(title_text="AR", row = 1, col = 1)
# Set y-axes titles
fig6.update_yaxes(title_text="# of NPs with this AR", row = 1, col = 1)
# Plot the figure
col2.plotly_chart(fig6, use_container_width=True)

st.write("""
          Shapes simulated with this widget
       """)
st.image("Shapes.jpg", width = 800)
