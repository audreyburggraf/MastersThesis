import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import seaborn as sns
import pandas as pd


import sys
import os

# Functions
# Define folder path
func_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FUNCTIONS/"

sys.path.append(func_folder_path)

from FITS_Image_Functions import *
from PolarizationFunctions import *
from PlottingWithFunction import *
from DataAnalysisFunctions import *
from GaussianFunctions import *

from custom_colormap import *

# Intro stuff
# distance_pc = 132


# Font sizes
# title_fs = 25
# axis_label_fs = 20
# axis_num_fs = 20
# legend_title_fs = 20
# legend_text_fs = 15
# cbar_fs = 20
# text_fs = 28


# ---------------------------------------------------------------------------------------------------------
# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")  # Replace with the actual path

# Now you can import constants.py
import constants

# Use the variable from constants.py
reference_length_AU = constants.reference_length_AU
distance_pc = constants.distance_pc
# ---------------------------------------------------------------------------------------------------------





# # Define folder path
folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND6/"

# # Define file paths
StokesI_file                = folder_path + "c2d_989_StokesI_233GHz.fits"
# StokesIerr_file             = folder_path + "c2d_989_StokesIerr_233GHz.fits"
# StokesQ_file                = folder_path + "c2d_989_StokesQ_233GHz.fits"
# StokesU_file                = folder_path + "c2d_989_StokesU_233GHz.fits"
# PolarizationAngle_file      = folder_path + "c2d_989_POLA_233GHz.fits"
# PolarizationAngle_err_file  = folder_path + "c2d_989_POLAerr_233GHz.fits"
# PolarizationFraction_file   = folder_path + "c2d_989_POLF_233GHz.fits"
# PolarizedIntensity_file     = folder_path + "c2d_989_POLI_233GHz.fits"
# PolarizedIntensity_err_file = folder_path + "c2d_989_POLIerr_233GHz.fits"

# # Load Stokes I data
StokesI_header, StokesI_data_4d_Jy, StokesI_data_2d_Jy, StokesI_wcs = read_in_file(StokesI_file)
# StokesI_data_2d_mJy = StokesI_data_2d_Jy * 1000

# # Define min/max for plotting
# StokesI_vmin = -0.2
# StokesI_vmax = 89.85

# # Clip and stretch data
# StokesI_clipped = np.clip(StokesI_data_2d_mJy, StokesI_vmin, StokesI_vmax)
# StokesI_stretched = stretch(StokesI_clipped, base=100, vmin=StokesI_vmin, vmax=StokesI_vmax)

# # Color bar ticks
# normalized_cbar_ticks = np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
# StokesI_unstretched_cbar_ticks = unstretch(normalized_cbar_ticks, StokesI_vmin, StokesI_vmax)

# # Load other FITS files
# _, _, StokesIerr_data_2d_Jy, StokesIerr_wcs = read_in_file(StokesIerr_file, dimensions=2)
# StokesIerr_data_2d_mJy = StokesIerr_data_2d_Jy * 1000

# StokesQ_header, _, StokesQ_data_2d_Jy, StokesQ_wcs = read_in_file(StokesQ_file)
# StokesQ_data_2d_mJy = StokesQ_data_2d_Jy * 1000

# StokesU_header, _, StokesU_data_2d_Jy, StokesU_wcs = read_in_file(StokesU_file)
# StokesU_data_2d_mJy = StokesU_data_2d_Jy * 1000

# PolarizationAngle_header, _, PolarizationAngle_data_2d_deg, PolarizationAngle_wcs = read_in_file(PolarizationAngle_file)
# PolarizationAngle_data_2d_rad = PolarizationAngle_data_2d_deg * (np.pi / 180)

# PolarizationAngle_err_header, _, PolarizationAngle_err_data_2d_deg, PolarizationAngle_err_wcs = read_in_file(PolarizationAngle_err_file)

# PolarizationFraction_header, _, PolarizationFraction_data_2d, PolarizationFraction_wcs = read_in_file(PolarizationFraction_file)

# PolarizedIntensity_header, _, PolarizedIntensity_data_2d_Jy, PolarizedIntensity_wcs = read_in_file(PolarizedIntensity_file)
# PolarizedIntensity_data_2d_mJy = PolarizedIntensity_data_2d_Jy * 1000

# PolarizedIntensity_err_header, _, PolarizedIntensity_err_data_2d_Jy, PolarizedIntensity_err_wcs = read_in_file(PolarizedIntensity_err_file)
# PolarizedIntensity_err_data_2d_mJy = PolarizedIntensity_err_data_2d_Jy * 1000

# Get beam info
beam_info = get_beam_info(StokesI_header)
BMAJ_deg, BMIN_deg, BMAJ_pix, BMIN_pix, BPA_deg_astronomy, BPA_deg_cartesian = beam_info

BPA_rad_astronomy = np.radians(BPA_deg_astronomy)
BPA_rad_cartesian = np.radians(BPA_deg_cartesian)
# print(f'The angle in astronomy angles is {BPA_astronomy_deg:.3f} degrees')
# print(f'The angle in cartesian angles is {BPA_deg_cartesian:.3f} degrees')

# Compute reference length
# reference_length_AU = 100
reference_length_pix = length_in_pixels(reference_length_AU, distance_pc, StokesI_header)
# print(f'The reference length is: {reference_length_AU} AU ({reference_length_pix:.3f} pix)')

# Compute polarization properties
# calculated_polarization_angle_rad = calculate_polarization_angle(StokesQ_data_2d_Jy, StokesU_data_2d_Jy)
# calc_polarized_frac = calculate_polarized_fraction(StokesQ_data_2d_mJy, StokesU_data_2d_mJy, StokesI_data_2d_mJy)
# calculated_polarized_intensity = calculate_polarized_intensity(StokesQ_data_2d_mJy, StokesU_data_2d_mJy)

# Find pixel values for sky positions
# min_str = 'J163135.570-240131.200'
# max_str = 'J163135.750-240128.700'
# centre_str = 'J163135.657-240129.935'

RA_min_pix, Dec_min_pix       = string_to_pixel(constants.min_str_band6, StokesI_wcs)
RA_max_pix, Dec_max_pix       = string_to_pixel(constants.max_str_band6, StokesI_wcs)
RA_centre_pix, Dec_centre_pix = string_to_pixel(constants.centre_str_band6, StokesI_wcs)
# print(f'StokesI_min_RA_pix: {RA_min_pix}, StokesI_min_Dec_pix: {Dec_min_pix}')
# print(f'StokesI_max_RA_pix: {RA_max_pix}, StokesI_max_Dec_pix: {Dec_max_pix}')

# Define plotting limits
xmin, xmax = RA_max_pix, RA_min_pix
ymin, ymax = Dec_min_pix, Dec_max_pix



# Minor and major angles
# --------------------------------------------------------------------------------------------------  
# PA_deg = 148

# # The positional angle is the major axis 
# major_angle_deg = PA_deg

# # The minor axis is 90 degrees from the major
# minor_angle_deg = major_angle_deg - 90
# minor_angle_rad_astronomy = np.radians(minor_angle_deg)

# # We have to offset the minor and major angle by 90 degrees if we want to plot it:
# major_angle_deg_cartesian = major_angle_deg - 90
# minor_angle_deg_cartesian = minor_angle_deg - 90

# # We can turn the plotted angles into radians
# minor_angle_rad_cartesian = np.radians(minor_angle_deg_cartesian)
# major_angle_rad_cartesian = np.radians(major_angle_deg_cartesian)
# --------------------------------------------------------------------------------------------------  





# General vector stuff
# --------------------------------------------------------------------------------------------------  
# vector_length_pix_const = 4

# max_length_pix = 400  # Maximum length of the vector in pixels for 100% polarization
# reference_fraction = 0.03
nx, ny = StokesI_data_2d_mJy.shape
# step = 7 
# --------------------------------------------------------------------------------------------------  


# Actual vectors with scaled length
# --------------------------------------------------------------------------------------------------  
# List to store vector data
vector_data_actual_cartesian = []

# Loop over values in x and y
for x in range(0, nx, step):
    for y in range(0, ny, step):
        # Check if the conditions are met
        if (StokesI_data_2d_mJy[y, x] / StokesIerr_data_2d_mJy[y, x] > 3 and 
            calculated_polarized_intensity[y, x] / PolarizedIntensity_err_data_2d_mJy[y, x] > 3 and 
            PolarizationAngle_err_data_2d_deg[y, x] < 10): 
            
            # Get the polarization angle at this pixel
            angle_rad_cartesian = calculated_polarization_angle_rad[y, x] + np.pi/2
            
            # Get the polarization fraction at this pixel
            polarization_fraction = calc_polarized_frac[y, x]  # Add your data source
            
            # Scale vector length by polarization fraction
            vector_length_pix = max_length_pix * polarization_fraction
            
            # Compute the vector components
            dx = vector_length_pix * np.cos(angle_rad_cartesian)
            dy = vector_length_pix * np.sin(angle_rad_cartesian)
            
            # Store vector data in a row
            vector_data_actual_cartesian.append([x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2])
# --------------------------------------------------------------------------------------------------              
            


# Actual vectors with the same length
# --------------------------------------------------------------------------------------------------      
# List to store vector data
vector_data_actual_same_length_cartesian = []
vector_angle_actual_same_length_astronomy = []

# Loop over values in x and y
for x in range(0, nx, step):
    for y in range(0, ny, step):
        # Check if the conditions are met
        if (StokesI_data_2d_mJy[y, x] / StokesIerr_data_2d_mJy[y, x] > 3 and 
            calculated_polarized_intensity[y, x] / PolarizedIntensity_err_data_2d_mJy[y, x] > 3 and 
            PolarizationAngle_err_data_2d_deg[y, x] < 10): 
            
            # Get the polarization angle at this pixel
            
            angle_rad_cartesian = calculated_polarization_angle_rad[y, x] + np.pi / 2
            
            # Compute the vector components
            dx = vector_length_pix_const * np.cos(angle_rad_cartesian)
            dy = vector_length_pix_const * np.sin(angle_rad_cartesian)
            
            # Store vector data in a row
            vector_data_actual_same_length_cartesian.append([x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2])
            vector_angle_actual_same_length_astronomy.append(calculated_polarization_angle_rad[y, x])
# --------------------------------------------------------------------------------------------------              
            
 # Position angle grids for uniform and azimuthal
# --------------------------------------------------------------------------------------------------  
# PA_grid_100Uniform   = np.zeros((ny, nx))
# PA_grid_100Azimuthal = np.zeros((ny, nx))


# # Loop over values in x and y
# for x in range(0, nx):
#     for y in range(0, ny): 
        
#         PA_grid_100Uniform[y, x] = minor_angle_rad_astronomy
        
#         # Calculate the azimuthal angle (in radians) based on the position
#         dx_center = x - RA_centre_pix  # Center of the disk (or image)
#         dy_center = y - Dec_centre_pix  # Center of the disk (or image)
            
#         # want tangent to radial but in sky coordates (+ 90 and - 90 cancel)
#         azimuthal_angle_rad_astronomy = np.arctan2(dy_center, dx_center)  # Calculate the angle
        
#         # azimuthal_angle_rad_astronomy = castesian_to_astronomy_rad(azimuthal_angle_rad_cartesian)
        
#         PA_grid_100Azimuthal[y, x] = azimuthal_angle_rad_astronomy
    
    
PA_grid_100Uniform = make_PA_grid_100Uniform(ny, nx, minor_angle_rad_sky_band6)
PA_grid_100Azimuthal = make_PA_grid_100Azimuthal(ny, nx, RA_centre_pix, Dec_centre_pix)        
        
        
        
# --------------------------------------------------------------------------------------------------              
           
            
# Uniform vectors
# --------------------------------------------------------------------------------------------------  
# List to store vector data
vector_data_uniform_cartesian  = []
vector_angle_uniform_astronomy = []

# Loop over values in x and y
for x in range(0, nx, step):
    for y in range(0, ny, step):
        # Check if the conditions are met
        if (StokesI_data_2d_mJy[y, x] / StokesIerr_data_2d_mJy[y, x] > 3 and 
            calculated_polarized_intensity[y, x] / PolarizedIntensity_err_data_2d_mJy[y, x] > 3 and 
            PolarizationAngle_err_data_2d_deg[y, x] < 10):
            
            uniform_angle_rad_astronomy = PA_grid_100Uniform[y,x]
            
            # Compute the vector components
            dx = vector_length_pix_const * np.cos(uniform_angle_rad_astronomy + np.pi/2)
            dy = vector_length_pix_const * np.sin(uniform_angle_rad_astronomy + np.pi/2)
            
            # Store vector data in a row
            vector_data_uniform_cartesian.append([x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2])
            vector_angle_uniform_astronomy.append(uniform_angle_rad_astronomy)
# --------------------------------------------------------------------------------------------------              



# Azimuthal vectors 
# --------------------------------------------------------------------------------------------------                       
# List to store vector data for azimuthal polarization
vector_data_azimuthal_cartesian = []
vector_angle_azimuthal_astronomy = []

# Loop over values in x and y
for x in range(0, nx, step):
    for y in range(0, ny, step):
        # Check if the conditions are met
        if (StokesI_data_2d_mJy[y, x] / StokesIerr_data_2d_mJy[y, x] > 3 and 
            calculated_polarized_intensity[y, x] / PolarizedIntensity_err_data_2d_mJy[y, x] > 3 and 
            PolarizationAngle_err_data_2d_deg[y, x] < 10): 
            
            azimuthal_angle_rad_astronomy = PA_grid_100Azimuthal[y,x]

            # Compute the vector components using the azimuthal angle
            dx = vector_length_pix_const * np.cos(azimuthal_angle_rad_astronomy + np.pi/2)
            dy = vector_length_pix_const * np.sin(azimuthal_angle_rad_astronomy + np.pi/2)
            
            # Store vector data in a row
            vector_data_azimuthal_cartesian.append([x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2])
            
            vector_angle_azimuthal_astronomy.append(azimuthal_angle_rad_astronomy)
# --------------------------------------------------------------------------------------------------                       
 
StokesQ_grid_100Uniform,   StokesU_grid_100Uniform   = recover_StokesQU(PA_grid_100Uniform,   StokesI_data_2d_mJy, ny, nx)
StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal = recover_StokesQU(PA_grid_100Azimuthal, StokesI_data_2d_mJy, ny, nx)


TestingQU_angles_labels = ['0 deg', '45 deg', '90 deg', '135 deg']
TestingQU_angles_deg_ast = [0, 45, 90, 135]
TestingQU_angles_rad_ast = np.radians(TestingQU_angles_deg_ast)