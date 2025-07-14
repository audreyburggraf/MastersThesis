import numpy as np 
import seaborn as sns
import sys
import pandas as pd



# Font sizes
title_fs = 25
axis_label_fs = 20
axis_num_fs = 20
legend_title_fs = 20
legend_text_fs = 15
cbar_fs = 20
text_fs = 28
cbar_num_fs = 20



normalized_cbar_ticks = np.array([0, 0.2, 0.4, 0.6, 0.8, 1])


# ALMA Band wavelengths
# ---------------------------------------------
# lambda_band1_mm  = 7.0
# lambda_band1_mm  = 3.1
# lambda_band4_mm  = 2.0
# lambda_band6_mm  = 1.3
# lambda_band7_mm  = 0.87
# lambda_band10_mm = 0.34


lambda_mm = {
    1: 7.0,
    3: 3.1,
    4: 2.0,
    6: 1.3,
    7: 0.87,
    10: 0.34
}

# lambda_band4_micron  = 2.0  * 1000  # mm to microns
# lambda_band6_micron  = 1.3  * 1000  # mm to microns
# lambda_band7_micron  = 0.87 * 1000  # mm to microns
# ---------------------------------------------



# Constants for the disk
distance_pc = 132


# File locations
band4_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND4/"
band6_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND6/"
band7_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND7/"
band7_v0_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND7_v0/"
band7_redo_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND7_redo/"
band7_nterms2_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND7_nterms2/"
band7_uncalibrated_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND7_uncalibrated/"

band4_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND4/"
band6_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND6/"
band7_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND7/"
band7_redo_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND7_redo/"
band7_nterms2_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND7_nterms2/"
band7_uncalibrated_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND7_uncalibrated/"


functions_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FUNCTIONS/"

image_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/WRITEUP_AND_IMAGES/IMAGES/"

writeup_image_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/WRITEUP_AND_IMAGES/IMAGES/WRITEUP/"
slideshow_image_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/WRITEUP_AND_IMAGES/IMAGES/SLIDESHOW/"
poster_image_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/WRITEUP_AND_IMAGES/IMAGES/POSTER/"

m_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/M FILES/"

# Load in necessary functions
# ---------------------------------------------
sys.path.append(functions_folder_path)

from IntroductionFunctions import *
from FITS_Image_Functions import * 
# ---------------------------------------------


# Minimum, maximum and centre constants
# ---------------------------------------------
# Band 6 
# ---------------------------------------------
min_str_band6    = 'J163135.570-240131.200'
max_str_band6    = 'J163135.750-240128.700'
centre_str_band6 = 'J163135.657-240129.935'
# ---------------------------------------------
# Band 4
# ---------------------------------------------
# min_str_band4    = 'J163135.570-240131.200'                  
# max_str_band4    = 'J163135.750-240128.700'
centre_str_band4 = 'J163135.656-240129.992' # This value is from CASA            
# ---------------------------------------------
# Band 7
# ---------------------------------------------
centre_str_band7 = 'J163135.656-240130.089' # This value is from CASA
# ---------------------------------------------



# Position, Minor and Major angles
# -------------------------------------------------------------------------------------------------------- 
# Function
# -------------------------------------------------------------------------------------------------------- 
def get_minor_major_angles(PA_deg_sky):
    
    # The major angle in sky coordinates is the position angle
    major_angle_deg_sky = PA_deg_sky
    
    # The minor angle in sky coordinates is 90 degrees less then major
    minor_angle_deg_sky = major_angle_deg_sky - 90
    
    # Convert from sky to cartesian corrdinates 
    major_angle_deg_cartesian = astronomy_to_cartesian(major_angle_deg_sky)
    minor_angle_deg_cartesian = astronomy_to_cartesian(minor_angle_deg_sky)
    
    # Convert from degrees to radians
    minor_angle_rad_sky = np.radians(minor_angle_deg_sky)
    major_angle_rad_sky = np.radians(major_angle_deg_sky)
    
    minor_angle_rad_cartesian = np.radians(minor_angle_deg_cartesian)
    major_angle_rad_cartesian = np.radians(major_angle_deg_cartesian)
    
    return minor_angle_rad_sky, major_angle_rad_sky, minor_angle_rad_cartesian, major_angle_rad_cartesian
# --------------------------------------------------------------------------------------------------------    
# Band 6 
# ---------------------------------------------
PA_deg_sky_band6 = 148

minor_angle_rad_sky_band6, major_angle_rad_sky_band6, minor_angle_rad_cartesian_band6, major_angle_rad_cartesian_band6 = get_minor_major_angles(PA_deg_sky_band6)
# ---------------------------------------------
# Band 4
# ---------------------------------------------
PA_deg_sky_band4 = 145.9 # This value is frmom CASA 

minor_angle_rad_sky_band4, major_angle_rad_sky_band4, minor_angle_rad_cartesian_band4, major_angle_rad_cartesian_band4 = get_minor_major_angles(PA_deg_sky_band4)
# ---------------------------------------------


# Band 7
# ---------------------------------------------
PA_deg_sky_band7 = 139.49 # this is from fittiting on casa Band 7 redo, p5

minor_angle_rad_sky_band7, major_angle_rad_sky_band7, minor_angle_rad_cartesian_band7, major_angle_rad_cartesian_band7 = get_minor_major_angles(PA_deg_sky_band7)
# ---------------------------------------------



# Rms errors
# ------------------------------------------------------------------------------------------
# Band 7
# ---------------------------------------------
# Get these values using STDEV on CARTA
StokesQ_err_mJy_band7_uncalibrated = 0.03
StokesU_err_mJy_band7_uncalibrated = 0.03
POLI_err_mJy_band7_uncalibrated = 0.05 
# ------------------------------------------------------------------------------------------
StokesI_err_mJy_band7_v0 = np.mean([0.0646, 0.0580, 0.0636, 0.056])
StokesQ_err_mJy_band7_v0 = np.mean([0.0317, 0.0311, 0.0398, 0.034])
StokesU_err_mJy_band7_v0 = np.mean([0.0247, 0.0249, 0.0255, 0.024])
# ------------------------------------------------------------------------------------------
# StokesI_err_mJy_band7_redo = np.mean([0.0398, 0.0599, 0.0370, 0.0478])
# StokesQ_err_mJy_band7_redo = np.mean([])
# StokesU_err_mJy_band7_redo = np.mean([])
# ------------------------------------------------------------------------------------------
StokesI_err_mJy_band7_nterms2 = np.mean([0.0318, 0.0314, 0.0380, 0.0285])
StokesQ_err_mJy_band7_nterms2 = np.mean([0.0342, 0.0372, 0.0321, 0.0363])
StokesU_err_mJy_band7_nterms2 = np.mean([0.0263, 0.0270, 0.0271, 0.0284])
# ------------------------------------------------------------------------------------------


# Vector Constants
# ------------------------------------------------------------------------------------------
reference_length_AU = 100
max_length_pix = 400  # Maximum length of the vector in pixels for 100% polarization
reference_fraction = 0.03
step_band6 = 6
step_band4 = 4
step_band7 = 6

vector_length_pix_const = 0 
vector_len_pix_band6 = 4
vector_len_pix_band4 = 2
vector_len_pix_band7 = 4
# ------------------------------------------------------------------------------------------




# Ratio Stuff
# ------------------------------------------------------------------------------------------
testing_ratios = [(1, 0), 
                  (0.9, 0.1), 
                  (0.8, 0.2), 
                  (0.7, 0.3), 
                  (0.6, 0.4), 
                  (0.5, 0.5), 
                  (0.4, 0.6), 
                  (0.3, 0.7), 
                  (0.2, 0.8), 
                  (0.1, 0.9), 
                  (0, 1)]
# ------------------------------------------------------------------------------------------