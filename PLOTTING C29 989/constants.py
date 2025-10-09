import numpy as np 
import seaborn as sns
import sys
import pandas as pd
from scipy.optimize import minimize


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
lambda_mm = {
    1: 7.0,
    3: 3.1,
    4: 2.0,
    5: 1.6, # check this 
    6: 1.3,
    7: 0.87,
    10: 0.34
}
# ---------------------------------------------



# Constants for the disk
distance_pc = 132


# File locations
# -------------------------------------------------------------------------------------------------------- 
# FITS Files
# -------------------------------------------------------------------------------------------------------- 
band4_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND4/"
band5_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND5/"
band6_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND6/"
band7_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND7/"
# band7_v0_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND7_v0/"
# band7_redo_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND7_redo/"
# band7_nterms2_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND7_nterms2/"
# band7_uncalibrated_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND7_uncalibrated/"
# -------------------------------------------------------------------------------------------------------- 
# CARTA
# -------------------------------------------------------------------------------------------------------- 
band4_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND4/"
band5_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND5/"
band6_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND6/"
band7_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND7/"
# band7_redo_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND7_redo/"
# band7_nterms2_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND7_nterms2/"
# band7_uncalibrated_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND7_uncalibrated/"
# -------------------------------------------------------------------------------------------------------- 
# Other
# -------------------------------------------------------------------------------------------------------- 
functions_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FUNCTIONS/"

image_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/WRITEUP_AND_IMAGES/IMAGES/"

writeup_image_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/WRITEUP_AND_IMAGES/IMAGES/WRITEUP/"
slideshow_image_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/WRITEUP_AND_IMAGES/IMAGES/SLIDESHOW/"
poster_image_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/WRITEUP_AND_IMAGES/IMAGES/POSTER/"

m_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/M FILES/"

P_omega_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/DUST MODEL NOTEBOOKS/P_omega_Data/"
# -------------------------------------------------------------------------------------------------------- 





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
# Band 5
# ---------------------------------------------
centre_str_band5 = 'J163135.656-240130.086' # This value is from CASA 
# ---------------------------------------------
# Band 4
# ---------------------------------------------
centre_str_band4 = 'J163135.656-240129.992' # This value is from CASA            
# ---------------------------------------------
# Band 7
# ---------------------------------------------
centre_str_band7 = 'J163135.656-240130.089' # This value is from CASA and was updated to new file
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



# ---------------------------------------------
# Band 6 
# ---------------------------------------------
PA_deg_sky_band6 = 148

minor_angle_rad_sky_band6, major_angle_rad_sky_band6, minor_angle_rad_cartesian_band6, major_angle_rad_cartesian_band6 = get_minor_major_angles(PA_deg_sky_band6)
# ---------------------------------------------
# Band 5
# ---------------------------------------------
PA_deg_sky_band5 = 140.25 # This value is from CASA ****** not updated to band 5(still?)

minor_angle_rad_sky_band5, major_angle_rad_sky_band5, minor_angle_rad_cartesian_band5, major_angle_rad_cartesian_band5 = get_minor_major_angles(PA_deg_sky_band5)
# ---------------------------------------------
# Band 4
# ---------------------------------------------
PA_deg_sky_band4 = 145.9 # This value is from CASA 

minor_angle_rad_sky_band4, major_angle_rad_sky_band4, minor_angle_rad_cartesian_band4, major_angle_rad_cartesian_band4 = get_minor_major_angles(PA_deg_sky_band4)
# ---------------------------------------------

# Band 7
# ---------------------------------------------
PA_deg_sky_band7 = 139.65 # This value is from fitting the new band 7 file (image component size, or 148.86 deconvolved)

minor_angle_rad_sky_band7, major_angle_rad_sky_band7, minor_angle_rad_cartesian_band7, major_angle_rad_cartesian_band7 = get_minor_major_angles(PA_deg_sky_band7)
# ---------------------------------------------



# Rms errors
# ------------------------------------------------------------------------------------------
# --------
# Band 5
# --------
StokesI_err_mJy_band5 = np.mean([0.0135, 0.0147, 0.0149, 0.0144])
StokesQ_err_mJy_band5 = np.mean([0.0127, 0.0137, 0.0129, 0.0128]) 
StokesU_err_mJy_band5 = np.mean([0.0150, 0.0155, 0.0157, 0.0148]) 
POLI_err_mJy_band5 = np.mean([StokesQ_err_mJy_band5, StokesU_err_mJy_band5])
# ------------------------------------------------------------------------------------------
# Band 7
# ---------------------------------------------
# Get these values using STDEV on CARTA
# StokesQ_err_mJy_band7_uncalibrated = 0.03
# StokesU_err_mJy_band7_uncalibrated = 0.03
# POLI_err_mJy_band7_uncalibrated = 0.05 
# ------------------------------------------------------------------------------------------
StokesI_err_mJy_band7 = np.mean([0.0340, 0.0382, 0.0375, 0.0310])
StokesQ_err_mJy_band7 = np.mean([0.0256, 0.0266, 0.0268, 0.0260]) 
StokesU_err_mJy_band7 = np.mean([0.0257, 0.0243, 0.0234, 0.0240]) 
POLI_err_mJy_band7 = np.mean([StokesQ_err_mJy_band7, StokesU_err_mJy_band7])
# ------------------------------------------------------------------------------------------

# Vector Constants
# ------------------------------------------------------------------------------------------
reference_length_AU = 100
max_length_pix = 400  # Maximum length of the vector in pixels for 100% polarization
reference_fraction = 0.03

step_band4 = 4
step_band5 = 8 
step_band6 = 6  
step_band7 = 6

vector_length_pix_const = 0 
vector_len_pix_band4 = 2
vector_len_pix_band5 = 6  
vector_len_pix_band6 = 4
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
