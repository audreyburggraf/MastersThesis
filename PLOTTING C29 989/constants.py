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
    # Band 1
    1: 7.0,
    # Band 3
    3: 3.1,
    # Band 4
    4: 2.0,
    42: 2.0,
    "Band 4": 2.0,
    "Band 4 nterms2": 2.0,
    "Band 4 nterms2 smooth": 2.0,
    "Band 4 nterms2 smooth B6": 2.0,
    "Band 4 nterms2 smooth B6 B7": 2.0,
    # Band 5
    5: 1.7, 
    50: 1.7,
    52: 1.7,
    "Band 5": 1.7,
    "Band 5 robust -2": 1.7, 
    "Band 5 robust -1": 1.7, 
    # Band 6
    6: 1.3,
    "Band 6": 1.3,
    "Band 6 smooth": 1.3,
    "Band 6 smooth B7": 1.3,
    # Band 7
    7: 0.87, 
    72: 0.87, 
    "Band 7": 0.87, 
    "Band 7 nterms2": 0.87, 
    "Band 7 nterms2 smooth": 0.87, 
    "Band 7 nterms2 smooth B6": 0.87, 
    # Band 10
    10: 0.34
}
# ---------------------------------------------
# Setting the colours for each band 
# The colour will be the same if its nterms = 1 or nterms = 2
# mw_colors = {
#     4:  "#ff1900",  # Red
#     42:  "#ff1900",  # Red
#     5: "#ffaf1b",  # Soft Orange
#     6: "#8ddafe",  # Sky Blue
#     7: "#c5b3fd",   # Soft Purple
#     72: "#c5b3fd"   # Soft Purple
# }


# Set the colours for each band 
alma_band_colors = {
    1: "#f63fe3",  # pink from Kataoka et al. 2015
    3: "#2321f6",  # blue from Kataoka et al. 2015
   
    #--------
    # Band 4
    #--------
    4: "red", 
    42: "red", 
    "Band 4": "red",
    "Band 4 nterms2": 'red',
    "Band 4 nterms2 smooth": 'red',
    "Band 4 nterms2 smooth B6": 'red',
    "Band 4 nterms2 smooth B6 B7": 'red',
    #--------
    # Band 5
    #--------
    5: "orange",  
    50: "orange",  
    52: "orange",  
    "Band 5": "orange",
    "Band 5 v0": "orange",
    "Band 5 nterms 2": "orange",
    "Band 5 robust -2": "orange",
    "Band 5 robust -1": "orange",
    #--------
    # Band 6
    #--------
    6: "#00bfc4", # Bright cyan
    "Band 6": "#00bfc4",
    "Band 6 smooth": "#00bfc4",
    "Band 6 smooth B7": "#00bfc4",
    #--------
    # Band 7
    #--------
    7: "darkviolet", 
    72: "darkviolet", 
    "Band 7":  "darkviolet",
    "Band 7 nterms2": "darkviolet",
    "Band 7 nterms2 smooth": "darkviolet",
    "Band 7 nterms2 smooth B6": "darkviolet",
    10:"#fb0d0d", # red from Kataoka et al. 2015
}







# Labels
# alma_quad_labels_mm = [lambda_mm[4], lambda_mm[5], lambda_mm[6], lambda_mm[7]]
# alma_quad_labels_band = ["Band 4", "Band 5", "Band 6", "Band 7"]
# alma_quad_labels = alma_quad_labels_mm
# ---------------------------------------------
# ---------------------------------------------



# Constants for the disk
distance_pc = 132


# File locations
# -------------------------------------------------------------------------------------------------------- 
base_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/"

plotting_path = base_path + "PLOTTING C29 989/"

writeup_path = base_path + "WRITEUP_AND_IMAGES/"

image_folder_path = writeup_path + "IMAGES/"


# -------------------------------------------------------------------------------------------------------- 
# FITS Files
# -------------------------------------------------------------------------------------------------------- 
fits_path = plotting_path + "FITS FILES/"

band4_data_folder_path = fits_path + "BAND4/"
band4_nterms2_data_folder_path = fits_path + "BAND4_nterms2/"
band4_nterms2_smooth_data_folder_path = fits_path + "BAND4_nterms2_smooth/"
band4_nterms2_smooth_B6_data_folder_path = fits_path + "BAND4_nterms2_smooth_B6/"
band4_nterms2_smooth_B6_B7_data_folder_path = fits_path + "BAND4_nterms2_smooth_B6_B7/"

band5_data_folder_path = fits_path + "BAND5/"
band5_v0_data_folder_path = fits_path + "BAND5_v0/"
band5_nterms2_data_folder_path = fits_path + "BAND5_nterms2/"
band5_robust_minus2_data_folder_path = fits_path + "BAND5_robust_minus2/"
band5_robust_minus1_data_folder_path = fits_path + "BAND5_robust_minus1/"

band6_data_folder_path = fits_path + "BAND6/"
band6_smooth_data_folder_path = fits_path + "BAND6_smooth/"
band6_smooth_B7_data_folder_path = fits_path + "BAND6_smooth_B7/"

band7_nterms2_data_folder_path = fits_path + "BAND7_nterms2/"
band7_nterms2_smooth_data_folder_path = fits_path + "BAND7_nterms2_smooth/"
band7_nterms2_smooth_B6_data_folder_path = fits_path + "BAND7_nterms2_smooth_B6/"
# -------------------------------------------------------------------------------------------------------- 
# CARTA
# -------------------------------------------------------------------------------------------------------- 
carta_path = base_path + "CARTA FILES/"

band4_carta_folder_path = carta_path + "BAND4/"
band4_nterms2_carta_folder_path = carta_path + "BAND4_nterms2/"
band4_nterms2_smooth_carta_folder_path = carta_path + "BAND4_nterms2_smooth/"
band4_nterms2_smooth_carta_folder_path = carta_path + "BAND4_nterms2_smooth_B6/"
band4_nterms2_smooth_carta_folder_path = carta_path + "BAND4_nterms2_smooth_B7/"

band5_carta_folder_path = carta_path + "BAND5/"
band5_v0_carta_folder_path = carta_path + "BAND5_v0/"
band5_nterms2_carta_folder_path = carta_path + "BAND5_nterms2/"
band5_robust_minus2_carta_folder_path = carta_path + "BAND5_robust_minus2/"
band5_robust_minus1_carta_folder_path = carta_path + "BAND5_robust_minus1/"

band6_carta_folder_path = carta_path + "BAND6/"
band6_smooth_carta_folder_path = carta_path + "BAND6_smooth/"
band6_smooth_carta_folder_path = carta_path + "BAND6_smooth_B6/"

band7_nterms2_carta_folder_path = carta_path + "BAND7_nterms2/"
band7_nterms2_smooth_carta_folder_path = carta_path + "BAND7_nterms2_smooth/"
band7_nterms2_smooth_carta_folder_path = carta_path + "BAND7_nterms2_smooth_B6/"
# -------------------------------------------------------------------------------------------------------- 
# Other
# -------------------------------------------------------------------------------------------------------- 

functions_folder_path = plotting_path + "FUNCTIONS/"
m_data_folder_path = plotting_path + "M FILES/"
P_omega_data_folder_path = plotting_path + "DUST MODEL NOTEBOOKS/P_omega_Data/"
# ------------------------------------------------------------------
# Image folders
# ------------------------------------------------------------------

writeup_image_folder_path = image_folder_path + "WRITEUP/"
slideshow_image_folder_path = image_folder_path + "SLIDESHOW/"
poster_image_folder_path = image_folder_path + "POSTER/"
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


# ------------------------------------------------------------------------------------------
# Image centre coordinates
# ------------------------------------------------------------------------------------------

centre_data = {

    # Band 6
    "band6": 'J163135.657-240129.935',
    "band6_smooth": 'J163135.657-240129.939',
    "band6_smooth_B7": None,

    # Band 5
    "band5_v0": 'J163135.656-240130.086',
    "band5": 'J163135.656-240130.086',
    "band5_nterms2": 'J163135.656-240130.086',
    "band5_robust_minus2": 'J163135.656-240130.086',
    "band5_robust_minus1": 'J163135.656-240130.086',  # update later

    # Band 4
    "band4": 'J163135.656-240129.992',
    "band4_nterms2": 'J163135.656-240130.006',
    "band4_nterms2_smooth": 'J163135.656-240130.004',
    "band4_nterms2_smooth_B6": 'J163135.656-240130.005',
    "band4_nterms2_smooth_B6_B7": 'J163135.656-240130.005',

    # Band 7
    "band7_nterms2": 'J163135.656-240130.089',
    "band7_nterms2_smooth": None,
    "band7_nterms2_smooth_B6": None,
}

# ------------------------------------------------------------------------------------------
# Create variables automatically (keeps compatibility with existing scripts)
# ------------------------------------------------------------------------------------------

for band, value in centre_data.items():
    globals()[f"centre_str_{band}"] = value
    
    
# centre_str_band6 = 'J163135.657-240129.935'
# centre_str_band6_smooth = 'J163135.657-240129.939'
# centre_str_band6_smooth_B7 = 
# # ---------------------------------------------
# # Band 5
# # ---------------------------------------------
# centre_str_band5_v0 = 'J163135.656-240130.086' # This value is from CASA 
# centre_str_band5 = 'J163135.656-240130.086'
# centre_str_band5_nterms2 = 'J163135.656-240130.086'
# centre_str_band5_robust_minus2 = 'J163135.656-240130.086'
# centre_str_band5_robust_minus1 = 'J163135.656-240130.086' # update this
# # ---------------------------------------------
# # Band 4
# # ---------------------------------------------
# centre_str_band4 = 'J163135.656-240129.992' # This value is from CASA            
# # ---------------------------------------------
# # Band 4 (nterms = 2)
# # ---------------------------------------------
# centre_str_band4_nterms2 = 'J163135.656-240130.006' # This is from CASA and has been updated for nterms 2
# centre_str_band4_nterms2_smooth = 'J163135.656-240130.004'
# centre_str_band4_nterms2_smooth_B6 = 
# centre_str_band4_nterms2_smooth_B6_B7 = 
# # ---------------------------------------------
# # Band 7
# # ---------------------------------------------
# # I think this might be nterms = 2
# # centre_str_band7 = 'J163135.656-240130.089' # This value is from CASA and was updated to new file
# # ---------------------------------------------
# # This value is correct for nterms = 2
# centre_str_band7_nterms2 = 'J163135.656-240130.089'
# centre_str_band7_nterms2_smooth = 
# centre_str_band7_nterms2_smooth_B6 = 
# # ---------------------------------------------


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
PA_values = {

    # Band 6
    "band6": 148,
    "band6_smooth": 149.13,
    "band6_smooth_B7": 0,

    # Band 5
    "band5_v0": 140.25,  # CASA value (may update)
    "band5": 148.79,
    "band5_nterms2": 148.83,
    "band5_robust_minus2": 148.41,
    "band5_robust_minus1": 148.41,

    # Band 4
    "band4": 145.9,  # CASA value
    "band4_nterms2": 139.31,  # updated for nterms = 2
    "band4_nterms2_smooth": 148.30,
    "band4_nterms2_smooth_B6": 147.23, 
    "band4_nterms2_smooth_B6_B7": 147.38,

    # Band 7
    "band7_nterms2": 139.65,
    "band7_nterms2_smooth": 0,
     "band7_nterms2_smooth_B6": 0,
}

# --------------------------------------------------------------------------------------------------------
# Automatically compute angles and create the SAME variable names as before
# --------------------------------------------------------------------------------------------------------

for band, PA in PA_values.items():

    minor_sky, major_sky, minor_cart, major_cart = get_minor_major_angles(PA)

    globals()[f"minor_angle_rad_sky_{band}"] = minor_sky
    globals()[f"major_angle_rad_sky_{band}"] = major_sky
    globals()[f"minor_angle_rad_cartesian_{band}"] = minor_cart
    globals()[f"major_angle_rad_cartesian_{band}"] = major_cart


# # This is the older way of doing things, keep it here for now:
# # ---------------------------------------------
# # Band 6 
# # ---------------------------------------------
# PA_deg_sky_band6 = 148

# minor_angle_rad_sky_band6, major_angle_rad_sky_band6, minor_angle_rad_cartesian_band6, major_angle_rad_cartesian_band6 = get_minor_major_angles(PA_deg_sky_band6)
# # ---------------------------------------------
# PA_deg_sky_band6_smooth = 149.13

# minor_angle_rad_sky_band6_smooth, major_angle_rad_sky_band6_smooth, minor_angle_rad_cartesian_band6_smooth, major_angle_rad_cartesian_band6_smooth = get_minor_major_angles(PA_deg_sky_band6_smooth)
# # ---------------------------------------------

# # Band 5
# # ---------------------------------------------
# PA_deg_sky_band5_v0 = 140.25 # This value is from CASA ****** not updated to band 5(still?)

# minor_angle_rad_sky_band5_v0, major_angle_rad_sky_band5_v0, minor_angle_rad_cartesian_band5_v0, major_angle_rad_cartesian_band5_v0 = get_minor_major_angles(PA_deg_sky_band5_v0)
# # ---------------------------------------------
# PA_deg_sky_band5 = 148.79 

# minor_angle_rad_sky_band5, major_angle_rad_sky_band5, minor_angle_rad_cartesian_band5, major_angle_rad_cartesian_band5 = get_minor_major_angles(PA_deg_sky_band5)
# # ---------------------------------------------
# PA_deg_sky_band5_nterms2 = 148.83

# minor_angle_rad_sky_band5_nterms2, major_angle_rad_sky_band5_nterms2, minor_angle_rad_cartesian_band5_nterms2, major_angle_rad_cartesian_band5_nterms2 = get_minor_major_angles(PA_deg_sky_band5_nterms2)
# # ---------------------------------------------
# PA_deg_sky_band5_robust_minus2 = 148.41 

# minor_angle_rad_sky_band5_robust_minus2, major_angle_rad_sky_band5_robust_minus2, minor_angle_rad_cartesian_band5_robust_minus2, major_angle_rad_cartesian_band5_robust_minus2 = get_minor_major_angles(PA_deg_sky_band5_robust_minus2)
# # ---------------------------------------------
# PA_deg_sky_band5_robust_minus1 = 148.41 # update this

# minor_angle_rad_sky_band5_robust_minus1, major_angle_rad_sky_band5_robust_minus1, minor_angle_rad_cartesian_band5_robust_minus1, major_angle_rad_cartesian_band5_robust_minus1 = get_minor_major_angles(PA_deg_sky_band5_robust_minus1)
# # ---------------------------------------------


# # Band 4
# # ---------------------------------------------
# PA_deg_sky_band4 = 145.9 # This value is from CASA 

# minor_angle_rad_sky_band4, major_angle_rad_sky_band4, minor_angle_rad_cartesian_band4, major_angle_rad_cartesian_band4 = get_minor_major_angles(PA_deg_sky_band4)
# # ---------------------------------------------
# # Band 4 (nterms = 2)
# # ---------------------------------------------
# PA_deg_sky_band4_nterms2 = 139.31 # This value is from CASA and has been updated for nterms = 2

# minor_angle_rad_sky_band4_nterms2, major_angle_rad_sky_band4_nterms2, minor_angle_rad_cartesian_band4_nterms2, major_angle_rad_cartesian_band4_nterms2 = get_minor_major_angles(PA_deg_sky_band4_nterms2)
# # ---------------------------------------------

# PA_deg_sky_band4_nterms2_smooth = 148.30

# minor_angle_rad_sky_band4_nterms2_smooth, major_angle_rad_sky_band4_nterms2_smooth, minor_angle_rad_cartesian_band4_nterms2_smooth, major_angle_rad_cartesian_band4_nterms2_smooth = get_minor_major_angles(PA_deg_sky_band4_nterms2_smooth)
# # ---------------------------------------------

# # Band 7
# # ---------------------------------------------
# # I believe this is for nterms = 2 
# # PA_deg_sky_band7 = 139.65 # This value is from fitting the new band 7 file (image component size, or 148.86 deconvolved)

# # minor_angle_rad_sky_band7, major_angle_rad_sky_band7, minor_angle_rad_cartesian_band7, major_angle_rad_cartesian_band7 = get_minor_major_angles(PA_deg_sky_band7)
# # ---------------------------------------------


# # Band 7
# # ---------------------------------------------
# PA_deg_sky_band7_nterms2 = 139.65

# minor_angle_rad_sky_band7_nterms2, major_angle_rad_sky_band7_nterms2, minor_angle_rad_cartesian_band7_nterms2, major_angle_rad_cartesian_band7_nterms2 = get_minor_major_angles(PA_deg_sky_band7_nterms2)
# # ---------------------------------------------
# PA_deg_sky_band7_nterms2_smooth = PA_deg_sky_band7_nterms2

# minor_angle_rad_sky_band7_nterms2_smooth, major_angle_rad_sky_band7_nterms2_smooth, minor_angle_rad_cartesian_band7_nterms2_smooth, major_angle_rad_cartesian_band7_nterms2_smooth = get_minor_major_angles(PA_deg_sky_band7_nterms2_smooth)
# # ---------------------------------------------




# Rms errors
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# RMS errors (mJy)
# ------------------------------------------------------------------------------------------

rms_data = {

    # Band 4
    "band4_nterms2": {
        "I": [0.0115, 0.0145, 0.0101, 0.0189],
        "Q": [0.0098, 0.0097, 0.0095, 0.0096],
        "U": [0.0102, 0.0105, 0.0098, 0.0100]
    },

    "band4_nterms2_smooth": {
        "I": [0.0119, 0.0121, 0.0119, 0.0139],
        "Q": [0.0110, 0.0109, 0.0111, 0.0113],
        "U": [0.0104, 0.0107, 0.0105, 0.0100]
    },
    
    "band4_nterms2_smooth_B6": {
        "I": [],
        "Q": [],
        "U": []
    },
    
    "band4_nterms2_smooth_B6_B7": {
        "I": [0.0120, 0.0142, 0.0126, 0.0134],
        "Q": [0.0110, 0.0114, 0.0116, 0.0134],
        "U": [0.0110, 0.0117, 0.0119, 0.0107]
    },

    # Band 5
    "band5_v0": {
        "I": [0.0135, 0.0147, 0.0149, 0.0144],
        "Q": [0.0127, 0.0137, 0.0129, 0.0128],
        "U": [0.0150, 0.0155, 0.0157, 0.0148]
    },

    "band5": {
        "I": [0.0138, 0.0139, 0.0142, 0.0136],
        "Q": [0.0147, 0.0131, 0.0133, 0.0137],
        "U": [0.0152, 0.0159, 0.0141, 0.0156]
    },

    "band5_nterms2": {
        "I": [0.0146, 0.0142, 0.0144, 0.0136],
        "Q": [0.0139, 0.0132, 0.0134, 0.0134],
        "U": [0.0146, 0.0144, 0.0143, 0.0146]
    },

    "band5_robust_minus2": {
        "I": [0.0889, 0.0955, 0.0956, 0.0960],
        "Q": [0.0815, 0.0830, 0.0803, 0.0842],
        "U": [0.0853, 0.0865, 0.0848, 0.0832]
    },

    # Band 6
    "band6_smooth": {
        "I": [0.2943, 0.2446, 0.2891, 0.5541],
        "Q": [0.0683, 0.0686, 0.0651, 0.0393],
        "U": [0.1024, 0.0695, 0.6244, 0.1065]
    },
    
    "band6_smooth_B7": {
        "I": [],
        "Q": [],
        "U": []
    },

    # Band 7
    "band7_nterms2": {
        "I": [0.0340, 0.0382, 0.0375, 0.0310],
        "Q": [0.0256, 0.0266, 0.0268, 0.0260],
        "U": [0.0257, 0.0243, 0.0234, 0.0240]
    }
    
#    "band7_nterms2_smooth": {
#         "I": [],
#         "Q": [],
#         "U": []
#     }
    
#   "band7_nterms2_smooth_B6": {
#         "I": [],
#         "Q": [],
#         "U": []
#     }
}

# ------------------------------------------------------------------------------------------
# Compute means and create variables
# ------------------------------------------------------------------------------------------

for band, stokes in rms_data.items():

    I_mean = np.mean(stokes["I"])
    Q_mean = np.mean(stokes["Q"])
    U_mean = np.mean(stokes["U"])

    POLI_mean = np.mean([Q_mean, U_mean])

    globals()[f"StokesI_err_mJy_{band}"] = I_mean
    globals()[f"StokesQ_err_mJy_{band}"] = Q_mean
    globals()[f"StokesU_err_mJy_{band}"] = U_mean
    globals()[f"POLI_err_mJy_{band}"] = POLI_mean
    
    
# Old way:
# # --------
# # Band 4 (nterms 2) 
# # --------
# StokesI_err_mJy_band4_nterms2 = np.mean([0.0115, 0.0145, 0.0101, 0.0189])
# StokesQ_err_mJy_band4_nterms2 = np.mean([0.0098, 0.0097, 0.0095, 0.0096]) 
# StokesU_err_mJy_band4_nterms2 = np.mean([0.0102, 0.0105, 0.0098, 0.0100]) 
# POLI_err_mJy_band4_nterms2 = np.mean([StokesQ_err_mJy_band4_nterms2, StokesU_err_mJy_band4_nterms2])
# # ------------------------------------------------------------------------------------------
# StokesI_err_mJy_band4_nterms2_smooth = np.mean([0.0119, 0.0121, 0.0119, 0.0139])
# StokesQ_err_mJy_band4_nterms2_smooth = np.mean([0.0110, 0.0109, 0.0111, 0.0113]) 
# StokesU_err_mJy_band4_nterms2_smooth = np.mean([0.0104, 0.0107, 0.0105, 0.0100]) 
# POLI_err_mJy_band4_nterms2_smooth = np.mean([StokesQ_err_mJy_band4_nterms2_smooth, StokesU_err_mJy_band4_nterms2_smooth])
# # ------------------------------------------------------------------------------------------
# # --------
# # Band 5
# # --------
# StokesI_err_mJy_band5_v0 = np.mean([0.0135, 0.0147, 0.0149, 0.0144])
# StokesQ_err_mJy_band5_v0 = np.mean([0.0127, 0.0137, 0.0129, 0.0128]) 
# StokesU_err_mJy_band5_v0 = np.mean([0.0150, 0.0155, 0.0157, 0.0148]) 
# POLI_err_mJy_band5_v0 = np.mean([StokesQ_err_mJy_band5_v0, StokesU_err_mJy_band5_v0])
# # ------------------------------------------------------------------------------------------
# StokesI_err_mJy_band5 = np.mean([0.0138, 0.0139, 0.0142, 0.0136])
# StokesQ_err_mJy_band5 = np.mean([0.0147, 0.0131, 0.0133, 0.0137]) 
# StokesU_err_mJy_band5 = np.mean([0.0152, 0.0159, 0.0141, 0.0156]) 
# POLI_err_mJy_band5 = np.mean([StokesQ_err_mJy_band5, StokesU_err_mJy_band5])
# # ------------------------------------------------------------------------------------------
# StokesI_err_mJy_band5_nterms2 = np.mean([0.0146, 0.0142, 0.0144, 0.0136])
# StokesQ_err_mJy_band5_nterms2 = np.mean([0.0139, 0.0132, 0.0134, 0.0134]) 
# StokesU_err_mJy_band5_nterms2 = np.mean([0.0146, 0.0144, 0.0143, 0.0146]) 
# POLI_err_mJy_band5_nterms2 = np.mean([StokesQ_err_mJy_band5_nterms2, StokesU_err_mJy_band5_nterms2])

# # ------------------------------------------------------------------------------------------
# StokesI_err_mJy_band5_robust_minus2 = np.mean([0.0889, 0.0955, 0.0956, 0.0960])
# StokesQ_err_mJy_band5_robust_minus2 = np.mean([0.0815, 0.0830, 0.0803, 0.0842]) 
# StokesU_err_mJy_band5_robust_minus2 = np.mean([0.0853, 0.0865, 0.0848, 0.0832]) 
# POLI_err_mJy_band5_robust_minus2 = np.mean([StokesQ_err_mJy_band5_robust_minus2, StokesU_err_mJy_band5_robust_minus2])
# # ------------------------------------------------------------------------------------------
# # StokesI_err_mJy_band5_robust_minus1 = np.mean([0.2081, 0.1469, 0.3065, 0.3852])
# # StokesQ_err_mJy_band5_robust_minus1 = np.mean([0.0667, 0.6199, 0.5987, 0.6216]) 
# # StokesU_err_mJy_band5_robust_minus1 = np.mean([0.0931, 0.1064, 0.0693, 0.0659]) 
# # POLI_err_mJy_band5_robust_minus1 = np.mean([StokesQ_err_mJy_band5_robust_minus1, StokesU_err_mJy_band5_robust_minus1])
# # ------------------------------------------------------------------------------------------
# # Band 6 smooth
# # ---------------------------------------------
# StokesI_err_mJy_band6_smooth = np.mean([0.2943, 0.2446, 0.2891, 0.5541])
# StokesQ_err_mJy_band6_smooth = np.mean([0.0683, 0.0686, 0.0651, 0.0393]) 
# StokesU_err_mJy_band6_smooth = np.mean([0.1024, 0.0695, 0.6244, 0.1065]) 
# POLI_err_mJy_band6_smooth = np.mean([StokesQ_err_mJy_band6_smooth, StokesU_err_mJy_band6_smooth])
# # ------------------------------------------------------------------------------------------
# # Band 7
# # ---------------------------------------------
# # Get these values using STDEV on CARTA
# # StokesQ_err_mJy_band7_uncalibrated = 0.03
# # StokesU_err_mJy_band7_uncalibrated = 0.03
# # POLI_err_mJy_band7_uncalibrated = 0.05 
# # ------------------------------------------------------------------------------------------
# # StokesI_err_mJy_band7 = np.mean([nan]) 
# # StokesQ_err_mJy_band7 = np.mean([nan]) 
# # StokesU_err_mJy_band7 = np.mean([nan]) 
# # POLI_err_mJy_band7 = np.mean([StokesQ_err_mJy_band7, StokesU_err_mJy_band7])
# # ------------------------------------------------------------------------------------------
# StokesI_err_mJy_band7_nterms2 = np.mean([0.0340, 0.0382, 0.0375, 0.0310])
# StokesQ_err_mJy_band7_nterms2 = np.mean([0.0256, 0.0266, 0.0268, 0.0260]) 
# StokesU_err_mJy_band7_nterms2 = np.mean([0.0257, 0.0243, 0.0234, 0.0240]) 
# POLI_err_mJy_band7_nterms2 = np.mean([StokesQ_err_mJy_band7_nterms2, StokesU_err_mJy_band7_nterms2])
# # ------------------------------------------------------------------------------------------

# Vector Constants
# ------------------------------------------------------------------------------------------
reference_length_AU = 100
max_length_pix = 400  # Maximum length of the vector in pixels for 100% polarization
reference_fraction = 0.03




# --------------------------------------------------------------------------------------------------------
# Vector step sizes
# --------------------------------------------------------------------------------------------------------

step_values = {
    "band4": 4,
    "band4_nterms2": 4,           # update later if needed
    "band4_nterms2_smooth": 8,
     "band4_nterms2_smooth_B6": 100, # This needs to be updated obviously
     "band4_nterms2_smooth_B6_B7": 100, # This needs to be updated obviously
    "band5": 8,
    "band6": 6,
    "band6_smooth": 8,
    "band6_smooth_B7": 8,  # This needs to be updated obviously
    "band7_nterms2": 6,
    "band7_nterms2_smooth": 100, # This needs to be updated obviously
    "band7_nterms2_smooth_B6": 100, # This needs to be updated obviously
}

for band, val in step_values.items():
    globals()[f"step_{band}"] = val


# --------------------------------------------------------------------------------------------------------
# Vector lengths (pixels)
# --------------------------------------------------------------------------------------------------------

vector_length_pix_const = 0

vector_length_values = {
    "band4": 2,
    "band4_nterms2": 2,           # update later if needed
    "band5": 6,
    "band6": 4,
    "band7_nterms2": 4
}

for band, val in vector_length_values.items():
    globals()[f"vector_len_pix_{band}"] = val
 # --------------------------------------------------------------------------------------------------------   
    
# Old way keep here for now:  
# step_band4 = 4
# step_band4_nterms2 = step_band4 # this needs updating
# step_band4_nterms2_smooth = 8 
# step_band5 = 8 
# step_band6 = 6  
# step_band6_smooth = 8
# # step_band7 = 6
# step_band7_nterms2 = 6

# vector_length_pix_const = 0 
# vector_len_pix_band4 = 2
# vector_len_pix_band4_nterms2 = vector_len_pix_band4 # This needs updating
# vector_len_pix_band5 = 6  
# vector_len_pix_band6 = 4
# # vector_len_pix_band7 = 4
# vector_len_pix_band7_nterms2 = 4
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
