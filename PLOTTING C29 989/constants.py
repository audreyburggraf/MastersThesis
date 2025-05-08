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

normalized_cbar_ticks = np.array([0, 0.2, 0.4, 0.6, 0.8, 1])



# Constants for the disk
distance_pc = 132


# File locations
band4_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND4/"
band6_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FITS FILES/BAND6/"

band4_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND4/"
band6_carta_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/CARTA FILES/BAND6/"

functions_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FUNCTIONS/"


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
 






# Vector Constants
# ------------------------------------------------------------------------------------------
reference_length_AU = 100
max_length_pix = 400  # Maximum length of the vector in pixels for 100% polarization
reference_fraction = 0.03
step_band6 = 6
step_band4 = 4

vector_length_pix_const = 0 
vector_len_pix_band6 = 4
vector_len_pix_band4 = 2
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