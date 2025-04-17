import numpy as np 
import seaborn as sns


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

functions_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FUNCTIONS/"


# Band 6 constants
min_str_band6 = 'J163135.570-240131.200'
max_str_band6 = 'J163135.750-240128.700'
centre_str_band6 = 'J163135.657-240129.935'


# Band 6 PA angle constants
PA_deg_sky_band6 = 148

# The positional angle is the major axis 
major_angle_deg_sky_band6 = PA_deg_sky_band6

# The minor axis is 90 degrees from the major
minor_angle_deg_sky_band6 = major_angle_deg_sky_band6 - 90

# The minor axis is 90 degrees from the major
major_angle_rad_sky_band6 = np.radians(major_angle_deg_sky_band6)
minor_angle_rad_sky_band6 = np.radians(minor_angle_deg_sky_band6)

# We have to offset the minor and major angle by 90 degrees if we want to plot it:
major_angle_deg_cartesian_band6 = major_angle_deg_sky_band6 - 90
minor_angle_deg_cartesian_band6 = minor_angle_deg_sky_band6 - 90

# We can turn the plotted angles into radians
minor_angle_rad_cartesian_band6 = np.radians(minor_angle_deg_cartesian_band6)
major_angle_rad_cartesian_band6 = np.radians(major_angle_deg_cartesian_band6)




# Refefence Lengths and such (for band 6 but same fort band 4?)

# Vector Constants
# ------------------------------------------------------------------------------------------
reference_length_AU = 100
max_length_pix = 400  # Maximum length of the vector in pixels for 100% polarization
reference_fraction = 0.03
step = 7 
vector_length_pix_const = 4
# ------------------------------------------------------------------------------------------
step_band4 = 7
vector_length_pix_const_band4 = 4
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



# Adding titles to each subplot
ratio_grid_plot_titles = [
    "100% Uniform 0% Azimuthal", "0% Uniform 100% Azimuthal", 
    "90% Uniform 10% Azimuthal", "10% Uniform 90% Azimuthal",
    "80% Uniform 20% Azimuthal", "20% Uniform 80% Azimuthal",
    "70% Uniform 30% Azimuthal", "30% Uniform 70% Azimuthal",
    "60% Uniform 40% Azimuthal", "40% Uniform 60% Azimuthal",
    "50% Uniform 50% Azimuthal", "50% Uniform 50% Azimuthal"
]
# ------------------------------------------------------------------------------------------