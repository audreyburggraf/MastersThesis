import sys

# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")  # Replace with the actual path

# Now you can import constants.py
import constants

# Use the variable from constants.py
band6_data_folder_path = constants.band6_data_folder_path
functions_folder_path = constants.functions_folder_path




# Load in the function:
# ------------------------------------------
from C29_functions import * 
from FITS_Image_Functions import * 
from PolarizationFunctions import *
from PlottingWithFunction import * 
from IntroductionFunctions import *
# ------------------------------------------


# Define file paths
StokesI_file      = band6_data_folder_path + "c2d_989_StokesI_233GHz.fits"
StokesI_err_file  = band6_data_folder_path + "c2d_989_StokesIerr_233GHz.fits"
StokesQ_file      = band6_data_folder_path + "c2d_989_StokesQ_233GHz.fits"
StokesU_file      = band6_data_folder_path + "c2d_989_StokesU_233GHz.fits"
PA_file           = band6_data_folder_path + "c2d_989_POLA_233GHz.fits"
PA_err_file       = band6_data_folder_path + "c2d_989_POLAerr_233GHz.fits"
POLF_file         = band6_data_folder_path + "c2d_989_POLF_233GHz.fits"
POLI_file         = band6_data_folder_path + "c2d_989_POLI_233GHz.fits"
POLI_err_file     = band6_data_folder_path + "c2d_989_POLIerr_233GHz.fits"


# Stokes I
# -------------------------------------------------------------------------------------------------------
StokesI_header, _, StokesI_Jy, StokesI_wcs = read_in_file(StokesI_file)
StokesI_mJy = convert_jy_to_mjy(StokesI_Jy)

# Define min/max for plotting
StokesI_custom_vmin_band6 = -0.2
StokesI_custom_vmax_band6 = 89.85


StokesI_stretched_mJy, StokesI_unstretched_cbar_ticks = normalize_stokesI_for_cmap(StokesI_mJy,
                                                                                  custom_min = StokesI_custom_vmin_band6, 
                                                                                  custom_max = StokesI_custom_vmax_band6)
# -------------------------------------------------------------------------------------------------------


# Stokes I error
# -------------------------------------------------------------------------------------------------------
_, _, StokesI_err_Jy, _ = read_in_file(StokesI_err_file, dimensions=2)
StokesI_err_mJy = convert_jy_to_mjy(StokesI_err_Jy)
# -------------------------------------------------------------------------------------------------------


# Stokes Q
# -------------------------------------------------------------------------------------------------------
_, _, StokesQ_Jy, _ = read_in_file(StokesQ_file)
StokesQ_mJy = convert_jy_to_mjy(StokesQ_Jy)
# -------------------------------------------------------------------------------------------------------


# Stokes U
# -------------------------------------------------------------------------------------------------------
_, _, StokesU_Jy, _ = read_in_file(StokesU_file)
StokesU_mJy = convert_jy_to_mjy(StokesU_Jy)
# -------------------------------------------------------------------------------------------------------


# Polarization Intensity
# -------------------------------------------------------------------------------------------------------
_, _, POLI_Jy, _ = read_in_file(POLI_file)
POLI_mJy = convert_jy_to_mjy(POLI_Jy)
# -------------------------------------------------------------------------------------------------------


# Polarization Intensity error
# -------------------------------------------------------------------------------------------------------
_, _, POLI_err_Jy, _ = read_in_file(POLI_err_file)
POLI_err_mJy = convert_jy_to_mjy(POLI_err_Jy)
# -------------------------------------------------------------------------------------------------------



# Polarization Angle
# -------------------------------------------------------------------------------------------------------
_, _, PA_deg, _ = read_in_file(PA_file)
PA_rad = np.radians(PA_deg)
# -------------------------------------------------------------------------------------------------------



# Polarization Angle error
# -------------------------------------------------------------------------------------------------------
_, _, PA_err_deg, _ = read_in_file(PA_err_file)
# -------------------------------------------------------------------------------------------------------


# Polarization Fraction
# -------------------------------------------------------------------------------------------------------
_, _, POLF, _ = read_in_file(POLF_file)
# -------------------------------------------------------------------------------------------------------


# Polarization Fraction error
# -------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------



PA_calc_rad = calculate_polarization_angle(StokesQ_Jy, StokesU_Jy)
POLF_calc = calculate_polarized_fraction(StokesQ_mJy, StokesU_mJy, StokesI_mJy)
POLI_calc = calculate_polarized_intensity(StokesQ_mJy, StokesU_mJy)

nx, ny = StokesI_mJy.shape





BMAJ_pix, BMIN_pix, BPA_deg_cartesian, reference_length_pix, RA_centre_pix, Dec_centre_pix, xmin, xmax, ymin, ymax = get_plotting_parameters(StokesI_header, StokesI_wcs, 6)




# Find the vectors
# -------------------------------------------------------------------------------------------------------
results = generate_polarization_vectors_band6(ny, nx,
                                              RA_centre_pix, Dec_centre_pix,
                                              constants.minor_angle_rad_sky_band6,
                                              StokesI_mJy, StokesI_err_mJy, 
                                              POLI_mJy, POLI_err_mJy,
                                              PA_rad, PA_err_deg)
# -------------------------------------------------------------------------------------------------------

# Accessing the actual vector data and anglesf
vector_data_actual_cartesian = results['vector_data_actual_cartesian']
vector_angle_actual_sky = results['vector_angle_actual_sky']


vector_data_100Uniform_cartesian = results['vector_data_100Uniform_cartesian']
vector_angle_100Uniform_sky = results['vector_angle_100Uniform_sky']


vector_data_100Azimuthal_cartesian = results['vector_data_100Azimuthal_cartesian']
vector_angle_100Azimuthal_sky      = results['vector_angle_100Azimuthal_sky']

StokesQ_grid_100Uniform  = results['StokesQ_grid_100Uniform']
StokesU_grid_100Uniform  = results['StokesU_grid_100Uniform']

StokesQ_grid_100Azimuthal = results['StokesQ_grid_100Azimuthal']
StokesU_grid_100Azimuthal = results['StokesU_grid_100Azimuthal']
# -------------------------------------------------------------------------------------------------------



