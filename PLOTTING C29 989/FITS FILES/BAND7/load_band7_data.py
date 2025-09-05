import sys
from astropy.io import fits

# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")  # Replace with the actual path

# Now you can import constants.py
import constants

# Use the variable from constants.py
band7_data_folder_path = constants.band7_data_folder_path
functions_folder_path = constants.functions_folder_path




# Load in the function:
# ------------------------------------------
from C29_functions import * 
from FITS_Image_Functions import * 
from PolarizationFunctions import *
from PlottingWithFunction import * 
from IntroductionFunctions import *
from DebiasFunctions import *
# ------------------------------------------


# Define file paths
# Define file paths
StokesI_file = band7_data_folder_path + "IRS63_BAND7_StokesI_clean_nopbcorr.fits"
StokesQ_file = band7_data_folder_path + "IRS63_BAND7_StokesQ_clean_nopbcorr.fits"
StokesU_file = band7_data_folder_path + "IRS63_BAND7_StokesU_clean_nopbcorr.fits"
POLI_file        = band7_data_folder_path + "POLI_mJy_calculated_BAND7.fits"
POLI_debias_file = band7_data_folder_path + "POLI_mJy_debias_BAND7.fits"


# Stokes I
# -------------------------------------------------------------------------------------------------------
StokesI_header, _, StokesI_Jy, StokesI_wcs = read_in_file(StokesI_file)
StokesI_mJy = convert_jy_to_mjy(StokesI_Jy)

# Stretch the Stokes I data and get the cbar ticks
StokesI_stretched_mJy, StokesI_unstretched_cbar_ticks = normalize_stokesI_for_cmap(StokesI_mJy)

nx, ny = StokesI_mJy.shape
# -------------------------------------------------------------------------------------------------------
# Stokes I error
# -------------------------------------------------------------------------------------------------------
StokesI_err_mJy = np.full((ny, nx), constants.StokesI_err_mJy_band7)
# -------------------------------------------------------------------------------------------------------


# Stokes Q
# -------------------------------------------------------------------------------------------------------
_, _, StokesQ_Jy, _ = read_in_file(StokesQ_file)
StokesQ_mJy = convert_jy_to_mjy(StokesQ_Jy)
# -------------------------------------------------------------------------------------------------------
# Stokes Q Error
# -------------------------------------------------------------------------------------------------------
StokesQ_err_mJy = np.full((ny, nx), constants.StokesQ_err_mJy_band7)
# -------------------------------------------------------------------------------------------------------


# Stokes U
# -------------------------------------------------------------------------------------------------------
_, _, StokesU_Jy, _ = read_in_file(StokesU_file)
StokesU_mJy = convert_jy_to_mjy(StokesU_Jy)
# -------------------------------------------------------------------------------------------------------
# Stokes U Error
# -------------------------------------------------------------------------------------------------------
StokesU_err_mJy = np.full((ny, nx), constants.StokesQ_err_mJy_band7)
# -------------------------------------------------------------------------------------------------------





# Polarization Intensity
# -------------------------------------------------------------------------------------------------------
# POLI_mJy = calculate_polarized_intensity(StokesQ_mJy, StokesU_mJy)



_, _, POLI_mJy, _ = read_in_file(POLI_file, dimensions = 2)
nx, ny = POLI_mJy.shape
_, _, POLI_debias_mJy, _ = read_in_file(POLI_debias_file, dimensions = 2)


# -------------------------------------------------------------------------------------------------------
# Polarization Intensity error
# -------------------------------------------------------------------------------------------------------
POLI_err_mJy = calculate_polarized_intensity_err(StokesQ_mJy, StokesU_mJy, StokesQ_err_mJy, StokesU_err_mJy)
# -------------------------------------------------------------------------------------------------------

# Polarization Angle
# -------------------------------------------------------------------------------------------------------
PA_calc_rad = calculate_polarization_angle(StokesQ_Jy, StokesU_Jy)
PA_rad = PA_calc_rad 
# -------------------------------------------------------------------------------------------------------
# Polarization Angle error
# -------------------------------------------------------------------------------------------------------
PA_err_rad = calculate_polarization_angle_error(StokesQ_mJy, StokesU_mJy, 
                                                StokesQ_err_mJy, StokesU_err_mJy)

PA_err_deg = np.degrees(PA_err_rad)
# -------------------------------------------------------------------------------------------------------


# Polarized Fraction
# -------------------------------------------------------------------------------------------------------
POLF = calculate_polarized_fraction(StokesQ_mJy, StokesU_mJy, StokesI_mJy)
# -------------------------------------------------------------------------------------------------------
# Polarized Fraction Error
# -------------------------------------------------------------------------------------------------------
POLF_err = calculate_polarized_fraction_err(StokesQ_mJy, StokesU_mJy, StokesI_mJy, 
                                            StokesQ_err_mJy, StokesU_err_mJy, StokesI_err_mJy)
# -------------------------------------------------------------------------------------------------------



BMAJ_deg, BMIN_deg, BMAJ_pix, BMIN_pix, BPA_deg_cartesian, reference_length_pix, RA_centre_pix, Dec_centre_pix, xmin, xmax, ymin, ymax = get_plotting_parameters(StokesI_header, StokesI_wcs, 7)




# Find the vectors
# -------------------------------------------------------------------------------------------------------
# in the file:MakingGridFunctions.py make_vectors_band47
results = generate_polarization_vectors(ny, nx,
                                        RA_centre_pix, Dec_centre_pix,
                                        constants.minor_angle_rad_sky_band7,
                                        StokesI_mJy, 
                                        POLI_debias_mJy, POLI_err_mJy,
                                        PA_rad, PA_err_deg,
                                        'Band 7')
# -------------------------------------------------------------------------------------------------------




# Accessing the actual vector data and angles
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




