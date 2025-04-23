import sys

# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")  # Replace with the actual path

# Now you can import constants.py
import constants

# Use the variable from constants.py
band4_data_folder_path = constants.band4_data_folder_path
functions_folder_path  = constants.functions_folder_path




# Load in the function:
# ------------------------------------------
from C29_functions import * 
from FITS_Image_Functions import * 
from PolarizationFunctions import *
from PlottingWithFunction import * 
from IntroductionFunctions import *
# ------------------------------------------

StokesI_file  = band4_data_folder_path + "IRS63_StokesI_clean_selfcal_pbcor_J2000.fits"
POLI_file     = band4_data_folder_path + "IRS63_selfcal_POLI_debiased_J2000.fits"
POLI_err_file = band4_data_folder_path + "IRS63_selfcal_POLIrms_debiased_J2000.fits"
PA_file       = band4_data_folder_path + "IRS63_selfcal_POLA_debiased_J2000.fits"
PA_err_file   = band4_data_folder_path + "IRS63_selfcal_POLArms_debiased_J2000.fits"



# Stokes I
# -------------------------------------------------------------------------------------------------------
StokesI_header, _, StokesI_Jy, StokesI_wcs = read_in_file(StokesI_file)

StokesI_mJy = convert_jy_to_mjy(StokesI_Jy)

# Stretch the Stokes I data and get the cbar ticks
StokesI_stretched_mJy, StokesI_unstretched_cbar_ticks = normalize_stokesI_for_cmap(StokesI_mJy)
# -------------------------------------------------------------------------------------------------------


# Polarization Intensity
# -------------------------------------------------------------------------------------------------------
POLI_header, _, POLI_Jy, _ = read_in_file(POLI_file, dimensions = 2)

POLI_mJy = convert_jy_to_mjy(POLI_Jy)
# -------------------------------------------------------------------------------------------------------


# Polarized Intensity Error
# -------------------------------------------------------------------------------------------------------
POLI_err_header, _, POLI_err_Jy, _ = read_in_file(POLI_err_file, dimensions = 2)

POLI_err_mJy = convert_jy_to_mjy(POLI_err_Jy)
# -------------------------------------------------------------------------------------------------------


# Polarization Angle
# -------------------------------------------------------------------------------------------------------
PA_header, _, PA_deg, _ = read_in_file(PA_file, dimensions = 2)

PA_rad = np.radians(PA_deg)
# -------------------------------------------------------------------------------------------------------


# Polarization Angle Error
# -------------------------------------------------------------------------------------------------------
PA_err_header, _, PA_err_deg, _ = read_in_file(PA_err_file, dimensions = 2)
# -------------------------------------------------------------------------------------------------------



BMAJ_pix, BMIN_pix, BPA_deg_cartesian, reference_length_pix, RA_centre_pix, Dec_centre_pix, xmin, xmax, ymin, ymax = get_plotting_parameters(StokesI_header, StokesI_wcs, 4)


nx, ny = StokesI_mJy.shape


angles_results = get_minor_major_angles(constants.PA_deg_sky_band4)

minor_angle_rad_sky_band4       = angles_results["minor_rad_sky"]
major_angle_rad_sky_band4       = angles_results["major_rad_sky"]
minor_angle_rad_cartesian_band4 = angles_results["minor_rad_cartesian"]
major_angle_rad_cartesian_band4 = angles_results["major_rad_cartesian"]

