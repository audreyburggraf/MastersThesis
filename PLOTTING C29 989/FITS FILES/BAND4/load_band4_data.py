import sys

# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")  # Replace with the actual path

# Now you can import constants.py
import constants

# Use the variable from constants.py
band4_data_folder_path = constants.band4_data_folder_path
functions_folder_path = constants.functions_folder_path




# Load in the function:
# ------------------------------------------
from FITS_Image_Functions import * 
from PolarizationFunctions import *
# ------------------------------------------


# Stokes I
StokesI_file = band4_data_folder_path + "IRS63_StokesI_clean_selfcal_pbcor_J2000.fits"

StokesI_header, StokesI_data_4d_Jy, StokesI_data_2d_Jy, StokesI_wcs = read_in_file(StokesI_file)

StokesI_data_2d_mJy = StokesI_data_2d_Jy * 1000




# Polarization Intensity Error


# Define the path to the Stokes I FITS file
StokesIerr_file = band4_data_folder_path + "IRS63_selfcal_POLIrms_debiased_J2000.fits"

StokesIerr_header, _ , StokesIerr_data_2d_Jy, StokesIerr_wcs = read_in_file(StokesIerr_file, dimensions = 2)

StokesIerr_data_2d_mJy = StokesIerr_data_2d_Jy * 1000




# Polarization Intensity

# Define the path to the Stokes I FITS file
PolarizedIntensity_file = band4_data_folder_path + "IRS63_selfcal_POLI_debiased_J2000.fits"

PolarizedIntensity_header, PolarizedIntensity_data_4d, PolarizedIntensity_data_2d_Jy, PolarizedIntensity_wcs = read_in_file(PolarizedIntensity_file, dimensions = 2)

PolarizedIntensity_data_2d_mJy = PolarizedIntensity_data_2d_Jy * 1000





# Polarized Intensity Error

# Define the path to the Stokes I FITS file
PolarizedIntensity_err_file = band4_data_folder_path + "IRS63_selfcal_POLI_debiased_J2000.fits"

PolarizedIntensity_err_header, PolarizedIntensity_err_data_4d, PolarizedIntensity_err_data_2d_Jy, PolarizedIntensity_err_wcs = read_in_file(PolarizedIntensity_err_file, dimensions = 2)

PolarizedIntensity_err_data_2d_mJy = PolarizedIntensity_err_data_2d_Jy * 1000





# Polarization Angle

# Define the path to the Stokes I FITS file
PolarizationAngle_file = band4_data_folder_path + "IRS63_selfcal_POLA_debiased_J2000.fits"

PolarizationAngle_header, PolarizationAngle_data_4d_deg, PolarizationAngle_data_2d_deg, PolarizationAngle_wcs = read_in_file(PolarizationAngle_file, dimensions = 2)

PolarizationAngle_data_2d_rad = PolarizationAngle_data_2d_deg * (np.pi / 180)





# Polarization Angle Error

# Define the path to the Stokes I FITS file
PolarizationAngle_err_file = band4_data_folder_path + "IRS63_selfcal_POLArms_debiased_J2000.fits"

PolarizationAngle_err_header, PolarizationAngle_err_data_4d_deg, PolarizationAngle_err_data_2d_deg, PolarizationAngle_err_wcs = read_in_file(PolarizationAngle_err_file, dimensions = 2)

# PolarizationAngle_err_data_2d_rad = PolarizationAngle_err_data_2d_deg * (np.pi / 180)