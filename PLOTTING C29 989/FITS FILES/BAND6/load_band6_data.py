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
from FITS_Image_Functions import * 
# ------------------------------------------


# Define file paths
StokesI_file                = band6_data_folder_path + "c2d_989_StokesI_233GHz.fits"
StokesIerr_file             = band6_data_folder_path + "c2d_989_StokesIerr_233GHz.fits"
StokesQ_file                = band6_data_folder_path + "c2d_989_StokesQ_233GHz.fits"
StokesU_file                = band6_data_folder_path + "c2d_989_StokesU_233GHz.fits"
PolarizationAngle_file      = band6_data_folder_path + "c2d_989_POLA_233GHz.fits"
PolarizationAngle_err_file  = band6_data_folder_path + "c2d_989_POLAerr_233GHz.fits"
PolarizationFraction_file   = band6_data_folder_path + "c2d_989_POLF_233GHz.fits"
PolarizedIntensity_file     = band6_data_folder_path + "c2d_989_POLI_233GHz.fits"
PolarizedIntensity_err_file = band6_data_folder_path + "c2d_989_POLIerr_233GHz.fits"

# Load Stokes I data
StokesI_header, StokesI_data_4d_Jy, StokesI_data_2d_Jy, StokesI_wcs = read_in_file(StokesI_file)
StokesI_data_2d_mJy = StokesI_data_2d_Jy * 1000

# Define min/max for plotting
StokesI_custom_vmin_band6 = -0.2
StokesI_custom_vmax_band6 = 89.85

# # Clip and stretch data
# StokesI_clipped = np.clip(StokesI_data_2d_mJy, StokesI_vmin, StokesI_vmax)
# StokesI_stretched = stretch(StokesI_clipped, base=100, vmin=StokesI_vmin, vmax=StokesI_vmax)

# # Color bar ticks
# normalized_cbar_ticks = np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
# StokesI_unstretched_cbar_ticks = unstretch(normalized_cbar_ticks, StokesI_vmin, StokesI_vmax)

# Load other FITS files
_, _, StokesIerr_data_2d_Jy, StokesIerr_wcs = read_in_file(StokesIerr_file, dimensions=2)
StokesIerr_data_2d_mJy = StokesIerr_data_2d_Jy * 1000

StokesQ_header, _, StokesQ_data_2d_Jy, StokesQ_wcs = read_in_file(StokesQ_file)
StokesQ_data_2d_mJy = StokesQ_data_2d_Jy * 1000

StokesU_header, _, StokesU_data_2d_Jy, StokesU_wcs = read_in_file(StokesU_file)
StokesU_data_2d_mJy = StokesU_data_2d_Jy * 1000

PolarizationAngle_header, _, PolarizationAngle_data_2d_deg, PolarizationAngle_wcs = read_in_file(PolarizationAngle_file)
PolarizationAngle_data_2d_rad = PolarizationAngle_data_2d_deg * (np.pi / 180)

PolarizationAngle_err_header, _, PolarizationAngle_err_data_2d_deg, PolarizationAngle_err_wcs = read_in_file(PolarizationAngle_err_file)

PolarizationFraction_header, _, PolarizationFraction_data_2d, PolarizationFraction_wcs = read_in_file(PolarizationFraction_file)

PolarizedIntensity_header, _, PolarizedIntensity_data_2d_Jy, PolarizedIntensity_wcs = read_in_file(PolarizedIntensity_file)
PolarizedIntensity_data_2d_mJy = PolarizedIntensity_data_2d_Jy * 1000

PolarizedIntensity_err_header, _, PolarizedIntensity_err_data_2d_Jy, PolarizedIntensity_err_wcs = read_in_file(PolarizedIntensity_err_file)
PolarizedIntensity_err_data_2d_mJy = PolarizedIntensity_err_data_2d_Jy * 1000