import sys
from astropy.io import fits

# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")  # Replace with the actual path

# Now you can import constants.py
import constants

# Use the variable from constants.py
band5_data_folder_path = constants.band5_data_folder_path
functions_folder_path = constants.functions_folder_path




# Load in the function:
# ------------------------------------------
from C29_functions import * 
from FITS_Image_Functions import * 
from PolarizationFunctions import *
from PlottingWithFunction import * 
from IntroductionFunctions import *
from POLF_Functions import *
# ------------------------------------------


# Define file paths
loc = band5_data_folder_path 
StokesI_file         = loc + "IRS63_BAND5_StokesI_redo_nterms1.fits"
StokesQ_file         = loc + "IRS63_BAND5_StokesQ_redo_nterms1.fits"
StokesU_file         = loc + "IRS63_BAND5_StokesU_redo_nterms1.fits"
POLI_biased_file     = loc + "POLI_biased_mJy_BAND5.fits"
POLI_err_file        = loc + "POLI_err_mJy_BAND5.fits"
POLI_debiased_file = loc + "POLI_debiased_mJy_BAND5.fits"


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
StokesI_err_mJy = np.full((ny, nx), constants.StokesI_err_mJy_band5)
# -------------------------------------------------------------------------------------------------------


# Stokes Q
# -------------------------------------------------------------------------------------------------------
_, _, StokesQ_Jy, _ = read_in_file(StokesQ_file)
StokesQ_mJy = convert_jy_to_mjy(StokesQ_Jy)
# -------------------------------------------------------------------------------------------------------
# Stokes Q Error
# -------------------------------------------------------------------------------------------------------
StokesQ_err_mJy = np.full((ny, nx), constants.StokesQ_err_mJy_band5)
# -------------------------------------------------------------------------------------------------------


# Stokes U
# -------------------------------------------------------------------------------------------------------
_, _, StokesU_Jy, _ = read_in_file(StokesU_file)
StokesU_mJy = convert_jy_to_mjy(StokesU_Jy)
# -------------------------------------------------------------------------------------------------------
# Stokes U Error
# -------------------------------------------------------------------------------------------------------
StokesU_err_mJy = np.full((ny, nx), constants.StokesU_err_mJy_band5)
# -------------------------------------------------------------------------------------------------------



# Polarization Intensity
# -------------------------------------------------------------------------------------------------------
# POLI_calc = calculate_polarized_intensity(StokesQ_mJy, StokesU_mJy)
# POLI_mJy = POLI_calc

# # Ensure POLI_mJy is shaped (ny, nx)
# POLI_mJy_resized = POLI_mJy[:ny, :nx]


# # Save data
# hdu = fits.PrimaryHDU(data=POLI_mJy, header=StokesI_header)
# hdu.writeto(band5_data_folder_path + "POLI_mJy_calculated_BAND5.fits", overwrite=True)

# # _, _, POLI_Jy, _ = read_in_file(POLI_file)
# # POLI_mJy = convert_jy_to_mjy(POLI_Jy)
# # -------------------------------------------------------------------------------------------------------
# # Polarization Intensity error
# # -------------------------------------------------------------------------------------------------------
# POLI_err_mJy = calculate_polarized_intensity_err(StokesQ_mJy, StokesU_mJy, StokesQ_err_mJy, StokesU_err_mJy)
# -------------------------------------------------------------------------------------------------------
# Polarization Intensity
# -------------------------------------------------------------------------------------------------------
_, _, POLI_biased_mJy, _ = read_in_file(POLI_biased_file, dimensions = 2)

_, _, POLI_debiased_mJy, _ = read_in_file(POLI_debiased_file, dimensions = 2)

_, _, POLI_err_mJy, _ = read_in_file(POLI_err_file, dimensions = 2)
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

# find_POLF_avg("Band 5 v0", POLF, StokesI_mJy, band5_v0_data_folder_path)
# -------------------------------------------------------------------------------------------------------
# Polarized Fraction Error
# -------------------------------------------------------------------------------------------------------
POLF_err = calculate_polarized_fraction_err(StokesQ_mJy, StokesU_mJy, StokesI_mJy, 
                                            StokesQ_err_mJy, StokesU_err_mJy, StokesI_err_mJy)
# -------------------------------------------------------------------------------------------------------



BMAJ_deg, BMIN_deg, BMAJ_pix, BMIN_pix, BPA_deg_cartesian, reference_length_pix, RA_centre_pix, Dec_centre_pix, xmin, xmax, ymin, ymax = get_plotting_parameters(StokesI_header, StokesI_wcs, 'Band 5')

delta = 100
xmin = xmin + delta
xmax = xmax - delta
ymin = ymin + delta
ymax = ymax - delta


# Remove stream
# -------------------------------------------------------------------------------------------------------
stream_xmin = 540
stream_xmax = 580
stream_ymin = 490
stream_ymax = 540


POLI_mJy_nostream = POLI_debiased_mJy.copy()
POLI_mJy_nostream[stream_ymin:stream_ymax+1, stream_xmin:stream_xmax+1] = np.nan
# -------------------------------------------------------------------------------------------------------


# Find the vectors
# -------------------------------------------------------------------------------------------------------
# in the file:MakingGridFunctions.py make_vectors_band47
results = generate_polarization_vectors(ny, nx,
                                        xmin, xmax, ymin, ymax, # This is for the nterms test
                                        RA_centre_pix, Dec_centre_pix,
                                        constants.minor_angle_rad_sky_band5_v0,
                                        StokesI_mJy, StokesI_err_mJy, 
                                        POLI_mJy_nostream, POLI_err_mJy,
                                        PA_rad, PA_err_deg,
                                        'Band 5')
# -------------------------------------------------------------------------------------------------------




# Accessing the actual vector data and anglesf
vector_data_actual_cartesian = results['vector_data_actual_cartesian']
vector_angle_actual_sky = results['vector_angle_actual_sky']
vector_angle_actual_sky_errors = results['vector_angle_actual_sky_errors']

vector_data_100Uniform_cartesian = results['vector_data_100Uniform_cartesian']
vector_angle_100Uniform_sky = results['vector_angle_100Uniform_sky']

vector_data_100Azimuthal_cartesian = results['vector_data_100Azimuthal_cartesian']
vector_angle_100Azimuthal_sky      = results['vector_angle_100Azimuthal_sky']

StokesQ_grid_100Uniform  = results['StokesQ_grid_100Uniform']
StokesU_grid_100Uniform  = results['StokesU_grid_100Uniform']

StokesQ_grid_100Azimuthal = results['StokesQ_grid_100Azimuthal']
StokesU_grid_100Azimuthal = results['StokesU_grid_100Azimuthal']
# -------------------------------------------------------------------------------------------------------


# Get stream vectors
# -------------------------------------------------------------------------------------------------------

POLI_mJy_stream = POLI_debiased_mJy.copy()
POLI_mJy_stream[:, :] = np.nan  # Set everything to NaN

# Fill in only the stream region
POLI_mJy_stream[stream_ymin:stream_ymax+1, stream_xmin:stream_xmax+1] = POLI_debiased_mJy[stream_ymin:stream_ymax+1, stream_xmin:stream_xmax+1]

results_stream = generate_polarization_vectors(ny, nx,
                                               xmin, xmax, ymin, ymax, # This is for the nterms test
                                               RA_centre_pix, Dec_centre_pix,
                                               constants.minor_angle_rad_sky_band5,
                                               StokesI_mJy, StokesI_err_mJy,
                                               POLI_mJy_stream, POLI_err_mJy,
                                               #PA_rad, PA_err_deg,
                                               PA_rad + np.pi/2, PA_err_deg,
                                               'Band 5')

stream_vectors = results_stream['vector_data_actual_cartesian']
# -------------------------------------------------------------------------------------------------------



    
    
# separate the data from the stream
# -------------------------------------------------------------------------------------------------------
# stream_vectors = []
# vector_data_actual_cartesian = []


# for row in vector_data_actual_cartesian_wstream:
#     x0, x1, y0, y1 = row

#     if (stream_xmin <= x0 <= stream_xmax and stream_xmin <= x1 <= stream_xmax and
#         stream_ymin <= y0 <= stream_ymax and stream_ymin <= y1 <= stream_ymax):
#         stream_vectors.append(row)
#     else:
#         vector_data_actual_cartesian.append(row)
# -------------------------------------------------------------------------------------------------------




# Save the beam info
# -----------------------------------------------------------------
save_beam_info(StokesI_header, 'Band 5', print_statements=False)
# -----------------------------------------------------------------