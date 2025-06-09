# Import constants
# ---------------------------------------------------------------------------------------
import sys

# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")

# Now you can import constants.py
import constants
# ---------------------------------------------------------------------------------------


# Load in the function:
# ------------------------------------------
from FITS_Image_Functions import * 
from PolarizationFunctions import *
# ------------------------------------------


def get_plotting_parameters(StokesI_header, StokesI_wcs, band):
    # Dictionary for band constants
    # Get the constants depending on what band of data we are using:
    if band == 6:
        centre_str = constants.centre_str_band6

    elif band == 4:
        centre_str = constants.centre_str_band4  # Changed to band 4   
    elif band == 7:
        centre_str = constants.centre_str_band7 
    else:
        return "Invalid band option"
    
    # Get beam information
    beam_info = get_beam_info(StokesI_header)
    BMAJ_deg, BMIN_deg, BMAJ_pix, BMIN_pix, BPA_deg_sky, BPA_deg_cartesian = beam_info
    BPA_rad_sky = np.radians(BPA_deg_sky)
    BPA_rad_cartesian = np.radians(BPA_deg_cartesian)
    
    # Get the length of the reference pixel
    reference_length_pix = length_in_pixels(constants.reference_length_AU, constants.distance_pc, StokesI_header)
    
    # Get the centre pixel values
    RA_centre_pix, Dec_centre_pix = string_to_pixel(centre_str, StokesI_wcs)
    
    min_str    = constants.min_str_band6
    max_str    = constants.max_str_band6
    RA_min_pix, Dec_min_pix = string_to_pixel(min_str, StokesI_wcs)
    RA_max_pix, Dec_max_pix = string_to_pixel(max_str, StokesI_wcs)

    # Define the plot boundaries (xmin, xmax, ymin, ymax)
    xmin, xmax = RA_max_pix, RA_min_pix
    ymin, ymax = Dec_min_pix, Dec_max_pix



#     if band == 6:
#         RA_min_pix, Dec_min_pix = string_to_pixel(min_str, StokesI_wcs)
#         RA_max_pix, Dec_max_pix = string_to_pixel(max_str, StokesI_wcs)
        
#         # Define the plot boundaries (xmin, xmax, ymin, ymax)
#         xmin, xmax = RA_max_pix, RA_min_pix
#         ymin, ymax = Dec_min_pix, Dec_max_pix
        
#     elif band == 4:
#         x_val = 40
#         y_val = 40
        
#         xmin = RA_centre_pix - x_val
#         xmax = RA_centre_pix + x_val
#         ymin = Dec_centre_pix - y_val
#         ymax = Dec_centre_pix + y_val
#     elif band == 7:
#         x_val = 40
#         y_val = 40
        
#         xmin = RA_centre_pix - x_val
#         xmax = RA_centre_pix + x_val
#         ymin = Dec_centre_pix - y_val
#         ymax = Dec_centre_pix + y_val
        
#     else:
#         return "Invalid band option"

#     x_val = 50
#     y_val = 50

#     xmin = RA_centre_pix - x_val
#     xmax = RA_centre_pix + x_val
#     ymin = Dec_centre_pix - y_val
#     ymax = Dec_centre_pix + y_val
    
    
    
   
    
    return BMAJ_deg, BMIN_deg, BMAJ_pix, BMIN_pix, BPA_deg_cartesian, reference_length_pix, RA_centre_pix, Dec_centre_pix, xmin, xmax, ymin, ymax




def generate_polarization_vectors_band47(ny, nx,
                                         RA_centre_pix, Dec_centre_pix,
                                         uniform_angle,
                                         StokesI_mJy, 
                                         POLI_mJy, POLI_err_mJy,
                                         PA_real_sky_rad, PA_err_deg,
                                         band):
    """
    Generates polarization vectors for different grid configurations and calculates the Stokes U and Q grids.

    Parameters:
    - ny, nx: Dimensions of the grid
    - RA_centre_pix, Dec_centre_pix: Centre coordinates in pixels
    - uniform_angle: The uniform angle for PA grid
    - POLI_mJy, POLI_err_mJy: Polarization values and their uncertainties
    - PA_real_sky_rad, PA_err_deg: Polarization angles (in radians) and their uncertainties

    Returns:
    - A dictionary containing the vector and angle data for different grid configurations, as well as the Stokes Q and U grids.
    """
    
    # Get vector and angles for the actual data
    vector_data_actual_cartesian, vector_angle_actual_sky = make_vectors_band47(ny, nx,  
                                                                                POLI_mJy, POLI_err_mJy,
                                                                                PA_real_sky_rad, PA_err_deg,
                                                                                band)
    
    # Make the PA grids for uniform and Azimuthal
    PA_grid_100Uniform   = make_PA_grid_100Uniform(ny,   nx, uniform_angle)
    PA_grid_100Azimuthal = make_PA_grid_100Azimuthal(ny, nx, RA_centre_pix, Dec_centre_pix)  
    
    # Get the vector and angle data for the 100 Uniform case 
    vector_data_100Uniform_cartesian, vector_angle_100Uniform_sky = make_vectors_band47(ny, nx,  
                                                                                       POLI_mJy, POLI_err_mJy,
                                                                                       PA_grid_100Uniform, PA_err_deg,
                                                                                       band)
    
    # Get the vector and angle data for the 100 Azimuthal case 
    vector_data_100Azimuthal_cartesian, vector_angle_100Azimuthal_sky = make_vectors_band47(ny, nx,  
                                                                                           POLI_mJy, POLI_err_mJy,
                                                                                           PA_grid_100Azimuthal, PA_err_deg,
                                                                                           band)
    
    # Get Stokes Q and U grids
    StokesQ_grid_100Uniform,   StokesU_grid_100Uniform   = recover_StokesQU(PA_grid_100Uniform,   StokesI_mJy, ny, nx)
    StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal = recover_StokesQU(PA_grid_100Azimuthal, StokesI_mJy, ny, nx)
    
    # Organize the results in a dictionary for easy access
    results = {
        'vector_data_actual_cartesian': vector_data_actual_cartesian,
        'vector_angle_actual_sky': vector_angle_actual_sky,
        'vector_data_100Uniform_cartesian': vector_data_100Uniform_cartesian,
        'vector_angle_100Uniform_sky': vector_angle_100Uniform_sky,
        'vector_data_100Azimuthal_cartesian': vector_data_100Azimuthal_cartesian,
        'vector_angle_100Azimuthal_sky': vector_angle_100Azimuthal_sky,
        'StokesQ_grid_100Uniform': StokesQ_grid_100Uniform,
        'StokesU_grid_100Uniform': StokesU_grid_100Uniform,
        'StokesQ_grid_100Azimuthal': StokesQ_grid_100Azimuthal,
        'StokesU_grid_100Azimuthal': StokesU_grid_100Azimuthal
    }
    
    return results






def generate_polarization_vectors_band6(ny, nx,
                                        RA_centre_pix, Dec_centre_pix,
                                        uniform_angle,
                                        StokesI_mJy, StokesI_err_mJy,
                                        POLI_mJy, POLI_err_mJy,
                                        PA_real_sky_rad, PA_err_deg):
    """
    Generates polarization vectors for different grid configurations and calculates the Stokes U and Q grids.

    Parameters:
    - ny, nx: Dimensions of the grid
    - RA_centre_pix, Dec_centre_pix: Centre coordinates in pixels
    - uniform_angle: The uniform angle for PA grid
    - StokesI_mJy, StokesI_err_mJy: Stokes I values and their uncertainties
    - POLI_mJy, POLI_err_mJy: Polarization values and their uncertainties
    - PA_real_sky_rad, PA_err_deg: Polarization angles (in radians) and their uncertainties

    Returns:
    - A dictionary containing the vector and angle data for different grid configurations, as well as the Stokes Q and U grids.
    """
    
    # Get vector and angles for the actual data
    vector_data_actual_cartesian, vector_angle_actual_sky = make_vectors_band6(ny, nx,  
                                                                                StokesI_mJy, StokesI_err_mJy,
                                                                                POLI_mJy, POLI_err_mJy,
                                                                                PA_real_sky_rad, PA_err_deg)
    
    # Make the PA grids for uniform and Azimuthal
    PA_grid_100Uniform = make_PA_grid_100Uniform(ny, nx, uniform_angle)
    PA_grid_100Azimuthal = make_PA_grid_100Azimuthal(ny, nx, RA_centre_pix, Dec_centre_pix)  
    
    # Get the vector and angle data for the 100 Uniform case 
    vector_data_100Uniform_cartesian, vector_angle_100Uniform_sky = make_vectors_band6(ny, nx,  
                                                                                         StokesI_mJy, StokesI_err_mJy,
                                                                                         POLI_mJy, POLI_err_mJy,
                                                                                         PA_grid_100Uniform, PA_err_deg)
    
    # Get the vector and angle data for the 100 Azimuthal case 
    vector_data_100Azimuthal_cartesian, vector_angle_100Azimuthal_sky = make_vectors_band6(ny, nx,  
                                                                                             StokesI_mJy, StokesI_err_mJy,
                                                                                             POLI_mJy, POLI_err_mJy,
                                                                                             PA_grid_100Azimuthal, PA_err_deg)
    
    # Get Stokes Q and U grids
    StokesQ_grid_100Uniform, StokesU_grid_100Uniform = recover_StokesQU(PA_grid_100Uniform, StokesI_mJy, ny, nx)
    StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal = recover_StokesQU(PA_grid_100Azimuthal, StokesI_mJy, ny, nx)
    
    # Organize the results in a dictionary for easy access
    results = {
        'vector_data_actual_cartesian': vector_data_actual_cartesian,
        'vector_angle_actual_sky': vector_angle_actual_sky,
        'vector_data_100Uniform_cartesian': vector_data_100Uniform_cartesian,
        'vector_angle_100Uniform_sky': vector_angle_100Uniform_sky,
        'vector_data_100Azimuthal_cartesian': vector_data_100Azimuthal_cartesian,
        'vector_angle_100Azimuthal_sky': vector_angle_100Azimuthal_sky,
        'StokesQ_grid_100Uniform': StokesQ_grid_100Uniform,
        'StokesU_grid_100Uniform': StokesU_grid_100Uniform,
        'StokesQ_grid_100Azimuthal': StokesQ_grid_100Azimuthal,
        'StokesU_grid_100Azimuthal': StokesU_grid_100Azimuthal
    }
    
    return results
 
    
    
    
    