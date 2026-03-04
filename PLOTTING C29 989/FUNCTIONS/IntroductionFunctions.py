# Import constants
# ---------------------------------------------------------------------------------------
import sys
import pandas as pd
from pathlib import Path

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
    if band == 'Band 6':
        centre_str = constants.centre_str_band6

    elif band == 'Band 4':
        centre_str = constants.centre_str_band4  # Changed to band 4   
        
    elif band == 'Band 4 nterms2':
        centre_str = constants.centre_str_band4_nterms2 
    
    elif band == 'Band 5':
        centre_str = constants.centre_str_band5 
        
        
    elif band == 'Band 7 nterms2':
        centre_str = constants.centre_str_band7_nterms2
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


# --------------------------------------------------------------------------------------------    
def reccomended_step_count(StokesI_header):
    # Niquist sampling says: sampling distance is no larger than half the size of the beam

    # Extract beam sizes in degrees from header and convert to arcseconds
    bmaj_arcsec = deg_to_arcsec(StokesI_header["BMAJ"])
    bmin_arcsec = deg_to_arcsec(StokesI_header["BMIN"])

    # Extract pixel scale (CDELT1) and convert to arcseconds/pixel
    pixscale_arcsec = deg_to_arcsec(abs(StokesI_header["CDELT1"]))  # arcsec/pixel

    # Convert beam size to pixels
    bmaj_pix = bmaj_arcsec / pixscale_arcsec
    bmin_pix = bmin_arcsec / pixscale_arcsec
    beam_pix = np.sqrt(bmaj_pix * bmin_pix)

    # Sampling recommendation
    sampling_step = round(beam_pix/2)
    print(f"Recommended sampling step: every {sampling_step} pixels")
# -------------------------------------------------------------------------------------------- 



def generate_polarization_vectors(ny, nx,
                                  xmin, xmax, ymin, ymax, 
                                  RA_centre_pix, Dec_centre_pix,
                                  uniform_angle,
                                  StokesI_mJy, 
                                  POLI_mJy, POLI_err_mJy,
                                  PA_real_sky_rad, PA_err_deg,
                                  band, 
                                  step = None, vector_len_pix = None):
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
    vector_data_actual_cartesian, vector_angle_actual_sky, vector_mask, in_plot_mask = make_vectors(ny, nx,  
                                                                                      xmin, xmax, ymin, ymax, 
                                                                                      POLI_mJy, POLI_err_mJy,
                                                                                      PA_real_sky_rad, PA_err_deg,
                                                                                      band, step, vector_len_pix)
    
    # Make the PA grids for uniform and Azimuthal
    PA_grid_100Uniform   = make_PA_grid_100Uniform(ny,   nx, uniform_angle)
    PA_grid_100Azimuthal = make_PA_grid_100Azimuthal(ny, nx, RA_centre_pix, Dec_centre_pix)  
    
    # Get the vector and angle data for the 100 Uniform case 
    # Currently not saving the mask here
    vector_data_100Uniform_cartesian, vector_angle_100Uniform_sky, _, _ = make_vectors(ny, nx,  
                                                                                    0, 0, 0, 0, # We dont need these here 
                                                                                    POLI_mJy, POLI_err_mJy,
                                                                                    PA_grid_100Uniform, PA_err_deg,
                                                                                    band, step, vector_len_pix)
    
    # Get the vector and angle data for the 100 Azimuthal case 
    # Currently not saving the mask here
    vector_data_100Azimuthal_cartesian, vector_angle_100Azimuthal_sky, _, _ = make_vectors(ny, nx, 
                                                                                        0, 0, 0, 0, # We dont need these here 
                                                                                        POLI_mJy, POLI_err_mJy,
                                                                                        PA_grid_100Azimuthal, PA_err_deg,
                                                                                        band, step, vector_len_pix)
    
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
        'StokesU_grid_100Azimuthal': StokesU_grid_100Azimuthal,
        'vector_mask': vector_mask,
        'in_plot_mask': in_plot_mask
    }
    
    return results





# --------------------------------------------------------------------------------------
def save_beam_info(header, band, print_statements=False):

    if band == 'Band 4':
        path = constants.band4_data_folder_path + "beam_BAND4.csv"
    elif band == 'Band 4 nterms2':
        path = constants.band4_nterms2_data_folder_path + "beam_BAND4_nterms2.csv"
    elif band == 'Band 5':
        path = constants.band5_data_folder_path + "beam_BAND5.csv"
    elif band == 'Band 6':
        path = constants.band6_data_folder_path + "beam_BAND6.csv"
    elif band == 'Band 7 nterms2':
        path = constants.band7_nterms2_data_folder_path + "beam_BAND7_nterms2.csv"
    else:
        raise ValueError("Unsupported band.")

    # Extract beam parameters (degrees)
    BMAJ_deg = header['BMAJ']
    BMIN_deg = header['BMIN']
    BPA_deg  = header['BPA']

    # Convert to arcsec
    BMAJ_arcsec = BMAJ_deg * 3600
    BMIN_arcsec = BMIN_deg * 3600

    # Pixel scale (arcsec/pixel)
    CDELT1_deg = abs(header['CDELT1'])
    pixel_scale_arcsec = CDELT1_deg * 3600
    
    # FWHM to sigma 
    # --------------------------------------------
    # FWHM = 2 sqrt(2 ln 2) sigma 
    # sigma = FWHM /2 sqrt(2 ln(2)
    # --------------------------------------------
    fwhm_to_sigma = 1 / (2 * np.sqrt(2 * np.log(2)))
    
    sigma = BMIN_arcsec * fwhm_to_sigma
    # --------------------------------------------
    

    data_dict = {
        "BMAJ_arcsec": [BMAJ_arcsec],
        "BMIN_arcsec": [BMIN_arcsec],
        "BPA_deg": [BPA_deg],
        "pixel_scale_arcsec": [pixel_scale_arcsec],
        "sigma": [sigma], 
    }

    df = pd.DataFrame(data_dict)
    df.to_csv(path, index=False)

    if print_statements:
        print(f"Saved beam info for {band}")
        print(f"BMAJ: {BMAJ_arcsec:.3f} arcsec")
        print(f"BMIN: {BMIN_arcsec:.3f} arcsec")
        print(f"BPA : {BPA_deg:.2f} deg")
        print(f"Pixel scale: {pixel_scale_arcsec:.4f} arcsec/pixel")
        print(f"sigma: {sigma:.4f} arcsec/pixel")

# --------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------
def load_beam_info(bands, print_things = True):

    results = {
        'Band': [],
        'BMAJ_arcsec': [],
        'BMIN_arcsec': [],
        'BPA_deg': [],
        'pixel_scale_arcsec': [],
        'sigma': []
    }
    
    
    

    for band in bands:
        if print_things:
            print(rf'Working on {band}')


        if band == "Band 4 nterms2" or band == 42:
            var_name = "band4_nterms2_data_folder_path"
            path = Path(getattr(constants, var_name)) / "beam_BAND4_nterms2.csv"

        elif band == "Band 7 nterms2" or band == 72:
            var_name = "band7_nterms2_data_folder_path"
            path = Path(getattr(constants, var_name)) / "beam_BAND7_nterms2.csv"
            
        elif band == "Band 5" or band == 5:
            var_name = "band5_data_folder_path"
            path = Path(getattr(constants, var_name)) / "beam_BAND5.csv"
        
        elif band == "Band 6" or band == 6:
            var_name = "band6_data_folder_path"
            path = Path(getattr(constants, var_name)) / "beam_BAND6.csv"


        else:
            return("band input not accepted")
        df = pd.read_csv(path)

        results['Band'].append(band)
        results['BMAJ_arcsec'].append(df.at[0, "BMAJ_arcsec"])
        results['BMIN_arcsec'].append(df.at[0, "BMIN_arcsec"])
        results['BPA_deg'].append(df.at[0, "BPA_deg"])
        results['pixel_scale_arcsec'].append(df.at[0, "pixel_scale_arcsec"])
        results['sigma'].append(df.at[0, "sigma"])

    df_beam = pd.DataFrame(results)

    return df_beam
# --------------------------------------------------------------------------------------



    
    

# def generate_polarization_vectors_band47(ny, nx,
#                                          RA_centre_pix, Dec_centre_pix,
#                                          uniform_angle,
#                                          StokesI_mJy, 
#                                          POLI_mJy, POLI_err_mJy,
#                                          PA_real_sky_rad, PA_err_deg,
#                                          band, 
#                                          step = None, vector_len_pix = None):
#     """
#     Generates polarization vectors for different grid configurations and calculates the Stokes U and Q grids.

#     Parameters:
#     - ny, nx: Dimensions of the grid
#     - RA_centre_pix, Dec_centre_pix: Centre coordinates in pixels
#     - uniform_angle: The uniform angle for PA grid
#     - POLI_mJy, POLI_err_mJy: Polarization values and their uncertainties
#     - PA_real_sky_rad, PA_err_deg: Polarization angles (in radians) and their uncertainties

#     Returns:
#     - A dictionary containing the vector and angle data for different grid configurations, as well as the Stokes Q and U grids.
#     """
    
    
    
    
#     # Get vector and angles for the actual data
#     vector_data_actual_cartesian, vector_angle_actual_sky = make_vectors_band47(ny, nx,  
#                                                                                 POLI_mJy, POLI_err_mJy,
#                                                                                 PA_real_sky_rad, PA_err_deg,
#                                                                                 band, step, vector_len_pix)
    
#     # Make the PA grids for uniform and Azimuthal
#     PA_grid_100Uniform   = make_PA_grid_100Uniform(ny,   nx, uniform_angle)
#     PA_grid_100Azimuthal = make_PA_grid_100Azimuthal(ny, nx, RA_centre_pix, Dec_centre_pix)  
    
#     # Get the vector and angle data for the 100 Uniform case 
#     vector_data_100Uniform_cartesian, vector_angle_100Uniform_sky = make_vectors_band47(ny, nx,  
#                                                                                        POLI_mJy, POLI_err_mJy,
#                                                                                        PA_grid_100Uniform, PA_err_deg,
#                                                                                        band, step, vector_len_pix)
    
#     # Get the vector and angle data for the 100 Azimuthal case 
#     vector_data_100Azimuthal_cartesian, vector_angle_100Azimuthal_sky = make_vectors_band47(ny, nx,  
#                                                                                            POLI_mJy, POLI_err_mJy,
#                                                                                            PA_grid_100Azimuthal, PA_err_deg,
#                                                                                            band, step, vector_len_pix)
    
#     # Get Stokes Q and U grids
#     StokesQ_grid_100Uniform,   StokesU_grid_100Uniform   = recover_StokesQU(PA_grid_100Uniform,   StokesI_mJy, ny, nx)
#     StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal = recover_StokesQU(PA_grid_100Azimuthal, StokesI_mJy, ny, nx)
    
#     # Organize the results in a dictionary for easy access
#     results = {
#         'vector_data_actual_cartesian': vector_data_actual_cartesian,
#         'vector_angle_actual_sky': vector_angle_actual_sky,
#         'vector_data_100Uniform_cartesian': vector_data_100Uniform_cartesian,
#         'vector_angle_100Uniform_sky': vector_angle_100Uniform_sky,
#         'vector_data_100Azimuthal_cartesian': vector_data_100Azimuthal_cartesian,
#         'vector_angle_100Azimuthal_sky': vector_angle_100Azimuthal_sky,
#         'StokesQ_grid_100Uniform': StokesQ_grid_100Uniform,
#         'StokesU_grid_100Uniform': StokesU_grid_100Uniform,
#         'StokesQ_grid_100Azimuthal': StokesQ_grid_100Azimuthal,
#         'StokesU_grid_100Azimuthal': StokesU_grid_100Azimuthal
#     }
    
#     return results






# def generate_polarization_vectors_band6(ny, nx,
#                                         RA_centre_pix, Dec_centre_pix,
#                                         uniform_angle,
#                                         StokesI_mJy, StokesI_err_mJy,
#                                         POLI_mJy, POLI_err_mJy,
#                                         PA_real_sky_rad, PA_err_deg, step = None):
#     """
#     Generates polarization vectors for different grid configurations and calculates the Stokes U and Q grids.

#     Parameters:
#     - ny, nx: Dimensions of the grid
#     - RA_centre_pix, Dec_centre_pix: Centre coordinates in pixels
#     - uniform_angle: The uniform angle for PA grid
#     - StokesI_mJy, StokesI_err_mJy: Stokes I values and their uncertainties
#     - POLI_mJy, POLI_err_mJy: Polarization values and their uncertainties
#     - PA_real_sky_rad, PA_err_deg: Polarization angles (in radians) and their uncertainties

#     Returns:
#     - A dictionary containing the vector and angle data for different grid configurations, as well as the Stokes Q and U grids.
#     """
    
#     # Get vector and angles for the actual data
#     vector_data_actual_cartesian, vector_angle_actual_sky = make_vectors_band6(ny, nx,  
#                                                                                 StokesI_mJy, StokesI_err_mJy,
#                                                                                 POLI_mJy, POLI_err_mJy,
#                                                                                 PA_real_sky_rad, PA_err_deg)
    
#     # Make the PA grids for uniform and Azimuthal
#     PA_grid_100Uniform = make_PA_grid_100Uniform(ny, nx, uniform_angle)
#     PA_grid_100Azimuthal = make_PA_grid_100Azimuthal(ny, nx, RA_centre_pix, Dec_centre_pix)  
    
#     # Get the vector and angle data for the 100 Uniform case 
#     vector_data_100Uniform_cartesian, vector_angle_100Uniform_sky = make_vectors_band6(ny, nx,  
#                                                                                          StokesI_mJy, StokesI_err_mJy,
#                                                                                          POLI_mJy, POLI_err_mJy,
#                                                                                          PA_grid_100Uniform, PA_err_deg)
    
#     # Get the vector and angle data for the 100 Azimuthal case 
#     vector_data_100Azimuthal_cartesian, vector_angle_100Azimuthal_sky = make_vectors_band6(ny, nx,  
#                                                                                              StokesI_mJy, StokesI_err_mJy,
#                                                                                              POLI_mJy, POLI_err_mJy,
#                                                                                              PA_grid_100Azimuthal, PA_err_deg)
    
#     # Get Stokes Q and U grids
#     StokesQ_grid_100Uniform, StokesU_grid_100Uniform = recover_StokesQU(PA_grid_100Uniform, StokesI_mJy, ny, nx)
#     StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal = recover_StokesQU(PA_grid_100Azimuthal, StokesI_mJy, ny, nx)
    
#     # Organize the results in a dictionary for easy access
#     results = {
#         'vector_data_actual_cartesian': vector_data_actual_cartesian,
#         'vector_angle_actual_sky': vector_angle_actual_sky,
#         'vector_data_100Uniform_cartesian': vector_data_100Uniform_cartesian,
#         'vector_angle_100Uniform_sky': vector_angle_100Uniform_sky,
#         'vector_data_100Azimuthal_cartesian': vector_data_100Azimuthal_cartesian,
#         'vector_angle_100Azimuthal_sky': vector_angle_100Azimuthal_sky,
#         'StokesQ_grid_100Uniform': StokesQ_grid_100Uniform,
#         'StokesU_grid_100Uniform': StokesU_grid_100Uniform,
#         'StokesQ_grid_100Azimuthal': StokesQ_grid_100Azimuthal,
#         'StokesU_grid_100Azimuthal': StokesU_grid_100Azimuthal
#     }
    
#     return results
 
    
    

    
    
    

   
    
    