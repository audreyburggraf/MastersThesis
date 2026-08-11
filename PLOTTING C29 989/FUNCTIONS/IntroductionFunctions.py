# Import constants
# ---------------------------------------------------------------------------------------
import sys
import pandas as pd
from pathlib import Path
import os

# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")

# Now you can import constants.py
import constants
# ---------------------------------------------------------------------------------------


# Load in the function:
# ------------------------------------------
from FITS_Image_Functions import * 
from PolarizationFunctions import *
from UnitConversion import * 
# ------------------------------------------


def get_plotting_parameters(StokesI_header, StokesI_wcs, band):
    # Dictionary for band constants
    # Get the constants depending on what band of data we are using:
    centre_lookup = {
        'Band 6': constants.centre_str_band6,
        'Band 6 smooth': constants.centre_str_band6_smooth,
        'Band 6 smooth B7': constants.centre_str_band6_smooth_B7,

        'Band 4': constants.centre_str_band4,
        'Band 4 nterms2': constants.centre_str_band4_nterms2,
        'Band 4 nterms2 robust -1': constants.centre_str_band4_nterms2_robust_minus1,
        'Band 4 nterms2 smooth': constants.centre_str_band4_nterms2_smooth,
        'Band 4 nterms2 smooth B6': constants.centre_str_band4_nterms2_smooth_B6,
        'Band 4 nterms2 smooth B6 B7': constants.centre_str_band4_nterms2_smooth_B6_B7,

        'Band 5 v0': constants.centre_str_band5_v0,
        'Band 5': constants.centre_str_band5,
        'Band 5 robust -2': constants.centre_str_band5_robust_minus2,
        'Band 5 robust -1': constants.centre_str_band5_robust_minus1,
        'Band 5 nterms2': constants.centre_str_band5_nterms2,

        'Band 7 nterms2': constants.centre_str_band7_nterms2,
        'Band 7 nterms2 smooth': constants.centre_str_band7_nterms2_smooth,
        'Band 7 nterms2 smooth B6': constants.centre_str_band7_nterms2_smooth
    }

    centre_str = centre_lookup.get(band)

    if centre_str is None:
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
    
    print("band =", band)
    print("max_str =", repr(max_str))
    print("StokesI_wcs CRPIX =", StokesI_wcs.wcs.crpix)
    print("StokesI_wcs CDELT =", StokesI_wcs.wcs.cdelt)
    
    pixel_scale_arcsec = abs(StokesI_wcs.wcs.cdelt[0]) * 3600
    print(f"Pixel scale = {pixel_scale_arcsec:.3f} arcsec/pixel")

    print("Dec_max_pix =", Dec_max_pix)
    print(" ")


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
def print_errors(band, StokesI_mJy, StokesI_err_mJy, StokesQ_err_mJy, StokesU_err_mJy, POLI_err_mJy, POLF_err, PA_err_deg):
    
    # Find max Stokes I 
    ymax, xmax = np.unravel_index(np.nanargmax(StokesI_mJy), StokesI_mJy.shape)
    
    const_err_vars = ['Stokes I', 'Stokes Q', 'Stokes U', 'POLI']
    for i, SV in enumerate([StokesI_err_mJy, StokesQ_err_mJy, StokesU_err_mJy, POLI_err_mJy]):

        # Stokes errors and debiased POLI are constant for Band 4, 5, 7
        if band in (4, 5, 7):
            SV_err_mJy = SV[0,0]
            SV_err_uJy = milli_to_micro(SV_err_mJy)
            print(rf'{const_err_vars[i]} err: {SV_err_mJy:.4f} mJy, or {SV_err_uJy:.2f} micro Jy')
         
        # Band 6 they are not all the same 
        else:
            SV_err_uJy = milli_to_micro(SV)
            print(rf'{const_err_vars[i]} error: mean: {np.nanmean(SV_err_uJy):.2f}, maximum: {np.nanmax(SV_err_uJy):.2f}, at Stokes I maximum: {SV_err_uJy[ymax, xmax]:.2f}')
        
    
    # Polarization fraction error is NOT the same 
    POLF_err_percent = POLF_err * 100
    print(rf'POLF error: mean: {np.nanmean(POLF_err_percent):.4f}, maximum: {np.nanmax(POLF_err_percent):.4f}, at Stokes I maximum: {POLF_err_percent[ymax, xmax]:.4f}, ALL MULTIPLIED BY 100')
    
    
    # Polarization fraction error is NOT the same 
    print(rf'PA error: mean: {np.nanmean(PA_err_deg):.4f}, maximum: {np.nanmax(PA_err_deg):.4f}, at Stokes I maximum: {PA_err_deg[ymax, xmax]:.4f}')
# -------------------------------------------------------------------------------------------- 



def generate_polarization_vectors(ny, nx,
                                  xmin, xmax, ymin, ymax, 
                                  RA_centre_pix, Dec_centre_pix,
                                  uniform_angle,
                                  StokesI_mJy, StokesI_err_mJy, 
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
    vector_data_actual_cartesian, vector_angle_actual_sky, vector_angle_actual_sky_errors, vector_mask, in_plot_mask = make_vectors_StokesIcutoff(ny, nx,  
                                                                                      xmin, xmax, ymin, ymax, 
                                                                                                                                                  StokesI_mJy, StokesI_err_mJy, 
                                                                                      POLI_mJy, POLI_err_mJy,
                                                                                      PA_real_sky_rad, PA_err_deg,
                                                                                      band, step, vector_len_pix)
    
    
    # Make the PA grids for uniform and Azimuthal
    PA_grid_100Uniform   = make_PA_grid_100Uniform(ny,   nx, uniform_angle)
    PA_grid_100Azimuthal = make_PA_grid_100Azimuthal(ny, nx, RA_centre_pix, Dec_centre_pix)  
    
    # Get the vector and angle data for the 100 Uniform case 
    # Currently not saving the mask here
    vector_data_100Uniform_cartesian, vector_angle_100Uniform_sky, _, _, _ = make_vectors_StokesIcutoff(ny, nx,  
                                                                                    0, 0, 0, 0, # We dont need these here 
                                                                                    StokesI_mJy, StokesI_err_mJy, 
                                                                                    POLI_mJy, POLI_err_mJy,
                                                                                    PA_grid_100Uniform, PA_err_deg,
                                                                                    band, step, vector_len_pix)
    
    # Get the vector and angle data for the 100 Azimuthal case 
    # Currently not saving the mask here
    vector_data_100Azimuthal_cartesian, vector_angle_100Azimuthal_sky,_,  _, _ = make_vectors_StokesIcutoff(ny, nx, 
                                                                                        0, 0, 0, 0, # We dont need these here 
                                                                                        StokesI_mJy, StokesI_err_mJy, 
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
        'vector_angle_actual_sky_errors': vector_angle_actual_sky_errors,
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

    beam_paths = {
        'Band 4': (constants.band4_data_folder_path, "beam_BAND4.csv"),
        'Band 4 nterms2': (constants.band4_nterms2_data_folder_path, "beam_BAND4_nterms2.csv"),
        'Band 4 nterms2 robust -1': (constants.band4_nterms2_robust_minus1_data_folder_path, 
                                     "beam_BAND4_nterms2_robust_minus1.csv"),
        'Band 4 nterms2 smooth': (constants.band4_nterms2_smooth_data_folder_path, "beam_BAND4_nterms2_smooth.csv"),
        'Band 4 nterms2 smooth B6': (constants.band4_nterms2_smooth_B6_data_folder_path, "beam_BAND4_nterms2_smooth_B6.csv"),
        'Band 4 nterms2 smooth B6 B7': (constants.band4_nterms2_smooth_B6_B7_data_folder_path, 
                                        "beam_BAND4_nterms2_smooth_B6_B7.csv"),

        'Band 5': (constants.band5_data_folder_path, "beam_BAND5.csv"),
        'Band 5 v0': (constants.band5_v0_data_folder_path, "beam_BAND5_v0.csv"),
        'Band 5 robust -2': (constants.band5_robust_minus2_data_folder_path, "beam_BAND5_robust_minus2.csv"),
        'Band 5 robust -1': (constants.band5_robust_minus1_data_folder_path, "beam_BAND5_robust_minus1.csv"),
        'Band 5 nterms2': (constants.band5_nterms2_data_folder_path, "beam_BAND5_nterms2.csv"),

        'Band 6': (constants.band6_data_folder_path, "beam_BAND6.csv"),
        'Band 6 smooth': (constants.band6_smooth_data_folder_path, "beam_BAND6_smooth.csv"),
        'Band 6 smooth B7': (constants.band6_smooth_B7_data_folder_path, "beam_BAND6_smooth_B7.csv"),

        'Band 7 nterms2': (constants.band7_nterms2_data_folder_path, "beam_BAND7_nterms2.csv"),
        'Band 7 nterms2 smooth': (constants.band7_nterms2_smooth_data_folder_path, "beam_BAND7_nterms2_smooth.csv"),
        'Band 7 nterms2 smooth B6': (constants.band7_nterms2_smooth_B6_data_folder_path, "beam_BAND7_nterms2_smooth_B6.csv"),
    }

    if band not in beam_paths:
        raise ValueError(f"Unsupported band: {band}")

    folder, filename = beam_paths[band]
    path = Path(folder) / filename

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
    
    
    beam_file_lookup = {
        "Band 4 nterms2": ("band4_nterms2_data_folder_path", "beam_BAND4_nterms2.csv"),
        "Band 4 nterms2 robust -1": ("band4_nterms2_robust_minus1_data_folder_path", "beam_BAND4_nterms2_robust_minus1.csv"),
        "Band 4 nterms2 smooth": ("band4_nterms2_smooth_data_folder_path", "beam_BAND4_nterms2_smooth.csv"),
        "Band 4 nterms2 smooth B6": ("band4_nterms2_smooth_B6_data_folder_path", "beam_BAND4_nterms2_smooth_B6.csv"),
        "Band 4 nterms2 smooth B6 B7": ("band4_nterms2_smooth_B6_B7_data_folder_path", "beam_BAND4_nterms2_smooth_B6_B7.csv"),

        "Band 7 nterms2": ("band7_nterms2_data_folder_path", "beam_BAND7_nterms2.csv"),
        "Band 7 nterms2 smooth": ("band7_nterms2_smooth_data_folder_path", "beam_BAND7_nterms2_smooth.csv"),
        "Band 7 nterms2 smooth B6": ("band7_nterms2_smooth_B6_data_folder_path", "beam_BAND7_nterms2_smooth_B6.csv"),

        "Band 5": ("band5_data_folder_path", "beam_BAND5.csv"),
        "Band 5 v0": ("band5_v0_data_folder_path", "beam_BAND5_v0.csv"),
        "Band 5 nterms2": ("band5_nterms2_data_folder_path", "beam_BAND5_nterms2.csv"),
        "Band 5 robust -1": ("band5_robust_minus1_data_folder_path", "beam_BAND5_robust_minus1.csv"),
        "Band 5 robust -2": ("band5_robust_minus2_data_folder_path", "beam_BAND5_robust_minus2.csv"),

        "Band 6": ("band6_data_folder_path", "beam_BAND6.csv"),
        "Band 6 smooth": ("band6_smooth_data_folder_path", "beam_BAND6_smooth.csv"),
        "Band 6 smooth B7": ("band6_smooth_B7_data_folder_path", "beam_BAND6_smooth_B7.csv"),
    }


    for band in bands:

        if print_things:
            print(f"Working on {band}")

        if band not in beam_file_lookup:
            raise ValueError(f"Band input not accepted: {band}")

        var_name, filename = beam_file_lookup[band]

        path = Path(getattr(constants, var_name)) / filename

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
 
    
    

# ------------------------------------------------------------
# Function to build band dictionary
# ------------------------------------------------------------
def BuildBandDictionary(
    band_data_folder_path,
    StokesI_wcs,
    xmin,
    xmax,
    ymin,
    ymax,
    reference_length_pix,
    reference_length_AU,
    BMAJ_pix,
    BMIN_pix,
    BPA_deg_cartesian,
    max_length_pix,
    reference_fraction,
    POLI_debiased_mJy,
    StokesI_mJy,
    RA_centre_pix, Dec_centre_pix,
):

    # --------------------------------------------------------
    # Load saved files
    # --------------------------------------------------------
    ratios = np.load(
        os.path.join(band_data_folder_path,
                     "GaussianUniformRatios_best.npy")
    )

    VectorsActual = np.load(
        os.path.join(band_data_folder_path,
                     "vector_data_actual_cartesian.npy")
    )
    VectorsGaussian = np.load(
        os.path.join(band_data_folder_path,
                     "vector_data_gaussian_best.npy")
    )
    
    VectorsRatio = np.load(
        os.path.join(band_data_folder_path,
                     "vector_data_ratio_best.npy")
    )

    params = pd.read_csv(
        os.path.join(band_data_folder_path,
                     "best_gaussian_params.csv")
    )
    
    slice_points = np.load(
        os.path.join(band_data_folder_path, "slice_points.npy"),
        allow_pickle=True,
    ).item()

    # --------------------------------------------------------
    # Make dictionary
    # --------------------------------------------------------
    BandDict = {

        'band_data_folder_path': band_data_folder_path,
        
        'StokesI_wcs': StokesI_wcs,

        'ratios': ratios,

        'xmin': xmin,
        'xmax': xmax,
        'ymin': ymin,
        'ymax': ymax,

        'reference_length_pix': reference_length_pix,
        'reference_length_AU': reference_length_AU,

        'BMAJ_pix': BMAJ_pix,
        'BMIN_pix': BMIN_pix,
        'BPA_deg_cartesian': BPA_deg_cartesian,

        'max_length_pix': max_length_pix,
        'reference_fraction': reference_fraction,

        'phi': params["phi"].iloc[0],
        'BMAJ': params["BMAJ"].iloc[0],
        'BMIN': params["BMIN"].iloc[0],
        'chi2': params["chi2"].iloc[0],

        'VectorsActual': VectorsActual,
        'VectorsGaussian': VectorsGaussian,
        'VectorsRatio': VectorsRatio,

        'POLI_debiased_mJy': POLI_debiased_mJy,
        'StokesI_mJy': StokesI_mJy,
        
        'minor_x': slice_points['minor_x'],
        'minor_y': slice_points['minor_y'], 
        
        'RA_centre_pix': RA_centre_pix, 
        'Dec_centre_pix': Dec_centre_pix, 
    }

    return BandDict    
    
    

   
    
    