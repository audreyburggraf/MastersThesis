import numpy as np

# Functions
from FITS_Image_Functions import *
from MakingGridFunctions import * 


from astropy.convolution import convolve, Gaussian2DKernel
# from radio_beam import EllipticalGaussian2DKernel
from tqdm.notebook import tqdm

def calculate_polarized_intensity(StokesQ, StokesU):
    """
    Calculate the polarized intensity.

    Parameters:
        StokesQ (2D array): Stokes Q data in mJy or similar units.
        StokesU (2D array): Stokes U data in mJy or similar units.

    Returns:
        2D array: Polarized intensity in the same units as the inputs.
    """
    return np.sqrt(StokesQ**2 + StokesU**2)


def calculate_polarized_fraction(StokesQ, StokesU, StokesI):
    """
    Calculate the polarized fraction.

    Parameters:
        StokesQ (2D array): Stokes Q data in mJy or similar units.
        StokesU (2D array): Stokes U data in mJy or similar units.
        StokesI (2D array): Stokes I data in mJy or similar units.

    Returns:
        2D array: Polarized fraction (unitless).
    """
    polarized_intensity = calculate_polarized_intensity(StokesQ, StokesU)
    
    return polarized_intensity / StokesI


def calculate_polarization_angle(StokesQ, StokesU):
    """
    Calculate the polarization angle in radians.

    Parameters:
        StokesQ (2D array): Stokes Q data in Jy or similar units.
        StokesU (2D array): Stokes U data in Jy or similar units.

    Returns:
        2D array: Polarization angle in radians.
    """
    return 0.5 * np.arctan2(StokesU, StokesQ)



def recover_StokesQU(PA_grid_rad_astronomy, StokesI_grid, ny, nx):
    """
    Compute Stokes Q and U given a polarization angle grid 

    Parameters:
    PA_grid (numpy.ndarray): 2D array of polarization angles in radians.

    Returns:
    tuple: (StokesQ_grid, StokesU_grid), both as 2D numpy arrays.
    """

    # Stokes Q and U calculations (vectorized)
    p_frac = 0.02  # Ensure this is the correct polarization fraction
    
    StokesQ_grid = p_frac * StokesI_grid * np.cos(2 * PA_grid_rad_astronomy)
    StokesU_grid = p_frac * StokesI_grid * np.sin(2 * PA_grid_rad_astronomy)

    return StokesQ_grid, StokesU_grid




def create_uniform_PA_grid(nx, ny, angle_rad_astronomy):
    """
    Create a 2D grid (ny x nx) where every element is assigned the same angle value.

    Parameters:
    nx (int): Number of columns (x-dimension).
    ny (int): Number of rows (y-dimension).
    angle_rad (float): The uniform angle value in radians.

    Returns:
    np.ndarray: A 2D array filled with the given angle.
    """
    return np.full((ny, nx), angle_rad_astronomy)






def mix_StokesQU_and_generate_vectors(Uniform_ratio, Azimuthal_ratio, 
                                      StokesQ_uniform, StokesU_uniform, 
                                      StokesQ_azimuthal, StokesU_azimuthal,
                                      ny, nx, 
                                      step, vector_length_pix_const, 
                                      StokesI_data_2d_mJy, StokesIerr_data_2d_mJy,
                                      calculated_polarized_intensity, PolarizedIntensity_err_data_2d_mJy,
                                      PolarizationAngle_err_data_2d_deg):
    """
    Blend Stokes Q and U grids using specified ratios, compute polarization angle (PA), 
    and generate mixed polarization vectors based on selection criteria.

    Parameters:
    Uniform_ratio (float): Fraction of uniform component (between 0 and 1).
    Azimuthal_ratio (float): Fraction of azimuthal component (between 0 and 1). 
                             Must satisfy Uniform_ratio + Azimuthal_ratio = 1.
    StokesQ_uniform (numpy.ndarray): 2D array of Stokes Q values from the uniform grid.
    StokesQ_azimuthal (numpy.ndarray): 2D array of Stokes Q values from the azimuthal grid.
    StokesU_uniform (numpy.ndarray): 2D array of Stokes U values from the uniform grid.
    StokesU_azimuthal (numpy.ndarray): 2D array of Stokes U values from the azimuthal grid.

    Returns:
    tuple: (PA_grid_mixed, StokesQ_grid_mixed, StokesU_grid_mixed, vector_mixed)
    """    
    # Ensure the ratios sum to 1 (for safety)
    #if not np.isclose(Uniform_ratio + Azimuthal_ratio, 1.0):
        #raise ValueError("Uniform_ratio and Azimuthal_ratio must sum to 1.")

    # Compute weighted sum of Stokes Q and U
    StokesQ_grid_mixed = Uniform_ratio * StokesQ_uniform + Azimuthal_ratio * StokesQ_azimuthal
    StokesU_grid_mixed = Uniform_ratio * StokesU_uniform + Azimuthal_ratio * StokesU_azimuthal

    # Compute polarization angle (theta = 1/2 * arctan(U/Q))
    PA_grid_mixed_rad_astronomy = calculate_polarization_angle(StokesQ_grid_mixed, StokesU_grid_mixed)
                              # = 0.5 * np.arctan2(StokesU_grid_mixed, StokesQ_grid_mixed)
    
    #PA_grid_mixed_rad_astronomy = castesian_to_astronomy_rad(PA_grid_mixed_rad_cartesian)

    # Generate polarization vectors based on selection criteria
    vector_mixed_cartesian = []
    vector_mixed_angle_rad_astronomy = []
    
    for x in range(0, nx, step):
        for y in range(0, ny, step):
            if (
                StokesI_data_2d_mJy[y, x] / StokesIerr_data_2d_mJy[y, x] > 3
                and calculated_polarized_intensity[y, x] / PolarizedIntensity_err_data_2d_mJy[y, x] > 3
                and PolarizationAngle_err_data_2d_deg[y, x] < 10
            ):
                # Extract the polarization angle at this location
                PA_angle_rad_astronomoy = PA_grid_mixed_rad_astronomy[y, x] 
                
                # Compute vector components
                dx = vector_length_pix_const * np.cos(PA_angle_rad_astronomoy + np.pi/2)
                dy = vector_length_pix_const * np.sin(PA_angle_rad_astronomoy + np.pi/2)

                # Append to the list in the format [x_start, x_end, y_start, y_end]
                vector_mixed_cartesian.append([x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2])
                
                vector_mixed_angle_rad_astronomy.append(PA_angle_rad_astronomoy)

    return PA_grid_mixed_rad_astronomy, StokesQ_grid_mixed, StokesU_grid_mixed, vector_mixed_cartesian, vector_mixed_angle_rad_astronomy






# def mix_StokesQU_and_generate_vectors_band4(Uniform_ratio, Azimuthal_ratio, 
#                                       StokesQ_uniform, StokesU_uniform, 
#                                       StokesQ_azimuthal, StokesU_azimuthal,
#                                       ny, nx, 
#                                       step, vector_length_pix_const, 
#                                       StokesI_data_2d_mJy, StokesIerr_data_2d_mJy,
#                                       #calculated_polarized_intensity, PolarizedIntensity_err_data_2d_mJy,
#                                       PolarizationAngle_err_data_2d_deg):
#     """
#     Blend Stokes Q and U grids using specified ratios, compute polarization angle (PA), 
#     and generate mixed polarization vectors based on selection criteria.

#     Parameters:
#     Uniform_ratio (float): Fraction of uniform component (between 0 and 1).
#     Azimuthal_ratio (float): Fraction of azimuthal component (between 0 and 1). 
#                              Must satisfy Uniform_ratio + Azimuthal_ratio = 1.
#     StokesQ_uniform (numpy.ndarray): 2D array of Stokes Q values from the uniform grid.
#     StokesQ_azimuthal (numpy.ndarray): 2D array of Stokes Q values from the azimuthal grid.
#     StokesU_uniform (numpy.ndarray): 2D array of Stokes U values from the uniform grid.
#     StokesU_azimuthal (numpy.ndarray): 2D array of Stokes U values from the azimuthal grid.

#     Returns:
#     tuple: (PA_grid_mixed, StokesQ_grid_mixed, StokesU_grid_mixed, vector_mixed)
#     """    
#     # Ensure the ratios sum to 1 (for safety)
#     #if not np.isclose(Uniform_ratio + Azimuthal_ratio, 1.0):
#         #raise ValueError("Uniform_ratio and Azimuthal_ratio must sum to 1.")

#     # Compute weighted sum of Stokes Q and U
#     StokesQ_grid_mixed = Uniform_ratio * StokesQ_uniform + Azimuthal_ratio * StokesQ_azimuthal
#     StokesU_grid_mixed = Uniform_ratio * StokesU_uniform + Azimuthal_ratio * StokesU_azimuthal

#     # Compute polarization angle (theta = 1/2 * arctan(U/Q))
#     PA_grid_mixed_rad_astronomy = calculate_polarization_angle(StokesQ_grid_mixed, StokesU_grid_mixed)
#                               # = 0.5 * np.arctan2(StokesU_grid_mixed, StokesQ_grid_mixed)
    
#     #PA_grid_mixed_rad_astronomy = castesian_to_astronomy_rad(PA_grid_mixed_rad_cartesian)

#     # Generate polarization vectors based on selection criteria
#     vector_mixed_cartesian = []
#     vector_mixed_angle_rad_astronomy = []
    
#     for x in range(0, nx, step):
#         for y in range(0, ny, step):
#             if (
#                 # Look at S/N of PI
#                 StokesI_data_2d_mJy[y, x] / StokesIerr_data_2d_mJy[y, x] > 3
#                 # and calculated_polarized_intensity[y, x] / PolarizedIntensity_err_data_2d_mJy[y, x] > 3
#                 and PolarizationAngle_err_data_2d_deg[y, x] < 10
#             ):
#                 # Extract the polarization angle at this location
#                 PA_angle_rad_astronomoy = PA_grid_mixed_rad_astronomy[y, x] 
                
#                 # Compute vector components
#                 dx = vector_length_pix_const * np.cos(PA_angle_rad_astronomoy + np.pi/2)
#                 dy = vector_length_pix_const * np.sin(PA_angle_rad_astronomoy + np.pi/2)

#                 # Append to the list in the format [x_start, x_end, y_start, y_end]
#                 vector_mixed_cartesian.append([x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2])
                
#                 vector_mixed_angle_rad_astronomy.append(PA_angle_rad_astronomoy)

#     return PA_grid_mixed_rad_astronomy, StokesQ_grid_mixed, StokesU_grid_mixed, vector_mixed_cartesian, vector_mixed_angle_rad_astronomy





def mix_StokesQU_and_generate_vectors_band4(Uniform_ratio, Azimuthal_ratio, 
                                            StokesQ_uniform, StokesU_uniform, 
                                            StokesQ_azimuthal, StokesU_azimuthal,
                                            ny, nx, 
                                            POLI_mJy, POLI_err_mJy,
                                            PA_err_deg):
    """
    Blend Stokes Q and U grids using specified ratios, compute polarization angle (PA), 
    and generate mixed polarization vectors based on selection criteria.

    Parameters:
    Uniform_ratio (float): Fraction of the uniform component (between 0 and 1).
    Azimuthal_ratio (float): Fraction of the azimuthal component (between 0 and 1). 
                             Must satisfy Uniform_ratio + Azimuthal_ratio = 1.
    StokesQ_uniform (np.ndarray): 2D array of Stokes Q values from the uniform grid.
    StokesU_uniform (np.ndarray): 2D array of Stokes U values from the uniform grid.
    StokesQ_azimuthal (np.ndarray): 2D array of Stokes Q values from the azimuthal grid.
    StokesU_azimuthal (np.ndarray): 2D array of Stokes U values from the azimuthal grid.
    ny (int): Number of pixels along the y-axis.
    nx (int): Number of pixels along the x-axis.
    POLI_mJy (np.ndarray): 2D array of polarization intensity in mJy.
    POLI_err_mJy (np.ndarray): 2D array of polarization intensity errors in mJy.
    PA_err_deg (np.ndarray): 2D array of polarization angle errors in degrees.

    Returns:
    tuple: 
        - PA_grid_mixed_rad_sky (np.ndarray): Polarization angle map (radians, sky frame).
        - StokesQ_grid_mixed (np.ndarray): Mixed Stokes Q map.
        - StokesU_grid_mixed (np.ndarray): Mixed Stokes U map.
        - vector_mixed_cartesian (list): List of vectors in [x0, x1, y0, y1] format.
        - vector_mixed_angle_rad_sky (list): List of vector angles in radians (sky frame).
    """

    # Compute weighted sum of Stokes Q and U
    StokesQ_grid_mixed = Uniform_ratio * StokesQ_uniform + Azimuthal_ratio * StokesQ_azimuthal
    StokesU_grid_mixed = Uniform_ratio * StokesU_uniform + Azimuthal_ratio * StokesU_azimuthal

    # Compute polarization angle (theta = 1/2 * arctan(U/Q))
    PA_grid_mixed_rad_sky = calculate_polarization_angle(StokesQ_grid_mixed, StokesU_grid_mixed)
                             
    # Generate polarization vectors
    vector_mixed_cartesian, vector_mixed_angle_rad_sky = make_vectors_band4(ny, nx, 
                                                                            POLI_mJy, POLI_err_mJy, 
                                                                            PA_grid_mixed_rad_sky, PA_err_deg)

    return PA_grid_mixed_rad_sky, StokesQ_grid_mixed, StokesU_grid_mixed, vector_mixed_cartesian, vector_mixed_angle_rad_sky








