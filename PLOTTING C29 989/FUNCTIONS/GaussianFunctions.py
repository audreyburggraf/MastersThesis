import numpy as np
from scipy.stats import norm


import sys
import os

# Functions
# Define folder path
func_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/FUNCTIONS/"

sys.path.append(func_folder_path)

from FITS_Image_Functions import *
from PolarizationFunctions import *
from PlottingWithFunction import *
from DataAnalysisFunctions import *
from GaussianFunctions import *

from custom_colormap import *


# ----------------------------------------------------------------------------------------
def find_complement_and_total(gaussian):
    """Calculate the complement of a Gaussian value and return (complement, total)."""
    
    complement = 1 - gaussian
    
    total = gaussian + complement
    
    return complement, total
# ----------------------------------------------------------------------------------------




# ----------------------------------------------------------------------------------------
def calculate_gaussian_1d(values, sigma, mean):
    
    #gaussian = norm.pdf(values)
    
    I_0 = 1 / (sigma * np.sqrt(2 * np.pi))
    
    values_shifted = values - mean
    
    gaussian = I_0 * np.exp(-0.5 * np.abs(values_shifted / sigma) ** 2)
    
    
    # Normalize the Gaussian to a maximum of 1
    gaussian /= np.max(gaussian)


    return gaussian
# ----------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------
def calculate_gaussian_2d(x_values, y_values, sigma_x, sigma_y, mean_x, mean_y):
    """
    Computes a 2D Gaussian function.
    
    Parameters:
    - x_values: Array representing the x-axis values
    - y_values: Array representing the y-axis values
    - sigma_x: Standard deviation in the x direction
    - sigma_y: Standard deviation in the y direction
    
    Returns:
    - Normalized 2D Gaussian values
    """
    
    # Create meshgrid for 2D
    X, Y = np.meshgrid(x_values, y_values)
    
    # Center coordinates (optional, assume mean-centered)
    X_shifted, Y_shifted = X - mean_x, Y - mean_y
    
    I_0 = 1 / (sigma_x * sigma_y * np.sqrt(2 * np.pi))
    
    # Compute 2D Gaussian function
    gaussian = I_0 * np.exp(-0.5 * ((X_shifted / sigma_x) ** 2 + (Y_shifted / sigma_y) ** 2))
    
    # Normalize the Gaussian to a maximum of 1
    gaussian /= np.max(gaussian)

    return gaussian
# ----------------------------------------------------------------------------------------





# ----------------------------------------------------------------------------------------
def calculate_gaussian_1d_flat(values, sigma, mean, phi):
    """
    Computes a flat-topped (super-Gaussian) 1D function.
    
    Parameters:
    - values: Array of x-values
    - sigma: Width parameter
    - phi: Shape parameter (higher = flatter top)
    
    Returns:
    - Normalized super-Gaussian values
    """
    
    I_0 = 1 / (sigma * np.sqrt(2 * np.pi))
    
    # Center coordinates (optional, assume mean-centered)
    values_shifted = values - mean
    
    gaussian = I_0 * np.exp(-0.5 * np.abs(values_shifted / sigma) ** phi)

    # Normalize to ensure peak value is 1
    gaussian /= np.max(gaussian)

    return gaussian

# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
def calculate_gaussian_2d_flat(x_values, y_values, sigma_x, sigma_y, mean_x, mean_y, phi):
    """
    Computes a 2D flat-topped (super-Gaussian) function.
    
    Parameters:
    - x_values: Array representing the x-axis values
    - y_values: Array representing the y-axis values
    - sigma_x: Standard deviation in the x direction
    - sigma_y: Standard deviation in the y direction
    - phi: Shape parameter (higher = flatter top)
    
    Returns:
    - Normalized 2D super-Gaussian v,alues
    """
    
    # Create meshgrid for 2D
    X, Y = np.meshgrid(x_values, y_values)
    
    # Center coordinates (optional, assume mean-centered)
    X_shifted, Y_shifted = X - mean_x, Y - mean_y
    
    # Compute the super-Gaussian function using absolute values and phi exponent
    I_0 = 1 / (sigma_x * sigma_y * np.sqrt(2 * np.pi))
    
    gaussian = I_0 * np.exp(-0.5 * (np.abs(X_shifted / sigma_x) ** phi + np.abs(Y_shifted / sigma_y) ** phi))
    
    # Normalize the Gaussian to a maximum of 1
    gaussian /= np.max(gaussian)

    return gaussian

# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
def calculate_gaussian_2d_flat_tilted(x_values, y_values, sigma_x, sigma_y, mean_x, mean_y, phi, theta_rad):
    """
    Computes a rotated 2D flat-topped (super-Gaussian) function.

    Parameters:
    - x_values: 1D array representing the x-axis values
    - y_values: 1D array representing the y-axis values
    - sigma_x: Standard deviation in the x direction
    - sigma_y: Standard deviation in the y direction
    - phi: Shape parameter (higher = flatter top)
    - theta_rad: Rotation angle in radians

    Returns:
    - 2D rotated super-Gaussian values
    """

    # Create a meshgrid
    X, Y = np.meshgrid(x_values, y_values)

    # Center coordinates (optional, assume mean-centered)
    X_shifted, Y_shifted = X - mean_x, Y - mean_y

    # Create the rotation matrix
    rotation_matrix = np.array([[np.cos(theta_rad), -np.sin(theta_rad)],
                                [np.sin(theta_rad), np.cos(theta_rad)]])

    # Rotate the coordinates
    coords = np.vstack([X_shifted.ravel(), Y_shifted.ravel()])
    rotated_coords = rotation_matrix @ coords
    X_rot, Y_rot = rotated_coords.reshape(2, *X.shape)  # Reshape back to 2D
    
    # Compute the super-Gaussian function using absolute values and phi exponent
    I_0 = 1 / (sigma_x * sigma_y * np.sqrt(2 * np.pi))
    
    gaussian = I_0 * np.exp(-0.5 * (np.abs(X_rot / sigma_x) ** phi + np.abs(Y_rot / sigma_y) ** phi))
    
    # Normalize the Gaussian to a maximum of 1
    gaussian /= np.max(gaussian)

    return gaussian, X_rot, Y_rot
# ----------------------------------------------------------------------------------------



def gaussian_2d_flat_topped_tilted_model(nx, ny, theta_rad, phi, BMAJ_pix, BMIN_pix, RA_centre_pix, Dec_centre_pix):
    """
    Computes a 2D Gaussian model with uniform and azimuthal ratio components.
    
    Parameters:
    X, Y : ndarray
        Meshgrid coordinate arrays.
    theta_rad : float
        Rotation angle in radians.
    std_x, std_y : float
        Standard deviations along x and y axes.
    mean_x, mean_y : float
        Mean values along x and y axes.
    flat_end_x, flat_end_y : float
        Boundaries for the flat region in x and y directions.

    Returns:
    GaussianUniformRatios : ndarray
        2D array representing the uniform Gaussian component.
    GaussianAzimuthalRatios : ndarray
        2D array representing the azimuthal Gaussian component.
    """
    
    # Generate x values
    x_values = np.linspace(0, nx, nx)

    # Generate y values
    y_values = np.linspace(0, ny, nx)
    
    # Calculate the sigma of minor and major axes
    # Get FWHM from BMAJ_pix and BMIN_pix
    beam_fwhm_major = BMAJ_pix 
    beam_fwhm_minor = BMIN_pix

    # Convert FWHM to standard deviation (σ = FWHM / (2√(2ln2)))
    beam_sigma_major = beam_fwhm_major / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    beam_sigma_minor = beam_fwhm_minor / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    
    # print(f'Sigma along the major axis is: {beam_sigma_major:.3f} pixels')
    # print(f'Sigma along the minor axis is: {beam_sigma_minor:.3f} pixels')



    GaussianUniformRatios, X_rot, Y_rot = calculate_gaussian_2d_flat_tilted(x_values, y_values, 
                                                                            beam_sigma_major, beam_sigma_minor, 
                                                                            RA_centre_pix, Dec_centre_pix,
                                                                            phi, theta_rad)
    
    GaussianAzimuthalRatios = 1 - GaussianUniformRatios
    
    # Ensure the sum of ratios is 1 everywhere
    assert np.allclose(GaussianUniformRatios + GaussianAzimuthalRatios, 1), "Error: The sum of ratios is not equal to 1."
    
    return GaussianUniformRatios, GaussianAzimuthalRatios, X_rot, Y_rot








# def mix_StokesQU_and_generate_vectors_gaussian(GaussianUniform_ratio, GaussianAzimuthal_ratio, 
#                                                StokesQ_uniform, StokesU_uniform, 
#                                                StokesQ_azimuthal, StokesU_azimuthal,
#                                                ny, nx, 
#                                                step, vector_length_pix_const, 
#                                                StokesI_data_2d_mJy, StokesIerr_data_2d_mJy,
#                                                calculated_polarized_intensity, PolarizedIntensity_err_data_2d_mJy,
#                                                PolarizationAngle_err_data_2d_deg):
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
# #     # Ensure the ratios sum to 1 (for safety)
# #     if not np.isclose(GaussianUniform_ratio + GaussianAzimuthal_ratio, 1.0):
# #         raise ValueError("Uniform_ratio and Azimuthal_ratio must sum to 1.")

#     # Compute weighted sum of Stokes Q and U
#     StokesQ_grid_mixed = GaussianUniform_ratio * StokesQ_uniform + GaussianAzimuthal_ratio * StokesQ_azimuthal
#     StokesU_grid_mixed = GaussianUniform_ratio * StokesU_uniform + GaussianAzimuthal_ratio * StokesU_azimuthal

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
#                 StokesI_data_2d_mJy[y, x] / StokesIerr_data_2d_mJy[y, x] > 3
#                 and calculated_polarized_intensity[y, x] / PolarizedIntensity_err_data_2d_mJy[y, x] > 3
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



