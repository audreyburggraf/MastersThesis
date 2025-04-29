import numpy as np
from scipy.stats import norm
from collections import defaultdict



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





# ----------------------------------------------------------------------------------------
def run_gaussian_model_band4(theta_rad, phi_values, BMAJ_values_pix, BMIN_values_pix, RA_centre_pix, Dec_centre_pix, 
                             StokesQ_grid_100Uniform, StokesU_grid_100Uniform,
                             StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal,
                             vector_angle_actual_sky,
                             ny, nx, 
                             POLI_mJy, POLI_err_mJy, 
                             PA_err_deg,
                             print_statements = True):
    
    # Dictionary to store results
    results = {}
    values = []

    for phi_val in phi_values: 
        for BMAJ_val in BMAJ_values_pix:
            for BMIN_val in BMIN_values_pix:

                # Run the gaussian model function
                # -------------------------------------------------------------------------------------------------------------
                GaussianUniformRatios, GaussianAzimuthalRatios, _, _ = gaussian_2d_flat_topped_tilted_model(nx, ny, 
                                                                                                            theta_rad, phi_val, 
                                                                                                            BMIN_val,
                                                                                                            BMAJ_val, 
                                                                                                            RA_centre_pix,
                                                                                                            Dec_centre_pix)
                # -------------------------------------------------------------------------------------------------------------


                # Recover the Q, U and vector angle
                # -------------------------------------------------------------------------------------------------------------
                _, _, _, vectors_data, vectors_angle = mix_StokesQU_and_generate_vectors_band4(GaussianUniformRatios,
                                                                                               GaussianAzimuthalRatios, 
                                                                                               StokesQ_grid_100Uniform, 
                                                                                               StokesU_grid_100Uniform,
                                                                                               StokesQ_grid_100Azimuthal, 
                                                                                               StokesU_grid_100Azimuthal,
                                                                                               ny, nx, 
                                                                                               POLI_mJy, POLI_err_mJy, 
                                                                                               PA_err_deg)
                # -------------------------------------------------------------------------------------------------------------
                # Create a key for the dictionary based on the values
                value_str = f"{int(phi_val)}_{int(BMAJ_val)}_{int(BMIN_val)}"  

                # Save the results in the dictionary
                results[f"vectors_data_{value_str}"] = vectors_data
                # -------------------------------------------------------------------------------------------------------------



                # Calculate and append chi squared 
                # --------------------------------------------------------------
                chi_squared = calculate_chi_squared_v2(vectors_angle, vector_angle_actual_sky) # (observed, expected)
                # --------------------------------------------------------------
                
                values.append((phi_val, BMAJ_val, BMIN_val, chi_squared))
                
    
    
    # Find the index of the minimum chi-squared value
    chi_values = [entry[3] for entry in values]
    min_index = chi_values.index(min(chi_values))

    # Extract the best values
    best_phi, best_BMAJ, best_BMIN, best_chi_squared = values[min_index]

    # Print the results
    if print_statements:
        print(f'The lowest chi-squared value is: χ² = {best_chi_squared:.3f} when')
        print(f'    phi  = {best_phi:.2f}')
        print(f'    BMAJ = {best_BMAJ:.2f}')
        print(f'    BMIN = {best_BMIN:.2f}')
    
    # Access the best vector data from the dictionary
    best_key = f"vectors_data_{int(best_phi)}_{int(best_BMAJ)}_{int(best_BMIN)}"
    vector_data_gaussian_best = results[best_key]

    return values, results, vector_data_gaussian_best
    
    
# ----------------------------------------------------------------------------------------






# ----------------------------------------------------------------------------------------
def run_gaussian_model_band6(theta_rad, phi_values, BMAJ_values_pix, BMIN_values_pix, RA_centre_pix, Dec_centre_pix, 
                             StokesQ_grid_100Uniform, StokesU_grid_100Uniform,
                             StokesQ_grid_100Azimuthal, StokesU_grid_100Azimuthal,
                             vector_angle_actual_sky,
                             ny, nx, 
                             StokesI_mJy, StokesI_err_mJy,
                             POLI_mJy, POLI_err_mJy, 
                             PA_err_deg,
                             print_statements = True):
    
    # Dictionary to store results
    results = {}
    values = []

    for phi_val in phi_values: 
        for BMAJ_val in BMAJ_values_pix:
            for BMIN_val in BMIN_values_pix:

                # Run the gaussian model function
                # -------------------------------------------------------------------------------------------------------------
                GaussianUniformRatios, GaussianAzimuthalRatios, _, _ = gaussian_2d_flat_topped_tilted_model(nx, ny, 
                                                                                                            theta_rad, phi_val, 
                                                                                                            BMIN_val,
                                                                                                            BMAJ_val, 
                                                                                                            RA_centre_pix,
                                                                                                            Dec_centre_pix)
                # -------------------------------------------------------------------------------------------------------------


                # Recover the Q, U and vector angle
                # -------------------------------------------------------------------------------------------------------------
                _, _, _, vectors_data, vectors_angle = mix_StokesQU_and_generate_vectors_band6(GaussianUniformRatios,
                                                                                               GaussianAzimuthalRatios, 
                                                                                               StokesQ_grid_100Uniform, 
                                                                                               StokesU_grid_100Uniform,
                                                                                               StokesQ_grid_100Azimuthal, 
                                                                                               StokesU_grid_100Azimuthal,
                                                                                               ny, nx, 
                                                                                               StokesI_mJy, StokesI_err_mJy,
                                                                                               POLI_mJy, POLI_err_mJy, 
                                                                                               PA_err_deg)
                # -------------------------------------------------------------------------------------------------------------
                # Create a key for the dictionary based on the values
                value_str = f"{int(phi_val)}_{int(BMAJ_val)}_{int(BMIN_val)}"  

                # Save the results in the dictionary
                results[f"vectors_data_{value_str}"] = vectors_data
                # -------------------------------------------------------------------------------------------------------------



                # Calculate and append chi squared 
                # --------------------------------------------------------------
                chi_squared = calculate_chi_squared_v2(vectors_angle, vector_angle_actual_sky) # (observed, expected)
                # --------------------------------------------------------------
                
                values.append((phi_val, BMAJ_val, BMIN_val, chi_squared))
                
    
    
    # Find the index of the minimum chi-squared value
    chi_values = [entry[3] for entry in values]
    min_index = chi_values.index(min(chi_values))

    # Extract the best values
    best_phi, best_BMAJ, best_BMIN, best_chi_squared = values[min_index]

    # Print the results
    if print_statements:
        print(f'The lowest chi-squared value is: χ² = {best_chi_squared:.3f} when')
        print(f'    phi  = {best_phi:.2f}')
        print(f'    BMAJ = {best_BMAJ:.2f}')
        print(f'    BMIN = {best_BMIN:.2f}')
    
    # Access the best vector data from the dictionary
    best_key = f"vectors_data_{int(best_phi)}_{int(best_BMAJ)}_{int(best_BMIN)}"
    vector_data_gaussian_best = results[best_key]

    return values, results, vector_data_gaussian_best
    
    
# ----------------------------------------------------------------------------------------












# ----------------------------------------------------------------------------------------
def make_gaussian_grids(gaussian_values, phi_values, BMAJ_values_pix, BMIN_values_pix):
    """
    Parameters:
    -----------
    gaussian_values : list of tuples
        Each tuple should contain (phi, BMAJ, BMIN, chi^2).
    phi_values : array-like
        Unique values of phi used in the fit grid.
    BMAJ_values_pix : array-like
        Unique values of BMAJ used in the fit grid.
    BMIN_values_pix : array-like
        Unique values of BMIN used in the fit grid.
    """
    
    # Unpack the fit results into individual lists
    phi_list, BMAJ_list, BMIN_list, chi_list = zip(*gaussian_values)

    # Reshape the lists into 3D grids for later use (optional)
    phi_grid  = np.array(phi_list).reshape(len(phi_values),  len(BMAJ_values_pix), len(BMIN_values_pix))
    BMAJ_grid = np.array(BMAJ_list).reshape(len(phi_values), len(BMAJ_values_pix), len(BMIN_values_pix))
    BMIN_grid = np.array(BMIN_list).reshape(len(phi_values), len(BMAJ_values_pix), len(BMIN_values_pix))
    chi_grid  = np.array(chi_list).reshape(len(phi_values),  len(BMAJ_values_pix), len(BMIN_values_pix))
    
    return phi_grid, BMAJ_grid, BMIN_grid, chi_grid
# ----------------------------------------------------------------------------------------





# # ----------------------------------------------------------------------------------------
# def analyze_gaussian_averages_vs_chi2(gaussian_values, phi_values, BMAJ_pix_values, BMIN_pix_values):
#     """
#     Analyzes chi-squared values from Gaussian fitting results and plots their averages
#     as functions of phi, BMAJ, and BMIN.

#     Parameters:
#     -----------
#     gaussian_values : list of tuples
#         Each tuple should contain (phi, BMAJ, BMIN, chi^2).
#     phi_values : array-like
#         Unique values of phi used in the fit grid.
#     BMAJ_pix_values : array-like
#         Unique values of BMAJ used in the fit grid.
#     BMIN_pix_values : array-like
#         Unique values of BMIN used in the fit grid.
#     """

#     phi_grid, BMAJ_grid, BMIN_grid, chi_grid = make_gaussian_grids(gaussian_values, phi_values, BMAJ_pix_values, BMIN_pix_values)

#     # Create dictionaries to group chi^2 values by phi, BMAJ, and BMIN
#     phi_dict  = defaultdict(list)
#     bmaj_dict = defaultdict(list)
#     bmin_dict = defaultdict(list)

#     # Populate the dictionaries
#     for phi, BMAJ, BMIN, chi in gaussian_values:
#         phi_dict[phi].append(chi)
#         bmaj_dict[BMAJ].append(chi)
#         bmin_dict[BMIN].append(chi)

#     # Compute average chi^2 values for each parameter value
#     phi_avg  = {phi: np.mean(chi_values) for phi, chi_values in phi_dict.items()}
#     bmaj_avg = {BMAJ: np.mean(chi_values) for BMAJ, chi_values in bmaj_dict.items()}
#     bmin_avg = {BMIN: np.mean(chi_values) for BMIN, chi_values in bmin_dict.items()}

#     # Plot the averaged results
#     def plot_gaussian_averages(phi_avg, bmaj_avg, bmin_avg):
#     # Create figure
#     fig, ax = plt.subplots(1, 3, figsize=(18, 5))

#     # Plot phi vs chi^2
#     ax[0].plot(list(phi_avg.keys()), list(phi_avg.values()), marker='o', linestyle='-', color = 'blue', lw = 4, ms = 15)
#     ax[0].set_xlabel(r'$\phi$', fontsize = axis_label_fs)
#     ax[0].set_title(r'$\phi$ vs $\chi^2$', fontsize = title_fs)

#     # Plot BMAJ vs chi^2
#     ax[1].plot(list(bmaj_avg.keys()), list(bmaj_avg.values()), marker='o', linestyle='-', color = 'red', lw = 4, ms = 15)
#     ax[1].set_xlabel('BMAJ', fontsize = axis_label_fs)
#     ax[1].set_title(r'BMAJ vs $\chi^2$', fontsize = title_fs)

#     # Plot BMIN vs chi^2
#     ax[2].plot(list(bmin_avg.keys()), list(bmin_avg.values()), marker='o', linestyle='-', color = 'forestgreen', lw = 4, ms = 15)
#     ax[2].set_xlabel('BMIN', fontsize = axis_label_fs)
#     ax[2].set_title(r'BMIN vs $\chi^2$', fontsize = title_fs)

#     for i in range (3):
#         # they all have the same y-axis label as chi^2
#         ax[i].set_ylabel(r'$\chi^2$', fontsize = axis_label_fs)

#         # Set minor and major ticks
#         ax[i].minorticks_on()
#         ax[i].tick_params(axis="x", which="major", direction="in", length=7, labelsize=axis_num_fs)
#         ax[i].tick_params(axis="y", which="major", direction="in", length=7, labelsize=axis_num_fs)

#     plt.tight_layout()
# # ----------------------------------------------------------------------------------------
