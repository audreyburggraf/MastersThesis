import numpy as np



# Import the constants
# -----------------------------------------------------------------------------------------
import sys

sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/") 

import constants

step = constants.step
vector_length_pix_const = constants.vector_length_pix_const
# -----------------------------------------------------------------------------------------


def make_PA_grid_100Uniform(ny, nx, uniform_angle_rad_sky):
    """Generate a grid filled with a uniform position angle (in radians)."""
    
    PA_grid_100Uniform = np.full((ny, nx), uniform_angle_rad_sky)
    
    return PA_grid_100Uniform


def make_PA_grid_100Azimuthal(ny, nx, RA_centre_pix, Dec_centre_pix):
    """Generate a grid of azimuthal angles (in radians) centered at a given pixel."""
    PA_grid_100Azimuthal = np.zeros((ny, nx))

    for x in range(nx):
        for y in range(ny):
            dx_center = x - RA_centre_pix
            dy_center = y - Dec_centre_pix

            # Calculate azimuthal angle in sky coordinates (tangent to radial)
            azimuthal_angle_rad_astronomy = np.arctan2(dy_center, dx_center)

            PA_grid_100Azimuthal[y, x] = azimuthal_angle_rad_astronomy

    return PA_grid_100Azimuthal



# --------------------------------------------------------------------------
def compute_polarization_vector(x, y, PA_grid):
    """
    Compute the vector components for polarization at the given (x, y) position.

    Parameters:
    x, y: pixel indices
    PA_grid: 2D array of polarization angles in radians

    Returns:
    A list containing [x_start, x_end, y_start, y_end] for the polarization vector,
    and the polarization angle in radians.
    """
    # Extract the polarization angle at this location
    PA_rad_sky = PA_grid[y, x] 

    # Compute vector components
    dx = vector_length_pix_const * np.cos(PA_rad_sky + np.pi/2)
    dy = vector_length_pix_const * np.sin(PA_rad_sky + np.pi/2)
    
    # Vector in Cartesian coordinates
    vector_cartesian = [x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2]
    
    return vector_cartesian, PA_rad_sky
# --------------------------------------------------------------------------

def make_vectors_band4(ny, nx, POLI_mJy, POLI_err_mJy, PA_grid, PA_err_deg):
    """
    Generate vectors for Band 4 polarization data.
    
    Parameters:
    ny, nx: Dimensions of the grid
    POLI_mJy: Polarization intensity
    POLI_err_mJy: Error on polarization intensity
    PA_grid: Polarization angle grid
    PA_err_deg: Polarization angle error
    
    Returns:
    vectors_cartesian: List of vectors in Cartesian coordinates
    vector_angles_sky: List of polarization angles
    """
    vectors_cartesian = []
    vector_angles_sky = []
    
    for x in range(0, nx, step):
        for y in range(0, ny, step):
            if (POLI_mJy[y, x] / POLI_err_mJy[y, x] > 3.3
                and PA_err_deg[y, x] < 10):
                # Use the helper function to compute the vector
                vector_cartesian, PA_rad_sky = compute_polarization_vector(x, y, PA_grid)
                vectors_cartesian.append(vector_cartesian)
                vector_angles_sky.append(PA_rad_sky)
    
    return vectors_cartesian, vector_angles_sky
# --------------------------------------------------------------------------
def make_vectors_band6(ny, nx, 
                       StokesI_mJy, StokesI_err_mJy, 
                       POLI_mJy, POLI_err_mJy, 
                       PA_grid, PA_err_deg):
    """
    Generate vectors for Band 6 polarization data.
    
    Parameters:
    ny, nx: Dimensions of the grid
    StokesI_mJy: Stokes I intensity
    StokesI_err_mJy: Error on Stokes I intensity
    POLI_mJy: Polarization intensity
    POLI_err_mJy: Error on polarization intensity
    PA_grid: Polarization angle grid
    PA_err_deg: Polarization angle error
    
    Returns:
    vectors_cartesian: List of vectors in Cartesian coordinates
    vector_angles_sky: List of polarization angles
    """
    vectors_cartesian = []
    vector_angles_sky = []
    
    for x in range(0, nx, step):
        for y in range(0, ny, step):
            if (StokesI_mJy[y, x] / StokesI_err_mJy[y, x] > 3 and 
                POLI_mJy[y, x] / POLI_err_mJy[y, x] > 3 and 
                PA_err_deg[y, x] < 10):
                # Use the helper function to compute the vector
                vector_cartesian, PA_rad_sky = compute_polarization_vector(x, y, PA_grid)
                vectors_cartesian.append(vector_cartesian)
                vector_angles_sky.append(PA_rad_sky)
    
    return vectors_cartesian, vector_angles_sky
# --------------------------------------------------------------------------