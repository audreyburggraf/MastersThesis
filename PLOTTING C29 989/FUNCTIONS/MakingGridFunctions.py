import numpy as np



# Import the constants
# -----------------------------------------------------------------------------------------
import sys

sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/") 

import constants

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
def compute_polarization_vector(x, y, PA_grid, band, vector_len_pix = None):

    """
    Compute the vector components for polarization at the given (x, y) position.

    Parameters:
    x, y: pixel indices
    PA_grid: 2D array of polarization angles in radians

    Returns:
    A list containing [x_start, x_end, y_start, y_end] for the polarization vector,
    and the polarization angle in radians.
    """
    
    vector_length_defaults = {
        'Band 4': constants.vector_len_pix_band4,
        'Band 4 nterms2': constants.vector_len_pix_band4_nterms2,
        'Band 4 nterms2 robust -1': constants.vector_len_pix_band4_nterms2_robust_minus1,
        'Band 4 nterms2 smooth': constants.vector_len_pix_band4_nterms2_smooth,
        'Band 4 nterms2 smooth B6': constants.vector_len_pix_band4_nterms2_smooth_B6,
        'Band 4 nterms2 smooth B6 B7': constants.vector_len_pix_band4_nterms2_smooth_B6_B7,
        
        'Band 5': constants.vector_len_pix_band5,
        'Band 5 v0': constants.vector_len_pix_band5_v0,
        'Band 5 robust -2': constants.vector_len_pix_band5_robust_minus2,
        'Band 5 robust -1': constants.vector_len_pix_band5_robust_minus2,
        
        'Band 6': constants.vector_len_pix_band6,
        'Band 6 smooth': constants.vector_len_pix_band6_smooth,
        'Band 6 smooth B7': constants.vector_len_pix_band6_smooth_B7,
        
        'Band 7 nterms2': constants.vector_len_pix_band7_nterms2,
        'Band 7 nterms2 smooth': constants.vector_len_pix_band7_nterms2_smooth,
        'Band 7 nterms2 smooth B6': constants.vector_len_pix_band7_nterms2_smooth_B6,
    }

    if vector_len_pix is None:
        if band not in vector_length_defaults:
            raise ValueError(f"Unsupported band: {band}")

        vector_len_pix = vector_length_defaults[band]

    
    
    # Extract the polarization angle at this location
    PA_rad_sky = PA_grid[y, x] 
    
#     print(rf'PA_rad_sky = {PA_rad_sky}')

    # Compute vector components
    dx = vector_len_pix * np.cos(PA_rad_sky + np.pi/2)
    dy = vector_len_pix * np.sin(PA_rad_sky + np.pi/2)
    
    # Vector in Cartesian coordinates
    vector_cartesian = [x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2]
    
    return vector_cartesian, PA_rad_sky
# --------------------------------------------------------------------------
def make_vectors(ny, nx, 
                 xmin, xmax, ymin, ymax, 
                 POLI_mJy, POLI_err_mJy, PA_grid, PA_err_deg, band, 
                 step = None, vector_len_pix = None):
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
    
    band_parameters = {
        'Band 4': (constants.step_band4, 4),
        'Band 4 nterms2': (constants.step_band4_nterms2, 4),
        'Band 4 nterms2 robust -1': (constants.step_band4_nterms2_robust_minus1, 4),
        'Band 4 nterms2 smooth': (constants.step_band4_nterms2_smooth, 4),
        'Band 4 nterms2 smooth B6': (constants.step_band4_nterms2_smooth_B6, 4),
        'Band 4 nterms2 smooth B6 B7': (constants.step_band4_nterms2_smooth_B6_B7, 4),

        'Band 5': (constants.step_band5, 4),
        'Band 5 v0': (constants.step_band5_v0, 4),
        'Band 5 robust -2': (constants.step_band5_robust_minus2, 4),
        'Band 5 robust -1': (constants.step_band5_robust_minus1, 4),

        'Band 6': (constants.step_band6, 3),
        'Band 6 smooth': (constants.step_band6_smooth, 3),
        'Band 6 smooth B7': (constants.step_band6_smooth_B7, 3),

        'Band 7 nterms2': (constants.step_band7_nterms2, 4),
        'Band 7 nterms2 smooth': (constants.step_band7_nterms2_smooth, 4),
        'Band 7 nterms2 smooth B6': (constants.step_band7_nterms2_smooth_B6, 4),
    }

    if band not in band_parameters:
        raise ValueError(f"Unsupported band: {band}")

    default_step, POLI_err_cutoff = band_parameters[band]

    # If step wasn't passed, use the band default
    if step is None:
        step = default_step
        
    # Make some empty grid that we will fill with the mask on if there is a vector or not
    vector_mask = np.zeros((ny, nx), dtype=bool)
    in_plot_mask = np.zeros((ny, nx), dtype=bool)
    
    # Empty arrays
    vectors_cartesian = []
    vector_angles_sky = []
    
    for x in range(0, nx, step):
        for y in range(0, ny, step):
            if (POLI_mJy[y, x] / POLI_err_mJy[y, x] > POLI_err_cutoff
                and PA_err_deg[y, x] < 10):
                # Use the helper function to compute the vector
                vector_cartesian, PA_rad_sky = compute_polarization_vector(x, y, PA_grid, band, vector_len_pix)
                vectors_cartesian.append(vector_cartesian)
                vector_angles_sky.append(PA_rad_sky)
                
                # Save the mask
                vector_mask[y, x] = True
                
                # Check that the vector is in the plotting range 
                in_plot = (xmin <= x <= xmax) and (ymin <= y <= ymax)

                if in_plot:
                    in_plot_mask[y, x] = True
    
    
    return vectors_cartesian, vector_angles_sky, vector_mask, in_plot_mask
# --------------------------------------------------------------------------



# def make_vectors_band47(ny, nx, POLI_mJy, POLI_err_mJy, PA_grid, PA_err_deg, band, 
#                         step = None, vector_len_pix = None):
#     """
#     Generate vectors for Band 4 polarization data.
    
#     Parameters:
#     ny, nx: Dimensions of the grid
#     POLI_mJy: Polarization intensity
#     POLI_err_mJy: Error on polarization intensity
#     PA_grid: Polarization angle grid
#     PA_err_deg: Polarization angle error
    
#     Returns:
#     vectors_cartesian: List of vectors in Cartesian coordinates
#     vector_angles_sky: List of polarization angles
#     """
    
#     if step is None:
#         if band == 'band 4':
#             step = constants.step_band4
#         elif band == 'band 5':
#             step = constants.step_band5
#         elif band == 'band 7':
#             step = constants.step_band7
#         else:
#             return('Currently only accepting Band 4 and Band 7')
        
#     vectors_cartesian = []
#     vector_angles_sky = []
    
#     for x in range(0, nx, step):
#         for y in range(0, ny, step):
#             if (POLI_mJy[y, x] / POLI_err_mJy[y, x] > 4
#                 and PA_err_deg[y, x] < 10):
#                 # Use the helper function to compute the vector
#                 vector_cartesian, PA_rad_sky = compute_polarization_vector(x, y, PA_grid, band, vector_len_pix)
#                 vectors_cartesian.append(vector_cartesian)
#                 vector_angles_sky.append(PA_rad_sky)
    
#     return vectors_cartesian, vector_angles_sky

# # --------------------------------------------------------------------------
# def make_vectors_band6(ny, nx, 
#                        StokesI_mJy, StokesI_err_mJy, 
#                        POLI_mJy, POLI_err_mJy, 
#                        PA_grid, PA_err_deg,
#                        step = None):
#     """
#     Generate vectors for Band 6 polarization data.
    
#     Parameters:
#     ny, nx: Dimensions of the grid
#     StokesI_mJy: Stokes I intensity
#     StokesI_err_mJy: Error on Stokes I intensity
#     POLI_mJy: Polarization intensity
#     POLI_err_mJy: Error on polarization intensity
#     PA_grid: Polarization angle grid
#     PA_err_deg: Polarization angle error
    
#     Returns:
#     vectors_cartesian: List of vectors in Cartesian coordinates
#     vector_angles_sky: List of polarization angles
#     """
    
#     if step is None:
#         step = constants.step_band6
        
        
#     vectors_cartesian = []
#     vector_angles_sky = []
    
#     for x in range(0, nx, step):
#         for y in range(0, ny, step):
#             if (StokesI_mJy[y, x] / StokesI_err_mJy[y, x] > 3 and 
#                 POLI_mJy[y, x] / POLI_err_mJy[y, x] > 3 and 
#                 PA_err_deg[y, x] < 10):
#                 # Use the helper function to compute the vector
#                 vector_cartesian, PA_rad_sky = compute_polarization_vector(x, y, PA_grid, band = 'band 6')
#                 vectors_cartesian.append(vector_cartesian)
#                 vector_angles_sky.append(PA_rad_sky)
    
#     return vectors_cartesian, vector_angles_sky
# # --------------------------------------------------------------------------