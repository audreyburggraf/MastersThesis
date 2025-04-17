import numpy as np



# Import the constants
# -----------------------------------------------------------------------------------------
import sys

sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/") 

import constants

step_band4 = constants.step_band4
vector_length_pix_const_band4 = constants.vector_length_pix_const_band4
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


def make_vectors_band4(ny, nx,  
                       POLI_mJy, POLI_err_mJy,
                       PA_grid, PA_err_data_2d_deg):
    
    vectors_cartesian_band4 = []
    vector_angles_sky_band4 = []
    
    for x in range(0, nx, step_band4):
        for y in range(0, ny, step_band4):
            if (POLI_mJy[y, x] / POLI_err_mJy[y, x] > 3
                and PA_err_data_2d_deg[y, x] < 10
            ):
                # Extract the polarization angle at this location
                PA_rad_sky = PA_grid[y, x] 

                # Compute vector components
                dx = vector_length_pix_const_band4 * np.cos(PA_rad_sky + np.pi/2)
                dy = vector_length_pix_const_band4 * np.sin(PA_rad_sky + np.pi/2)

                # Append to the list in the format [x_start, x_end, y_start, y_end]
                vectors_cartesian_band4.append([x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2])

                vector_angles_sky_band4.append(PA_rad_sky)
    
    return vectors_cartesian_band4, vector_angles_sky_band4