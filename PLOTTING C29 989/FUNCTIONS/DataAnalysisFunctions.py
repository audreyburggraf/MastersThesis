import numpy as np 

# https://radio-beam.readthedocs.io/en/latest/api/radio_beam.EllipticalGaussian2DKernel.html#radio_beam.EllipticalGaussian2DKernel

#https://astro-docs.readthedocs.io/en/latest/api/astropy.convolution.Gaussian2DKernel.html

def smooth_stokes_maps(Q_map, U_map, BMAJ_pix, BMIN_pix, BPA_rad_cartesian):
    """
    Apply Gaussian smoothing to Stokes Q and U maps based on the beam size.

    Parameters:
    - Q_map (numpy.ndarray): Stokes Q map.
    - U_map (numpy.ndarray): Stokes U map.
    - beam_fwhm_major (float): Beam major axis in pixels.
    - beam_fwhm_minor (float): Beam minor axis in pixels.

    Returns:
    - Q_smooth (numpy.ndarray): Smoothed Stokes Q map.
    - U_smooth (numpy.ndarray): Smoothed Stokes U map.
    """
    
    # Get FWHM from BMAJ_pix and BMIN_pix
    beam_fwhm_major = BMAJ_pix 
    beam_fwhm_minor = BMIN_pix
    
    # Convert FWHM to standard deviation (σ = FWHM / (2√(2ln2)))
    beam_sigma_major = beam_fwhm_major / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    beam_sigma_minor = beam_fwhm_minor / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    
    # This is my older kernel that just takes the average
    # |
    # Use the average beam sigma for a symmetric Gaussian kernel
    # sigma_average = np.mean([beam_sigma_major, beam_sigma_minor])
    # psf = Gaussian2DKernel(sigma_average)
    # |
    
    # Now I will try to use an elliptical kernel
    # Gaussian2DKernel(x_std, y_std, theta, x_size, y_size)
    kernel = Gaussian2DKernel(beam_sigma_major, beam_sigma_minor, -BPA_rad_cartesian)

    # Apply Gaussian smoothing
    Q_smooth = convolve(Q_map, kernel)
    U_smooth = convolve(U_map, kernel)

    return Q_smooth, U_smooth



# def calculate_chi_squared(observed_value, expected_value):
#     # Compute the Chi-squared value for a single observation (angle differences)
#     # Calculate the difference between observed and expected values
#     chi_squared = np.sum((observed_value - expected_value) ** 2 / expected_value)
    
#     return chi_squared




def normalize_angle(angle_rad):
    
    if angle_rad < 0:
        
        angle_rad = angle_rad + np.pi
        
    return angle_rad



def calculate_chi_squared_v2(observed_value, expected_value, tolerance=1e-6):
    
    # Normalize both observed and expected angles
    observed_value = np.array([normalize_angle(angle) for angle in observed_value])
    expected_value = np.array([normalize_angle(angle) for angle in expected_value])
    
    # Calculate the squared differences between normalized observed and expected values
    chi_squared = np.sum((observed_value - expected_value) ** 2)
    
    
    return chi_squared