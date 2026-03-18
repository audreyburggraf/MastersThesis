# Import Functions
from FITS_Image_Functions import *
from DataAnalysisFunctions import *
from scipy.stats import norm

# Import the constants
# -----------------------------------------------------------------------------------------
import sys

sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/") 

import constants
# -----------------------------------------------------------------------------------------

import pandas as pd

def extract_axis_data(axis_x, axis_y, data_2d, minor_or_major, centre_pix, gridsize, header, print_statement = False):
    """
    Extracts data along a specified axis (major or minor) and returns the data along with valid positions.

    Parameters:
    - axis_x: Array of x positions along the axis (e.g., major_x or minor_x)
    - axis_y: Array of y positions along the axis (e.g., major_y or minor_y)
    - data_2d: 2D array representing the data (e.g., StokesI_data_2d_mJy)
    - centre_pix: Tuple (RA_centre_pix, Dec_centre_pix) specifying the center pixel position
    - gridsize: Tuple representing the dimensions of the grid (data_2d.shape)
    - header: FITS header to extract the plate scale for arcsecond conversion

    Returns:
    - axis_data: List of extracted data values
    - offset_pixels: List of pixel offsets (positive and negative)
    - offset_arcsec: List of arcsecond offsets (positive and negative)
    """
   
    axis_data = []
    offset_pixels = []
    offset_arcsec = []

    RA_centre_pix, Dec_centre_pix = centre_pix  # Unpack center coordinates

    # Extract pixel scale from header (arcsec per pixel)
    if 'CDELT1' in header:
        pixel_scale_arcsec = abs(header['CDELT1']) * 3600  # Convert degrees to arcsec
    elif 'CD1_1' in header:
        pixel_scale_arcsec = abs(header['CD1_1']) * 3600  # Alternative plate scale
    else:
        raise ValueError("Could not determine plate scale from FITS header.")

    for i in range(len(axis_x)):
        # Round to nearest integer for pixel indexing
        xi, yi = int(round(axis_x[i])), int(round(axis_y[i]))

        # Ensure indices are within bounds
        if 0 <= xi < gridsize[1] and 0 <= yi < gridsize[0]:
            # Extract data
            axis_data.append(data_2d[yi, xi])

            # Compute pixel offsets
            delta_x = xi - RA_centre_pix
            delta_y = yi - Dec_centre_pix
            offset_pixel = np.sqrt(delta_x**2 + delta_y**2)
            
            # Determine offset based on 'major' or 'minor'
            if minor_or_major == 'major':
                # If both delta_x and delta_y are negative, the offset is negative
                if delta_x < 0: #and delta_y < 0:
                    offset_pixels.append(offset_pixel)  # Negative offset in pixels
                    offset_arcsec.append(offset_pixel * pixel_scale_arcsec)  # Negative offset in arcseconds
                else:
                    offset_pixels.append(-offset_pixel)  # Positive offset in pixels
                    offset_arcsec.append(-offset_pixel * pixel_scale_arcsec)  # Positive offset in arcseconds

            elif minor_or_major == 'minor':
                # If delta_x is negative, the offset is negative
                if delta_x < 0:
                    offset_pixels.append(offset_pixel)  # Negative offset in pixels
                    offset_arcsec.append(offset_pixel * pixel_scale_arcsec)  # Negative offset in arcseconds
                else:
                    offset_pixels.append(-offset_pixel)  # Positive offset in pixels
                    offset_arcsec.append(-offset_pixel * pixel_scale_arcsec)  # Positive offset in arcseconds

    if print_statement: 
        print(' ')
        print('in extract_axis_data')
        print('gridsize:', gridsize)
        print('length of axis_data:', len(axis_data))
        print('length of offset_pixels:', len(offset_pixels))
        print('length of offset_arcsec:', len(offset_arcsec))
        print(' ')

    return axis_data, offset_pixels, offset_arcsec




# -----------------------------------------------------------------------------------------
def fit_slices_slope(beam_area_arcsec2, intensity_mJy_beam, offsets_arcsec): 
    
    intensity_mJy_arcsec2 = np.array(intensity_mJy_beam) /beam_area_arcsec2
    
    # Convert inputs to numpy arrays
    # ------------------------------------------------------
    offsets_arcsec = np.array(offsets_arcsec)
    intensity_mJy_arcsec2 = np.array(intensity_mJy_arcsec2)
    # ------------------------------------------------------
    
    # Mask to separate positive and negative side
    mask_pos = (intensity_mJy_arcsec2 > 0) & (offsets_arcsec > 0)  # Positive side
    mask_neg = (intensity_mJy_arcsec2 > 0) & (offsets_arcsec < 0)  # Negative side


    # Log-log fit
    log_r_pos = np.log10(np.abs(offsets_arcsec[mask_pos]))
    log_I_pos = np.log10(intensity_mJy_arcsec2[mask_pos])
    
    log_r_neg = np.log10(np.abs(offsets_arcsec[mask_neg]))
    log_I_neg = np.log10(intensity_mJy_arcsec2[mask_neg])
    
    # Linear fit in log-log space
    slope_pos, intercept_pos = np.polyfit(log_r_pos, log_I_pos, 1)
    slope_neg, intercept_neg = np.polyfit(log_r_neg, log_I_neg, 1)
    
    return -slope_pos, intercept_pos, -slope_neg, intercept_neg
# -----------------------------------------------------------------------------------------




# -----------------------------------------------------------------------------------------

def run_slices(data, StokesI_header, StokesI_wcs, carta_minor_data, carta_major_data,
               carta_minor_offset, carta_major_offset, band, print_statement = False):
    """
    Generate and compare slice data along the major and minor axes of an image
    for a given ALMA band, optionally computing chi-squared values between the 
    slices and reference CARTA data.

    Parameters
    ----------
    data : 2D numpy array
        Image data (e.g., Stokes I).
    StokesI_header : FITS header
        Header from the Stokes I FITS image.
    StokesI_wcs : astropy.wcs.WCS
        World Coordinate System for the image.
    carta_minor_data : 1D array-like
        CARTA-provided intensity values along the minor axis.
    carta_major_data : 1D array-like
        CARTA-provided intensity values along the major axis.
    xmin, xmax, ymin, ymax : int
        Bounds of the image (in pixels).
    band : int
        ALMA band number. Currently supports band 4 and 6.

    Returns
    -------
    major_data : 1D numpy array
        Extracted pixel values along the major axis.
    major_offset_arcsec : 1D numpy array
        Offset positions along the major axis (in arcseconds).
    minor_data : 1D numpy array
        Extracted pixel values along the minor axis.
    minor_offset_arcsec : 1D numpy array
        Offset positions along the minor axis (in arcseconds).
    """

    band_parameters = {
        
       'Band 4': (
            constants.centre_str_band4,
            constants.major_angle_rad_cartesian_band4,
            constants.minor_angle_rad_cartesian_band4,
            1.4
        ),

        'Band 4 nterms2': (
            constants.centre_str_band4_nterms2,
            constants.major_angle_rad_cartesian_band4_nterms2,
            constants.minor_angle_rad_cartesian_band4_nterms2,
            1.4
        ),

        'Band 4 nterms2 smooth': (
            constants.centre_str_band4_nterms2_smooth,
            constants.major_angle_rad_cartesian_band4_nterms2_smooth,
            constants.minor_angle_rad_cartesian_band4_nterms2_smooth,
            1.4
        ),
        
       'Band 4 nterms2 smooth B6': (
            constants.centre_str_band4_nterms2_smooth_B6,
            constants.major_angle_rad_cartesian_band4_nterms2_smooth_B6,
            constants.minor_angle_rad_cartesian_band4_nterms2_smooth_B6,
            1.4
        ),
        
       'Band 4 nterms2 smooth B6 B7': (
            constants.centre_str_band4_nterms2_smooth_B6_B7,
            constants.major_angle_rad_cartesian_band4_nterms2_smooth_B6_B7,
            constants.minor_angle_rad_cartesian_band4_nterms2_smooth_B6_B7,
            1.4
        ),

        'Band 5 v0': (
            constants.centre_str_band5_v0,
            constants.major_angle_rad_cartesian_band5_v0,
            constants.minor_angle_rad_cartesian_band5_v0,
            1.8
        ),

        'Band 5': (
            constants.centre_str_band5,
            constants.major_angle_rad_cartesian_band5,
            constants.minor_angle_rad_cartesian_band5,
            1.8
        ),

        'Band 5 robust -2': (
            constants.centre_str_band5_robust_minus2,
            constants.major_angle_rad_cartesian_band5_robust_minus2,
            constants.minor_angle_rad_cartesian_band5_robust_minus2,
            1.8
        ),

        'Band 5 robust -1': (
            constants.centre_str_band5_robust_minus1,
            constants.major_angle_rad_cartesian_band5_robust_minus1,
            constants.minor_angle_rad_cartesian_band5_robust_minus1,
            1.8
        ),

        'Band 5 nterms2': (
            constants.centre_str_band5_nterms2,
            constants.major_angle_rad_cartesian_band5_nterms2,
            constants.minor_angle_rad_cartesian_band5_nterms2,
            1.8
        ),
        
        'Band 6': (
            constants.centre_str_band6,
            constants.major_angle_rad_cartesian_band6,
            constants.minor_angle_rad_cartesian_band6,
            1.4
        ),

        'Band 6 smooth': (
            constants.centre_str_band6_smooth,
            constants.major_angle_rad_cartesian_band6_smooth,
            constants.minor_angle_rad_cartesian_band6_smooth,
            1.4
        ),

        'Band 6 smooth B7': (
            constants.centre_str_band6_smooth_B7,
            constants.major_angle_rad_cartesian_band6_smooth_B7,
            constants.minor_angle_rad_cartesian_band6_smooth_B7,
            1.4
        ),

        'Band 7 nterms2': (
            constants.centre_str_band7_nterms2,
            constants.major_angle_rad_cartesian_band7_nterms2,
            constants.minor_angle_rad_cartesian_band7_nterms2,
            2.0
        ),

        'Band 7 nterms2 smooth': (
            constants.centre_str_band7_nterms2_smooth, 
            constants.major_angle_rad_cartesian_band7_nterms2_smooth,
            constants.minor_angle_rad_cartesian_band7_nterms2_smooth,
            2.0
        ),
        
       'Band 7 nterms2 smooth B6': (
            constants.centre_str_band7_nterms2_smooth_B6,  
            constants.major_angle_rad_cartesian_band7_nterms2_smooth_B6,
            constants.minor_angle_rad_cartesian_band7_nterms2_smooth_B6,
            2.0
        )
    }


    if band not in band_parameters:
        raise ValueError(f"Unsupported band: {band}")

    centre_str, major_angle_rad_cartesian, minor_angle_rad_cartesian, line_length_arcsec = band_parameters[band]


    # Normalize angles to standard range
    major_angle_rad_cartesian = normalize_angle(major_angle_rad_cartesian)
    minor_angle_rad_cartesian = normalize_angle(minor_angle_rad_cartesian)


    # Convert centre position from string to pixel coordinates
    centre_pix = list(string_to_pixel(centre_str, StokesI_wcs))


    # Convert desired line length from arcsec to pixels
    line_length_pix = arcsec_to_pixels(StokesI_header, line_length_arcsec)

    # Sampling points (based on CARTA slice length)
    num_points = len(carta_major_data)
    gridsize = data.shape

    # Define coordinates along the major axis
    delta = np.linspace(-line_length_pix / 2, line_length_pix / 2, num_points)
    major_x = centre_pix[0] + delta * np.cos(major_angle_rad_cartesian)
    major_y = centre_pix[1] + delta * np.sin(major_angle_rad_cartesian)

    # Define coordinates along the minor axis
    minor_x = centre_pix[0] + delta * np.cos(minor_angle_rad_cartesian)
    minor_y = centre_pix[1] + delta * np.sin(minor_angle_rad_cartesian)
    
    if print_statement:
        print(' ')
        print('in run_slices')
        print('The length of major_x, _y, minor_x, _ is:', len(major_x), len(major_y), len(minor_x), len(minor_y))

    # Extract image values along the defined axes
    major_data, _, major_offset_arcsec = extract_axis_data(major_x, major_y, 
                                                           data, 'major', 
                                                           centre_pix, gridsize, 
                                                           StokesI_header,
                                                           print_statement)
    
    
    minor_data, _, minor_offset_arcsec = extract_axis_data(minor_x, minor_y, 
                                                           data, 'minor', 
                                                           centre_pix, gridsize, 
                                                           StokesI_header,
                                                           print_statement)
    
    points = [major_x, major_y, minor_x, minor_y]
    
    if print_statement:
        print(rf"The major angle (cartesian) is: {major_angle_rad_cartesian * 180/np.pi:.1f} degrees")
        print(rf"The minor angle (cartesian) is: {minor_angle_rad_cartesian * 180/np.pi:.1f} degrees")
        print(" ")
        print(f'Major start: ({major_x[-1]:.3f}, {major_y[-1]:.3f}) pixels')
        print(f'Major end:   ({major_x[0]:.3f}, {major_y[0]:.3f}) pixels')
        print(f'Minor start: ({minor_x[-1]:.3f}, {minor_y[-1]:.3f}) pixels')
        print(f'Minor end:   ({minor_x[0]:.3f}, {minor_y[0]:.3f}) pixels')
        print(rf" ")
        print("RA_centre_pix = ", centre_pix[0])
        print("Dec_centre_pix = ", centre_pix[1])

    return major_data, major_offset_arcsec, minor_data, minor_offset_arcsec, points
# -----------------------------------------------------------------------------------------



# def plot_gaussian_beam(ax, band, offsets_arcsec, color = 'black', 
#                        normalize=True, scale_to_data=None):
#     """
#     Plots projected 1D Gaussian beam along a given slice angle.

#     Parameters
#     ----------
#     ax : matplotlib axis
#     band : str
#     offsets_arcsec : array
#         X-axis offsets (arcsec) of the slice
#     normalize : bool
#         If True, normalize beam to peak = 1
#     scale_to_data : float or None
#         If provided, scales beam peak to this value
#     """

#     # ----------------------------
#     # Select correct constants file
#     # ----------------------------
#     if band == 'Band 6' or 6:
#         file = constants.band6_data_folder_path + 'beam_BAND6.csv'
#         minor_angle_rad_cartesian = constants.minor_angle_rad_cartesian_band6

#     elif band == 'Band 5' or 5:
#         file = constants.band5_data_folder_path + 'beam_BAND5.csv'
#         minor_angle_rad_cartesian = constants.minor_angle_rad_cartesian_band5

#     elif band == 'Band 4' or 4:
#         file = constants.band4_data_folder_path + 'beam_BAND4.csv'
#         minor_angle_rad_cartesian = constants.minor_angle_rad_cartesian_band4

#     elif band == 'Band 4 nterms2' or 42:
#         file = constants.band4_nterms2_data_folder_path + 'beam_BAND4_nterms2.csv'
#         minor_angle_rad_cartesian = constants.minor_angle_rad_cartesian_band4_nterms2

#     elif band == 'Band 7 nterms2 or 72':
#         file = constants.band7_nterms2_data_folder_path + 'beam_BAND7_nterms2.csv'
#         minor_angle_rad_cartesian = constants.minor_angle_rad_cartesian_band7_nterms2

#     else:
#         raise ValueError("Unsupported band.")

#     df = pd.read_csv(file)

#     BMAJ_arcsec = df.at[0, "BMAJ_arcsec"]
#     BMIN_arcsec = df.at[0, "BMIN_arcsec"]
#     BPA_deg     = df.at[0, "BPA_deg"]

#     # We only care about the minor axis, and the minor FWHM is BMIN_arcsec
    
    
#     # FWHM to sigma 
#     # --------------------------------------------
#     # FWHM = 2 sqrt(2 ln 2) sigma 
#     # sigma = FWHM /2 sqrt(2 ln(2)
#     # --------------------------------------------
#     fwhm_to_sigma = 1 / (2 * np.sqrt(2 * np.log(2)))
    
#     sigma_maj = BMAJ_arcsec * fwhm_to_sigma
#     sigma_min = BMIN_arcsec * fwhm_to_sigma
#     # --------------------------------------------

#     # Convert BPA to radians
#     BPA_rad = np.deg2rad(BPA_deg)

#     # ----------------------------
#     # Effective sigma along slice
#     # ----------------------------
#     delta_theta = minor_angle_rad_cartesian  - BPA_rad

#     sigma_eff = np.sqrt(
#         (sigma_maj * np.cos(delta_theta))**2 +
#         (sigma_min * np.sin(delta_theta))**2
#     )

#     # Calculate the 1D gaussian using scipy.stats
#     beam_profile = norm.pdf(offsets_arcsec, loc=0, scale=sigma_eff)

#     if normalize:
#         beam_profile /= beam_profile.max()

#     if scale_to_data is not None:
#         beam_profile *= scale_to_data

#     # ----------------------------
#     # Plot
#     # ----------------------------
#     ax.plot(offsets_arcsec,
#             beam_profile,
#             linestyle='--',
#             color=color,
#             linewidth=2,
#             label='Projected Beam')

#     #return sigma_eff

    
    
    

