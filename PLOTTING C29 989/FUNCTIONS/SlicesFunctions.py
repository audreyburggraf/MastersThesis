# Import Functions
from FITS_Image_Functions import *
from DataAnalysisFunctions import *

# Import the constants
# -----------------------------------------------------------------------------------------
import sys

sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/") 

import constants
# -----------------------------------------------------------------------------------------



def extract_axis_data(axis_x, axis_y, data_2d, minor_or_major, centre_pix, gridsize, header):
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



    return axis_data, offset_pixels, offset_arcsec



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

    # Select band-specific parameters
    if band == 6:
        centre_str = constants.centre_str_band6
        major_angle_rad_cartesian = constants.major_angle_rad_cartesian_band6
        minor_angle_rad_cartesian = constants.minor_angle_rad_cartesian_band6
        line_length_arcsec = 1.4
    elif band == 4:
        centre_str = constants.centre_str_band4
        major_angle_rad_cartesian = constants.major_angle_rad_cartesian_band4
        minor_angle_rad_cartesian = constants.minor_angle_rad_cartesian_band4
        line_length_arcsec = 1.4  # Adjust this if needed for Band 4
    else:
        raise ValueError("Unsupported band. Only Band 4 and Band 6 are currently supported.")
        
    

    # Normalize angles to standard range
    major_angle_rad_cartesian = normalize_angle(major_angle_rad_cartesian)
    minor_angle_rad_cartesian = normalize_angle(minor_angle_rad_cartesian)


    # Convert centre position from string to pixel coordinates
    centre_pix = list(string_to_pixel(centre_str, StokesI_wcs))
    
#     centre_pix = [513.472, 510.353]


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

    # Extract image values along the defined axes
    major_data, _, major_offset_arcsec = extract_axis_data(major_x, major_y, 
                                                           data, 'major', 
                                                           centre_pix, gridsize, 
                                                           StokesI_header)
    
    minor_data, _, minor_offset_arcsec = extract_axis_data(minor_x, minor_y, 
                                                           data, 'minor', 
                                                           centre_pix, gridsize, 
                                                           StokesI_header)
    
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


    return major_data, major_offset_arcsec, minor_data, minor_offset_arcsec
# -----------------------------------------------------------------------------------------

