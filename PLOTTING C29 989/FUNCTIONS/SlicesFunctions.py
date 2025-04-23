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




def run_slices(data_2d_mJy, StokesI_header, StokesI_wcs, length_of_carta, xmin, xmax, ymin, ymax, band):
    """
    Generate slice data along the major and minor axes for a given ALMA band.
    
    Parameters:
    - StokesI_header: FITS header of the Stokes I image
    - StokesI_wcs: WCS of the Stokes I image
    - length_of_carta: Number of points to sample along the axes (from CARTA)
    - xmin, xmax, ymin, ymax: Image bounds in pixels
    - band: ALMA band number (4 or 6)

    Returns:
    - major_data: Values along the major axis
    - major_offset_arcsec: Offsets along the major axis (arcsec)
    - minor_data: Values along the minor axis
    - minor_offset_arcsec: Offsets along the minor axis (arcsec)
    """

    # Get band-specific constants
    if band == 6:
        centre_str  = constants.centre_str_band6
        major_angle_rad_cartesian = constants.major_angle_rad_cartesian_band6
        minor_angle_rad_cartesian = constants.minor_angle_rad_cartesian_band6
        line_length_arcsec = 1.4

    elif band == 4:
        centre_str  = constants.centre_str_band4
        major_angle_rad_cartesian = constants.major_angle_rad_cartesian_band4
        minor_angle_rad_cartesian = constants.minor_angle_rad_cartesian_band4
        line_length_arcsec = 1.4  # <-- TODO: change this if needed

    else:
        raise ValueError("Unsupported band. Only Band 4 and Band 6 are currently supported.")
        
    # Normalize angles
    major_angle_rad_cartesian = normalize_angle(major_angle_rad_cartesian)
    minor_angle_rad_cartesian = normalize_angle(minor_angle_rad_cartesian)

    # Convert centre from string to pixel coordinates
    centre_pix = list(string_to_pixel(centre_str, StokesI_wcs))
#     print(f'My centre point: ({centre_pix[0]:.3f} px, {centre_pix[1]:.3f} px)')
#     print(f'Carta centre point: (858.123 px, 827.024 px)')

    # Convert line length to pixels
    line_length_pix = arcsec_to_pixels(StokesI_header, line_length_arcsec)
#     print(f"Line length in pixels: {line_length_pix:.2f}")

    # Define the number of sampling points (assumes CARTA has already been loaded)
    num_points = length_of_carta
    gridsize = data_2d_mJy.shape
#     print("the gridside is:", gridsize)
    
#     print(" ")
    # print([major_angle_rad_cartesian * (180/np.pi), minor_angle_rad_cartesian* (180/np.pi)])

    # Define x and y components along the major axis
    delta = np.linspace(-line_length_pix / 2, line_length_pix / 2, num_points)
    major_x = centre_pix[0] + delta * np.cos(major_angle_rad_cartesian)
    major_y = centre_pix[1] + delta * np.sin(major_angle_rad_cartesian)

    # Define x and y components along the minor axis
    minor_x = centre_pix[0] + delta * np.cos(minor_angle_rad_cartesian)
    minor_y = centre_pix[1] + delta * np.sin(minor_angle_rad_cartesian)
    
#     print("centre_pix[0] = ", centre_pix[0])
#     print("centre_pix[1] = ", centre_pix[1])
    
#     print(" ")
#     print(f'Major start: ({major_x[-1]:.3f}, {major_y[-1]:.3f})')
#     print(f'Major end:   ({major_x[0]:.3f}, {major_y[0]:.3f})')
#     print(f'Minor start: ({minor_x[-1]:.3f}, {minor_y[-1]:.3f})')
#     print(f'Minor end:   ({minor_x[0]:.3f}, {minor_y[0]:.3f})')

    # Extract data along the major axis
    major_data, _, major_offset_arcsec = extract_axis_data(
        major_x, major_y, 
        data_2d_mJy, 
        'major', 
        centre_pix, 
        gridsize, 
        StokesI_header
    )

    # Extract data along the minor axis
    minor_data, _, minor_offset_arcsec = extract_axis_data(
        minor_x, minor_y, 
        data_2d_mJy, 
        'minor', 
        centre_pix, 
        gridsize, 
        StokesI_header
    )

    return major_data, major_offset_arcsec, minor_data, minor_offset_arcsec
