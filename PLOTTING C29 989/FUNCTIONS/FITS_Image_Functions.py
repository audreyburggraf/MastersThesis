import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

# ---------------------------------------------------------------------------------------------------------
from astropy.io import fits
from astropy.wcs import WCS

from astropy.io import fits
from astropy.wcs import WCS

def read_in_file(file_name, dimensions = 4):
    """
    Reads a FITS file and extracts the header, data, and WCS (World Coordinate System) information.

    Parameters:
    ----------
    file_name : str
        Name of the FITS file (e.g., "c2d_989_StokesI_233GHz.fits").
    dimensions : int
        The expected dimensionality of the data (e.g., 4 for 4D data or 2 for 2D data).

    Returns:
    -------
    tuple
        A tuple containing:
        - header (astropy.io.fits.Header): The FITS file header.
        - data_4d (numpy.ndarray or None): The 4D data array from the FITS file (if applicable).
        - data_2d (numpy.ndarray): The 2D data slice from the FITS file.
        - wcs (astropy.wcs.WCS): The WCS object containing coordinate transformation information.

    Notes:
    -----
    - For 4D data, the 2D slice corresponds to the first Stokes parameter and first frequency slice.
    - The WCS object is initialized for 2D data (RA and Dec).
    - If `dimensions` is 2, only the 2D data will be extracted, and `data_4d` will be `None`.
    """
    # Open the FITS file with memory mapping enabled for large files
    hdu_list = fits.open(file_name, memmap=True)
    
    # Extract the header from the primary HDU (Header Data Unit)
    header = hdu_list[0].header
    
    # Initialize variables
    data_4d = None
    data_2d = None

    # Extract data based on the specified dimensions
    if dimensions == 4:
        # Extract the full 4D data array
        data_4d = hdu_list[0].data
        
        # Extract the 2D data slice (first Stokes parameter, first frequency slice)
        data_2d = data_4d[0, 0, :, :]
    elif dimensions == 2:
        # Directly extract the 2D data array
        data_2d = hdu_list[0].data
    else:
        raise ValueError("Invalid dimensions. Expected 2 or 4.")
    
    # Initialize a WCS object for 2D data (RA-Dec plane)
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = [header['CTYPE1'], header['CTYPE2']]  # Coordinate types (e.g., RA---TAN, DEC--TAN)
    wcs.wcs.crpix = [header['CRPIX1'], header['CRPIX2']]  # Reference pixel positions
    wcs.wcs.crval = [header['CRVAL1'], header['CRVAL2']]  # Reference world coordinates
    wcs.wcs.cdelt = [header['CDELT1'], header['CDELT2']]  # Coordinate increments
    
    # Return the extracted components
    return header, data_4d, data_2d, wcs

# ---------------------------------------------------------------------------------------------------------







# ---------------------------------------------------------------------------------------------------------
def get_beam_info(header):
        
    # Extract pixel sizes (in degrees per pixel)
    pixel_size_x_deg = abs(header.get('CDELT2'))  # degrees per pixel (X-axis) ** CDELT1 NEAGTIVE BC RA
    pixel_size_y_deg = abs(header.get('CDELT2'))  # degrees per pixel (Y-axis)
    

    # Extract beam parameters (in degrees)
    BMAJ_deg = header.get('BMAJ')  # beam major axis in degrees
    BMIN_deg = header.get('BMIN')  # beam minor axis in degrees
    BPA_deg_astronomy  = header.get('BPA')   # beam position angle in degrees
    
    BPA_deg_cartesian = astronomy_to_cartesian(BPA_deg_astronomy)
    
    
    # Convert beam size from degrees to pixels
    BMAJ_pix = BMAJ_deg / pixel_size_x_deg # beam major axis in pixels
    BMIN_pix = BMIN_deg / pixel_size_y_deg  # beam minor axis in pixels
    
    
    # Return a dictionary with the values
    return BMAJ_deg, BMIN_deg, BMAJ_pix, BMIN_pix, BPA_deg_astronomy, BPA_deg_cartesian

# ---------------------------------------------------------------------------------------------------------





# ---------------------------------------------------------------------------------------------------------
def length_in_pixels(length_AU, distance_pc, header):
    """
    Converts a reference length in AU to pixels along the X-axis using the FITS header's pixel scale.

    Parameters:
    - length_AU (float): The reference length in AU.
    - distance_pc (float): The distance to the object in parsecs.
    - header (dict): The FITS header that contains the CDELT1 value for pixel scale.

    Returns:
    - length_pixels_x (float): The reference length in pixels along the X-axis.
    """
    # Extract the pixel scale from the header (in degrees per pixel for X-axis)
    pixel_size_x_deg = abs(header.get('CDELT1', 0))  # degrees per pixel (X-axis)

    # Convert reference length (AU) to parsecs (using the given distance)
    # length_AU is the physical length, we calculate its angular size at the given distance (distance_pc)
    length_pc = length_AU / 206265.0  # Convert AU to parsecs, where 1 AU = 1/206265 parsecs
    
    # Convert parsecs to angular size in degrees
    length_deg = length_pc / distance_pc * (180 / np.pi)  # in degrees
    
    # Convert angular size to pixels
    length_pixels_x = length_deg / pixel_size_x_deg  # in pixels

    
    return length_pixels_x
# ---------------------------------------------------------------------------------------------------------







# ---------------------------------------------------------------------------------------------------------
def astronomy_to_cartesian(astronomy_angle):

        
    return astronomy_angle + 90
# ---------------------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------------------
def castesian_to_astronomy_rad(cartesian_angle_rad):

        
    return np.pi/2 - cartesian_angle_rad
# ---------------------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------------------
def astronomy_to_cartesian_rad(astronomy_angle_rad):

        
    return astronomy_angle_rad + np.pi/2
# ---------------------------------------------------------------------------------------------------------








# -----------------------------------------------------------------------------------------------------------------------------
def string_to_int_symbol(string_symbol):
    """
    Converts a string representing a symbol ('+' or '-') to its corresponding integer value.

    Args:
        string_symbol (str): A string symbol, either '+' or '-'.

    Returns:
        int: -1 for '-' and 1 for '+'.

    Raises:
        ValueError: If the input is neither '+' nor '-'.
    """
    
    if string_symbol == '-':
        return int(-1)  # Return -1 for minus symbol
    elif string_symbol == '+':
        return int(1)  # Return 1 for plus symbol
    else:
        raise ValueError("Input must be either '+' or '-'")  # Raise error if invalid symbol
# -----------------------------------------------------------------------------------------------------------------------------
def convert_ra_parts_into_degrees(ra_hours, ra_min, ra_sec):
    """
    Converts Right Ascension (RA) from hours, minutes, and seconds format into degrees.

    Args:
        ra_hours (int): The number of hours in the RA (range 0-24).
        ra_min (int): The number of minutes in the RA (range 0-59).
        ra_sec (float): The number of seconds in the RA (range 0-60).

    Returns:
        float: The RA value in degrees.

    Notes:
        - 1 hour = 15 degrees.
        - 1 minute = 1/60 of an hour = 15/60 degrees = 0.25 degrees.
        - 1 second = 1/3600 of an hour = 15/3600 degrees ≈ 0.0041667 degrees.
    """
    
    # Convert RA hours to degrees
    ra_hours_deg = ra_hours * 15

    # Convert RA minutes to degrees
    ra_min_deg = ra_min * (15 / 60)
    
    # Convert RA seconds to degrees
    ra_sec_deg = ra_sec * (15 / 3600)

    # Total RA in degrees
    ra_total_degrees = ra_hours_deg + ra_min_deg + ra_sec_deg
    
    return ra_total_degrees

# -----------------------------------------------------------------------------------------------------------------------------
def convert_dec_parts_into_degrees(dec_sign, dec_degrees, dec_arcmin, dec_arcsec):
    """
    Converts Declination (DEC) from degrees, arcminutes, and arcseconds format into degrees.

    Args:
        dec_sign (int): The sign of the declination, either +1 or -1.
        dec_degrees (int): The integer degrees of declination.
        dec_arcmin (int): The arcminutes of declination.
        dec_arcsec (int): The arcseconds of declination.

    Returns:
        float: The declination value in degrees, including the sign.

    Notes:
        - 1 arcminute = 1/60 degrees.
        - 1 arcsecond = 1/3600 degrees.
    """
    
    # Conversion factor: 1 arcminute = 1/60 degree, 1 arcsecond = 1/3600 degree
    arcminute_to_deg = 1 / 60.0
    arcsecond_to_deg = 1 / 3600.0
    
    # Convert arcminutes and arcseconds to degrees
    dec_arcmin_deg = dec_arcmin * arcminute_to_deg
    dec_arcsec_deg = dec_arcsec * arcsecond_to_deg
    
    # Calculate total declination in degrees
    dec_total_deg = dec_degrees + dec_arcmin_deg + dec_arcsec_deg
    
    # Apply the sign to the declination
    dec_deg_with_sign = dec_total_deg * dec_sign
    
    return dec_deg_with_sign
# -----------------------------------------------------------------------------------------------------------------------------
def find_ra_and_dec(ID):
    """
    Extracts Right Ascension (RA) and Declination (Dec) from a formatted source identifier 
    and converts them into degrees using pre-defined conversion functions.

    Args:
        ID (str): A source identifier in the format:
                  "Jhhmmss.sss+/-ddmmss.sss"
                  where:
                  - 'hhmmss.sss' represents RA (hours, minutes, seconds).
                  - '+/-ddmmss.sss' represents Dec (degrees, arcminutes, arcseconds).

    Returns:
        tuple: A tuple containing:
               - ra (float): Right Ascension in degrees.
               - dec (float): Declination in degrees.

    Example:
        For the source ID "J105339.700-771233.000":
            RA: 10h 53m 39.700s
            Dec: -77° 12' 33.000"
            The function returns:
            (163.41541666666666, -77.20916666666667)
    """

    # Parse RA components
    RA_hours = int(ID[1:3])  # Extract hours from position 1-3
    RA_min = int(ID[3:5])    # Extract minutes from position 3-5
    RA_seconds = float(ID[5:11])  # Extract seconds including an additional decimal place (position 5-11)

    # Parse Dec components
    Dec_sign = 1 if ID[11] == "+" else -1  # Determine sign (+1 or -1)
    Dec_deg = int(ID[12:14])               # Extract degrees from position 12-14
    Dec_arcmin = int(ID[14:16])            # Extract arcminutes from position 14-16
    Dec_arcsec = float(ID[16:21])          # Extract arcseconds including an additional decimal place (position 16-21)

    # Calculate RA in degrees
    ra = convert_ra_parts_into_degrees(RA_hours, RA_min, RA_seconds)

    # Calculate Dec in degrees
    dec = convert_dec_parts_into_degrees(Dec_sign, Dec_deg, Dec_arcmin, Dec_arcsec)

    return ra, dec


# -----------------------------------------------------------------------------------------------------------------------------
# def find_ra(ID):
#     """
#     Extracts Right Ascension (RA) and Declination (Dec) from a formatted source identifier 
#     and converts them into degrees using pre-defined conversion functions.

#     Args:
#         ID (str): A source identifier in the format:
#                   "Jhhmmss.ss+/-ddmmss.ss"
#                   where:
#                   - 'hhmmss.ss' represents RA (hours, minutes, seconds).
#                   - '+/-ddmmss.ss' represents Dec (degrees, arcminutes, arcseconds).

#     Returns:
#         tuple: A tuple containing:
#                - ra (float): Right Ascension in degrees.
#                - dec (float): Declination in degrees.

#     Notes:
#         - RA is converted using `convert_ra_parts_into_degrees`.
#         - Dec is converted using `convert_dec_parts_into_degrees`.

#     Example:
#         For the source ID "J105339.70-771233.00":
#             RA: 10h 53m 39.70s
#             Dec: -77° 12' 33.00"
#             The function returns:
#             (163.415, -77.209167)
#     """

#     # Parse RA components
#     RA_hours = int(ID[1:3])                       # Extract hours from position 1-3
#     RA_min = int(ID[3:5])                         # Extract minutes from position 3-5
#     RA_seconds = int(ID[5:7]) + 0.1 * int(ID[8])  # Extract seconds including tenths

#     # Calculate RA in degrees
#     ra = convert_ra_parts_into_degrees(RA_hours, RA_min, RA_seconds)

#     return ra
# -----------------------------------------------------------------------------------------------------------------------------
# def find_dec(ID):
#     """
#     Extracts Right Ascension (RA) and Declination (Dec) from a formatted source identifier 
#     and converts them into degrees using pre-defined conversion functions.

#     Args:
#         ID (str): A source identifier in the format:
#                   "Jhhmmss.ss+/-ddmmss.ss"
#                   where:
#                   - 'hhmmss.ss' represents RA (hours, minutes, seconds).
#                   - '+/-ddmmss.ss' represents Dec (degrees, arcminutes, arcseconds).

#     Returns:
#         tuple: A tuple containing:
#                - ra (float): Right Ascension in degrees.
#                - dec (float): Declination in degrees.

#     Notes:
#         - RA is converted using `convert_ra_parts_into_degrees`.
#         - Dec is converted using `convert_dec_parts_into_degrees`.

#     Example:
#         For the source ID "J105339.70-771233.00":
#             RA: 10h 53m 39.70s
#             Dec: -77° 12' 33.00"
#             The function returns:
#             (163.415, -77.209167)
#     """

#     # Parse Dec components
#     Dec_sign = string_to_int_symbol(ID[9])  # Get sign as an integer (+1 or -1)
#     Dec_deg = int(ID[10:12])                # Extract degrees from position 10-12
#     Dec_arcmin = int(ID[12:14])            # Extract arcminutes from position 12-14
#     Dec_arcsec = int(ID[14:17])            # Extract arcseconds from position 14-17

#     # Calculate Dec in degrees
#     dec = convert_dec_parts_into_degrees(Dec_sign, Dec_deg, Dec_arcmin, Dec_arcsec)

#     return dec  
# -----------------------------------------------------------------------------------------------------------------------------





# -----------------------------------------------------------------------------------------------------------------------------
def degrees_to_pixels(ra_deg, dec_deg, head):
    """
    Convert Right Ascension (RA) and Declination (Dec) in degrees to pixel coordinates 
    in a FITS image based on the header.

    Args:
        ra_deg (float): The RA position in degrees.
        dec_deg (float): The Dec position in degrees.
        head (dict): The FITS header containing WCS information, specifically:
                     - 'CRPIX1' (reference pixel for RA)
                     - 'CRPIX2' (reference pixel for Dec)
                     - 'CRVAL1' (reference RA in degrees)
                     - 'CRVAL2' (reference Dec in degrees)
                     - 'CDELT1' (RA increment per pixel, in degrees)
                     - 'CDELT2' (Dec increment per pixel, in degrees)

    Returns:
        tuple: A tuple containing:
               - x_pixel (float): X-coordinate in pixels.
               - y_pixel (float): Y-coordinate in pixels.
    """
    # Extract header values
    ra_ref = head['CRVAL1']  # Reference RA in degrees
    dec_ref = head['CRVAL2']  # Reference Dec in degrees
    x_ref = head['CRPIX1']  # Reference X pixel
    y_ref = head['CRPIX2']  # Reference Y pixel
    ra_inc = np.abs(head['CDELT1'])  # Increment per pixel in RA (degrees/pixel)
    dec_inc = head['CDELT2']  # Increment per pixel in Dec (degrees/pixel)

    # Calculate pixel coordinates
    x_pixel = (ra_deg - ra_ref) / ra_inc + x_ref
    y_pixel = (dec_deg - dec_ref) / dec_inc + y_ref

    return x_pixel, y_pixel
# -----------------------------------------------------------------------------------------------------------------------------
def string_to_pixel(string, wcs):
    """
    Converts a string containing RA and Dec coordinates into pixel coordinates.

    Parameters:
    ----------
    string : str
        A string containing the RA and Dec coordinates in a recognizable format.
    wcs : astropy.wcs.WCS
        A WCS object used to transform world coordinates (RA, Dec) into pixel coordinates.

    Returns:
    -------
    tuple
        A tuple containing the pixel coordinates (RA_pix, Dec_pix) as integers.
    
    Notes:
    -----
    - This function uses the `find_ra_and_dec` function to parse the RA and Dec from the input string.
    - The RA and Dec are assumed to be in degrees and are converted into a SkyCoord object.
    - The provided `wcs` parameter is used for the coordinate transformation.
    """
    # Parse the RA and Dec values from the input string
    RA_deg, Dec_deg = find_ra_and_dec(string)
    
    # Convert the RA and Dec into a SkyCoord object
    sky_coord = SkyCoord(ra=RA_deg * u.deg, dec=Dec_deg * u.deg, frame='icrs')

    # Use the WCS object to convert world coordinates (RA, Dec) to pixel coordinates
    pix = wcs.world_to_pixel(sky_coord)

    # Extract and convert pixel values to integers for use in image plots
    RA_pix = int(pix[0])
    Dec_pix = int(pix[1])
    
    return RA_pix, Dec_pix

# ---------------------------------------------------------------------------------------------------------
def arcsec_to_pixels(header, line_length_arcsec):
    """
    Converts a length from arcseconds to pixels using the FITS header plate scale.

    Parameters:
        header (astropy.io.fits.Header): FITS header object.
        line_length_arcsec (float): Length in arcseconds.

    Returns:
        float: Length in pixels.
    """
    # Extract plate scale from header (assuming CDELT1 is in degrees per pixel)
    if "CDELT1" in header:
        plate_scale = abs(header["CDELT1"]) * 3600  # Convert to arcsec/pixel
    else:
        raise ValueError("CDELT1 keyword not found in FITS header.")

    # Convert line length to pixels
    line_length_pix = line_length_arcsec / plate_scale
    return line_length_pix




# log stretching
# -----------------------------------------------------------------------------------------------------------------------------
def stretch(data, base=100, vmin=None, vmax=None):
    """
    Apply a logarithmic stretch to data with a specified base.

    Parameters:
    data (numpy.ndarray): The input data array.
    base (float): The base of the logarithm. Default is 100.
    vmin (float): Minimum value for stretching. If None, data minimum is used.
    vmax (float): Maximum value for stretching. If None, data maximum is used.

    Returns:
    numpy.ndarray: The transformed data with log stretch applied.
    """
    # Set default vmin and vmax if not specified
    if vmin is None:
        vmin = np.nanmin(data)
    if vmax is None:
        vmax = np.nanmax(data)
    
    # Normalize the data between 0 and 1
    normalized_data = (data - vmin) / (vmax - vmin)
    
    # Apply the log stretch using log base 100
    log_stretched_data = np.log(normalized_data * (base - 1) + 1) / np.log(base)
    
    return log_stretched_data
# -----------------------------------------------------------------------------------------------------------------------------
def unstretch(stretched_value, vmin, vmax, base=100):
    """
    Recover the original value from a log-stretched value.

    Parameters:
    stretched_value (float): The value from the log-stretched data (between 0 and 1).
    vmin (float): The minimum value of the original data.
    vmax (float): The maximum value of the original data.
    base (float): The base of the logarithm used for stretching. Default is 100.

    Returns:
    float: The corresponding original value before the stretch.
    """
    # Reverse the log stretch to get the normalized value
    normalized_data = (base**stretched_value - 1) / (base - 1)
    
    # Map the normalized value back to the original value range
    original_value = vmin + normalized_data * (vmax - vmin)
    
    return original_value
# -----------------------------------------------------------------------------------------------------------------------------






# def extract_axis_data(axis_x, axis_y, data_2d, centre_pix, gridsize):
#     """
#     Extracts data along a specified axis (major or minor) and returns the data along with valid positions.

#     Parameters:
#     - axis_x: Array of x positions along the axis (e.g., major_x or minor_x)
#     - axis_y: Array of y positions along the axis (e.g., major_y or minor_y)
#     - data_2d: 2D array representing the data (e.g., StokesI_data_2d_mJy)
#     - centre_pix: Tuple (RA_centre_pix, Dec_centre_pix) specifying the center pixel position
#     - gridsize: Tuple representing the dimensions of the grid (data_2d.shape)

#     Returns:
#     - axis_data: List of extracted data values
#     - offset_pixels: List of pixel offsets (positive and negative)
#     - offset_arcsec: List of arcsecond offsets (positive and negative)
#     """
#     axis_data = []
#     offset_pixels = []
#     offset_arcsec = []

#     RA_centre_pix, Dec_centre_pix = centre_pix  # Unpack center coordinates

#     for i in range(len(axis_x)):
#         # Round to nearest integer for pixel indexing
#         xi, yi = int(round(axis_x[i])), int(round(axis_y[i]))

#         # Ensure the indices are within the data grid bounds
#         if 0 <= xi < gridsize[1] and 0 <= yi < gridsize[0]:
#             # Append data value at (xi, yi)
#             axis_data.append(data_2d[yi, xi])

#             # Compute pixel offsets relative to the center
#             delta_x = xi - RA_centre_pix
#             delta_y = yi - Dec_centre_pix

#             # Compute radial offset in pixels
#             offset_pixel = np.sqrt(delta_x**2 + delta_y**2)
#             offset_pixels.append(offset_pixel)

#             # Convert pixel offset to arcseconds (assuming 1 pixel = 1 arcsec for now)
#             # You may need to multiply by plate scale if different
#             offset_arcsec.append(offset_pixel)  

#     return axis_data, offset_pixels, offset_arcsec



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




