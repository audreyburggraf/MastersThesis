

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
                  "Jhhmmss.s+/-ddmmss"
                  where:
                  - 'hhmmss.s' represents RA (hours, minutes, seconds).
                  - '+/-ddmmss' represents Dec (degrees, arcminutes, arcseconds).

    Returns:
        tuple: A tuple containing:
               - ra (float): Right Ascension in degrees.
               - dec (float): Declination in degrees.

    Notes:
        - RA is converted using `convert_ra_parts_into_degrees`.
        - Dec is converted using `convert_dec_parts_into_degrees`.

    Example:
        For the source ID "J105339.7-771233":
            RA: 10h 53m 39.7s
            Dec: -77° 12' 33.0"
            The function returns:
            (163.415, -77.209167)
    """

    # Parse RA components
    RA_hours = int(ID[1:3])         # Extract hours from position 1-3
    RA_min = int(ID[3:5])          # Extract minutes from position 3-5
    RA_seconds = int(ID[5:7]) + 0.1 * int(ID[8])  # Extract seconds including tenths

    # Parse Dec components
    Dec_sign = string_to_int_symbol(ID[9])  # Get sign as an integer (+1 or -1)
    Dec_deg = int(ID[10:12])                # Extract degrees from position 10-12
    Dec_arcmin = int(ID[12:14])            # Extract arcminutes from position 12-14
    Dec_arcsec = int(ID[14:17])            # Extract arcseconds from position 14-17

    # Calculate RA in degrees
    ra = convert_ra_parts_into_degrees(RA_hours, RA_min, RA_seconds)

    # Calculate Dec in degrees
    dec = convert_dec_parts_into_degrees(Dec_sign, Dec_deg, Dec_arcmin, Dec_arcsec)

    return ra, dec  
# -----------------------------------------------------------------------------------------------------------------------------
def find_ra(ID):
    """
    Extracts Right Ascension (RA) and Declination (Dec) from a formatted source identifier 
    and converts them into degrees using pre-defined conversion functions.

    Args:
        ID (str): A source identifier in the format:
                  "Jhhmmss.s+/-ddmmss"
                  where:
                  - 'hhmmss.s' represents RA (hours, minutes, seconds).
                  - '+/-ddmmss' represents Dec (degrees, arcminutes, arcseconds).

    Returns:
        tuple: A tuple containing:
               - ra (float): Right Ascension in degrees.
               - dec (float): Declination in degrees.

    Notes:
        - RA is converted using `convert_ra_parts_into_degrees`.
        - Dec is converted using `convert_dec_parts_into_degrees`.

    Example:
        For the source ID "J105339.7-771233":
            RA: 10h 53m 39.7s
            Dec: -77° 12' 33.0"
            The function returns:
            (163.415, -77.209167)
    """

    # Parse RA components
    RA_hours = int(ID[1:3])                       # Extract hours from position 1-3
    RA_min = int(ID[3:5])                         # Extract minutes from position 3-5
    RA_seconds = int(ID[5:7]) + 0.1 * int(ID[8])  # Extract seconds including tenths

    # Calculate RA in degrees
    ra = convert_ra_parts_into_degrees(RA_hours, RA_min, RA_seconds)

    return ra
# -----------------------------------------------------------------------------------------------------------------------------
def find_dec(ID):
    """
    Extracts Right Ascension (RA) and Declination (Dec) from a formatted source identifier 
    and converts them into degrees using pre-defined conversion functions.

    Args:
        ID (str): A source identifier in the format:
                  "Jhhmmss.s+/-ddmmss"
                  where:
                  - 'hhmmss.s' represents RA (hours, minutes, seconds).
                  - '+/-ddmmss' represents Dec (degrees, arcminutes, arcseconds).

    Returns:
        tuple: A tuple containing:
               - ra (float): Right Ascension in degrees.
               - dec (float): Declination in degrees.

    Notes:
        - RA is converted using `convert_ra_parts_into_degrees`.
        - Dec is converted using `convert_dec_parts_into_degrees`.

    Example:
        For the source ID "J105339.7-771233":
            RA: 10h 53m 39.7s
            Dec: -77° 12' 33.0"
            The function returns:
            (163.415, -77.209167)
    """

    # Parse Dec components
    Dec_sign = string_to_int_symbol(ID[9])  # Get sign as an integer (+1 or -1)
    Dec_deg = int(ID[10:12])                # Extract degrees from position 10-12
    Dec_arcmin = int(ID[12:14])            # Extract arcminutes from position 12-14
    Dec_arcsec = int(ID[14:17])            # Extract arcseconds from position 14-17

    # Calculate Dec in degrees
    dec = convert_dec_parts_into_degrees(Dec_sign, Dec_deg, Dec_arcmin, Dec_arcsec)

    return dec  
# -----------------------------------------------------------------------------------------------------------------------------

