import numpy as np

def calculate_polarized_intensity(StokesQ, StokesU):
    """
    Calculate the polarized intensity.

    Parameters:
        StokesQ (2D array): Stokes Q data in mJy or similar units.
        StokesU (2D array): Stokes U data in mJy or similar units.

    Returns:
        2D array: Polarized intensity in the same units as the inputs.
    """
    return np.sqrt(StokesQ**2 + StokesU**2)


def calculate_polarized_fraction(StokesQ, StokesU, StokesI):
    """
    Calculate the polarized fraction.

    Parameters:
        StokesQ (2D array): Stokes Q data in mJy or similar units.
        StokesU (2D array): Stokes U data in mJy or similar units.
        StokesI (2D array): Stokes I data in mJy or similar units.

    Returns:
        2D array: Polarized fraction (unitless).
    """
    polarized_intensity = calculate_polarized_intensity(StokesQ, StokesU)
    
    return polarized_intensity / StokesI


def calculate_polarization_angle(StokesQ, StokesU):
    """
    Calculate the polarization angle in radians.

    Parameters:
        StokesQ (2D array): Stokes Q data in Jy or similar units.
        StokesU (2D array): Stokes U data in Jy or similar units.

    Returns:
        2D array: Polarization angle in radians.
    """
    return 0.5 * np.arctan2(StokesU, StokesQ)
