import numpy as np

# ------------------------------------------------------
# This function is based on equation 7 from Sadavoy 2019 
# (https://iopscience.iop.org/article/10.3847/1538-4365/ab4257)
# ------------------------------------------------------
def galaxy_probability(S, band):
    """
    Differential number counts per deg^2
    S in mJy
    """
    if band == 'Band 6':
        S0 = 1.7 # mJy
    else:
        return print('Currently function only supports Band 6')
    
    dN_dS = 1800 * (S/S0)**(-2.08) * np.exp(-S/S0)
    
    return dN_dS
# ---------------------------------------------------