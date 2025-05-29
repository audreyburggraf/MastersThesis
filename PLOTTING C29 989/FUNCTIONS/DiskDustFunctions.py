import numpy as np 

def get_average_POLF_where_scattered(UniformRatios, POLF, tolerance=0.001):
    # Convert inputs to numpy arrays in case they aren't already
    UniformRatios = np.array(UniformRatios)
    POLF = np.array(POLF)
    
    # Create a mask where the ratio is within ±tolerance of 1
    mask = np.abs(UniformRatios - 1) <= tolerance
    
    # Use the mask to filter POLF_mJy values
    matching_POLF = POLF[mask]
    
    if matching_POLF.size == 0:
        return np.nan  # or raise an exception if preferred
    
    # Compute and return the average
    POLF_average = np.mean(matching_POLF)
    return POLF_average