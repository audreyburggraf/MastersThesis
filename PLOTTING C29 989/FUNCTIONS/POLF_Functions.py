import sys
import numpy as np
import pandas as pd

# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")  # Replace with the actual path

# Now you can import constants.py
import constants

# Function to find 
# --------------------------------------------------------------------------------------
def find_gaussian_POLF(UniformRatios, POLF, tolerance=0.001):
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
    
    count = 0

    for row in mask:
        for val in row:
            if val == True:
                count += 1

    
    return POLF_average
# --------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------
def find_POLF_at_max_POLI(StokesI_mJy, POLI_mJy, POLF, print_statements = False):
    # Find index of max POLI, ignoring NaNs
    y_max, x_max = np.unravel_index(np.nanargmax(POLI_mJy), POLI_mJy.shape)

    # Values at that location
    max_POLI = POLI_mJy[y_max, x_max]
    corresponding_I = StokesI_mJy[y_max, x_max]
    calc_avg_POLF = max_POLI / corresponding_I
    POLF_maxPOLI = POLF[y_max, x_max]

    if print_statements:
        print(rf'The maximum POLI value is {max_POLI:.3f}, and it happens at the index ({y_max}, {x_max})')
        print(rf'At that index, the Stokes I value is {corresponding_I:.3f}')
        print(rf'At that index, the calculated POLF value is {calc_avg_POLF:.4f} or {calc_avg_POLF * 100:.2f}%')
        print(rf'At that index, the map POLF value is {POLF_maxPOLI:.4f}')
    
    return POLF_maxPOLI
# --------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------
def find_POLF_at_max_StokesI(StokesI_mJy, POLI_mJy, POLF, print_statements = False):

    
# Find index of max Stokes I, ignoring NaNs
    y_max, x_max = np.unravel_index(np.nanargmax(StokesI_mJy), StokesI_mJy.shape)

    # Extract values at that index
    max_StokesI = StokesI_mJy[y_max, x_max]
    corresponding_POLI = POLI_mJy[y_max, x_max]
    calc_avg_POLF = corresponding_POLI / max_StokesI
    POLF_maxStokesI = POLF[y_max, x_max]

    if print_statements:
        print(rf'The maximum Stokes I value is {max_StokesI:.3f}, at index ({y_max}, {x_max})')
        print(rf'At that index, the POLI value is {corresponding_POLI:.3f}')
        print(rf'Calculated POLF = {calc_avg_POLF:.4f} or {calc_avg_POLF * 100:.2f}%')
        print(rf'Map POLF value = {POLF_maxStokesI:.4f}')

    return POLF_maxStokesI
# --------------------------------------------------------------------------------------





# --------------------------------------------------------------------------------------
def get_all_POLF_and_save(StokesI_mJy, POLI_mJy, POLF, UniformRatios, band, gaussian_tolerance = 0.001, print_statements = False):
    if band == 'Band 4':
        path = constants.band4_data_folder_path + "constants_BAND4.csv"
    elif band == 'Band 5':
        path = constants.band6_data_folder_path + "constants_BAND6.csv"
    elif band == 'Band 6':
        path = constants.band6_data_folder_path + "constants_BAND6.csv"
    elif band == 'Band 7':
        path = constants.band7_data_folder_path + "constants_BAND7.csv"
    else:
        raise ValueError("Only 'Band 4', 'Band 5', 'Band 6', and 'Band 7' are supported.")
        
    
    POLF_Gaussian = find_gaussian_POLF(UniformRatios, POLF, gaussian_tolerance)
    
    POLF_maxPOLI = find_POLF_at_max_POLI(StokesI_mJy, POLI_mJy, POLF, print_statements)
    
    POLF_maxStokesI = find_POLF_at_max_StokesI(StokesI_mJy, POLI_mJy, POLF, print_statements)
    
    
    data_dict = {
        "POLF_Gaussian": [POLF_Gaussian],
        "POLF_maxPOLI": [POLF_maxPOLI],
        "POLF_maxStokesI": [POLF_maxStokesI],
    }

    # Convert to DataFrame
    df = pd.DataFrame(data_dict)

    # Save to CSV
    df.to_csv(path, index=False)


# --------------------------------------------------------------------------------------
