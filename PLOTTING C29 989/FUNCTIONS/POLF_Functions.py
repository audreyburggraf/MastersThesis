import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import minimize

# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")  # Replace with the actual path

# Now you can import constants.py
import constants


from DustModelFunctions import * 

POLF_columns = ['POLF_Gaussian', 'POLF_maxPOLI', 'POLF_maxStokesI']



# Calculate POLF average
# --------------------------------------------------------------------------------------
def find_POLF_avg(band, POLF, StokesI_mJy, path, sigma_cutoff = 5):

    if band == "Band 4":
        out_path = Path(path) / "constants_BAND4.csv"

    elif band == "Band 4 nterms2":
        out_path = Path(path) / "constants_BAND4_nterms2.csv"
        
    elif band == "Band 4 nterms2 smooth":
        out_path = Path(path) / "constants_BAND4_nterms2_smooth.csv"

    elif band == "Band 5":
        out_path = Path(path) / "constants_BAND5.csv"
    
    elif band == "Band 5 v0":
        out_path = Path(path) / "constants_BAND5_v0.csv"
    
    elif band == "Band 5 robust -2":
        out_path = Path(path) / "constants_BAND5_robust_minus2.csv"
    
    elif band == "Band 5 robust -1":
        out_path = Path(path) / "constants_BAND5_robust_minus1.csv"
        
    elif band == "Band 5 nterms2":
        out_path = Path(path) / "constants_BAND5_nterms2.csv"

    elif band == "Band 6":
        out_path = Path(path) / "constants_BAND6.csv"
    
    elif band == "Band 6 smooth":
        out_path = Path(path) / "constants_BAND6_smooth.csv"

    elif band == "Band 7 nterms2":
        out_path = Path(path) / "constants_BAND7_nterms2.csv"
    
    elif band == "Band 7 nterms2 smooth":
        out_path = Path(path) / "constants_BAND7_nterms2_smooth.csv"

    else:
        raise ValueError(
            "Function currently only accepts: "
            "'Band 4', 'Band 4 nterms2', 'Band 4 nterms2 smooth', 'Band 5', 'Band 5 v0', 'Band 5 robust -2', 'Band 5 robust -1', 'Band 5 nterms2', 'Band 6', 'Band 6 smooth', 'Band 7 nterms2', 'Band 7 nterms2 smooth'")

    
    threshold = sigma_cutoff * np.nanstd(StokesI_mJy)   # e.g. 5-sigma cut
    mask = StokesI_mJy > threshold

    POLF_region = POLF[mask]
    POLF_avg = np.nanmean(POLF_region)


    # Your computed value
    POLF_avg = float(POLF_avg)

    # Load existing file
    df = pd.read_csv(out_path)

    # Add / overwrite POLF_avg column
    df["POLF_avg"] = [POLF_avg]   # assumes 1-row CSV, like your others

    # Save back to disk
    df.to_csv(out_path, index=False)
# --------------------------------------------------------------------------------------



# Function to find 
# --------------------------------------------------------------------------------------
def find_gaussian_POLF(band, UniformRatios, POLF, tolerance=0.001):
    # Convert inputs to numpy arrays in case they aren't already
    UniformRatios = np.array(UniformRatios)
    POLF = np.array(POLF)
    
    # Create a mask where the ratio is within ±tolerance of 1
    mask = np.abs(UniformRatios - 1) <= tolerance
    
    # Use the mask to filter POLF_mJy values
    matching_POLF = POLF[mask]
    
    if matching_POLF.size == 0:
        return np.nan, mask  # or raise an exception if preferred
    
    # Compute and return the average
    POLF_average = np.mean(matching_POLF)
    
    count = 0

    for row in mask:
        for val in row:
            if val == True:
                count += 1
    
#     if band == "Band 4":
#         out_path = Path(path) / "POLF_Gaussisn_mask_BAND4.csv"
#         np.save(out_path, mask)


#     elif band == "Band 4 nterms2":
#         out_path = Path(path) / "POLF_Gaussisn_mask_BAND4_nterms2.csv"
#         np.save(out_path, mask)


#     elif band == "Band 5":
#         out_path = Path(path) / "POLF_Gaussisn_mask_BAND5.csv"
#         np.save(out_path, mask)


#     elif band == "Band 6":
#         out_path = Path(path) / "POLF_Gaussisn_mask_BAND6.csv"
#         np.save(out_path, mask)


#     elif band == "Band 7 nterms2":
#         out_path = Path(path) / "POLF_Gaussisn_mask_BAND7_nterms2.csv"
#         np.save(out_path, mask)


#     else:
#         raise ValueError(
#             "Function currently only accepts: "
#             "'Band 4', 'Band 4 nterms2', 'Band 5', 'Band 6', 'Band 7 nterms2'")
    

    
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
        
    elif band == 'Band 4 nterms2':
        path = constants.band4_nterms2_data_folder_path + "constants_BAND4_nterms2.csv"
    
    elif band == 'Band 4 nterms2 smooth':
        path = constants.band4_nterms2_smooth_data_folder_path + "constants_BAND4_nterms2_smooth.csv"
        
    elif band == 'Band 5':
        path = constants.band5_data_folder_path + "constants_BAND5.csv"
    elif band == 'Band 5 v0':
        path = constants.band5_v0_data_folder_path + "constants_BAND5_v0.csv"
    
    elif band == 'Band 5 robust -2':
        path = constants.band5_robust_minus2_data_folder_path + "constants_BAND5_robust_minus2.csv"
        
    elif band == 'Band 5 robust -1':
        path = constants.band5_robust_minus1_data_folder_path + "constants_BAND5_robust_minus1.csv"
        
    elif band == 'Band 5 nterms2':
        path = constants.band5_nterms2_data_folder_path + "constants_BAND5_nterms2.csv"
    
    elif band == 'Band 6':
        path = constants.band6_data_folder_path + "constants_BAND6.csv"
        
        
    elif band == 'Band 6 smooth':
        path = constants.band6_smooth_data_folder_path + "constants_BAND6_smooth.csv"
    
    elif band == 'Band 7 nterms2':
        path = constants.band7_nterms2_data_folder_path + "constants_BAND7_nterms2.csv"
        
    elif band == 'Band 7 nterms2 smooth':
        path = constants.band7_nterms2_smooth_data_folder_path + "constants_BAND7_nterms2_smooth.csv"
        
        
    else:
        raise ValueError("Only 'Band 4', 'Band 4 nterms2', 'Band 4 nterms2 smooth', 'Band 5', 'Band 5 nterms2', 'Band 5 robust -2', 'Band 5 robust -1', 'Band 6', 'Band 6 smooth', 'Band 7 nterms2' and 'Band 7 nterms2 smooth' are supported.")
        
    
    POLF_Gaussian = find_gaussian_POLF(band, UniformRatios, POLF, gaussian_tolerance)
    
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




# Load POLF
# --------------------------------------------------------------------------------------
def load_POLF(bands, print_things = False):
    # Initialize storage
    results = {
        'Band': [],
        'POLF_Gaussian': [],
        'POLF_maxPOLI': [],
        'POLF_maxStokesI': []
    }

    for band in bands:
        if print_things:
            print(rf'Band = {band}')
        # Case 1: Band 4 nterms2
        if band == "Band 4 nterms2":
            var_name = "band4_nterms2_data_folder_path"
            path = Path(getattr(constants, var_name)) / "constants_BAND4_nterms2.csv"
            
 
        elif band == "Band 4 nterms2 smooth":
            var_name = "band4_nterms2_smooth_data_folder_path"
            path = Path(getattr(constants, var_name)) / "constants_BAND4_nterms2_smooth.csv"

        # Case 2: Band 7 nterms2
        elif band == "Band 7 nterms2":
            var_name = "band7_nterms2_data_folder_path"
            path = Path(getattr(constants, var_name)) / "constants_BAND7_nterms2.csv"
            
        elif band == "Band 7 nterms2 smooth":
            var_name = "band7_nterms2_smooth_data_folder_path"
            path = Path(getattr(constants, var_name)) / "constants_BAND7_nterms2_smooth.csv"
            
        # Case 3: Band 5 v0
        elif band == "Band 5 v0":
            var_name = "band5_v0_data_folder_path"
            path = Path(getattr(constants, var_name)) / "constants_BAND5_v0.csv"
            
        # Case 4: Band 5 robust -2
        elif band == "Band 5 robust -2":
            var_name = "band5_robust_minus2_data_folder_path"
            path = Path(getattr(constants, var_name)) / "constants_BAND5_robust_minus1.csv"
        
        # Case 5: Band 5 robust -1
        elif band == "Band 5 robust -1":
            var_name = "band5_robust_minus2_data_folder_path"
            path = Path(getattr(constants, var_name)) / "constants_BAND5_robust_minus1.csv"
            
            
        elif band == "Band 6 smooth":
            var_name = "band6_smooth_data_folder_path"
            path = Path(getattr(constants, var_name)) / "constants_BAND6_smooth.csv"

        # Case 6: Normal bands: "Band 4", "Band 5", "Band 6", "Band 7", ...
        else:
            # Extract the number from "Band X"
            band_num = int(band.split()[-1])   # "Band 4" -> 4

            var_name = f"band{band_num}_data_folder_path"
            path = Path(getattr(constants, var_name)) / f"constants_BAND{band_num}.csv"
   
        

        df = pd.read_csv(path)

        results['Band'].append(band)
        results['POLF_Gaussian'].append(df.at[0, "POLF_Gaussian"] * 100)
        results['POLF_maxPOLI'].append(df.at[0, "POLF_maxPOLI"] * 100)
        results['POLF_maxStokesI'].append(df.at[0, "POLF_maxStokesI"] * 100)
    

    # Convert to DataFrame
    df_POLF = pd.DataFrame(results)

    df_POLF["POLF_mean"] = df_POLF[["POLF_Gaussian", "POLF_maxPOLI", "POLF_maxStokesI"]].mean(axis=1)

    
    return df_POLF
# --------------------------------------------------------------------------------------



# Find best sf
# --------------------------------------------------------------------------------------
def find_best_sf(a_max_f_dist_micron, P_times_omega, df_POLF, print_results = False):
    
    # Start a df and test
    df = pd.DataFrame({'a_max micron': a_max_f_dist_micron})
    
    
    
    # Loop over all columns of POLF
    for j in range(len(POLF_columns)):
        # Extract the column name
        col_name = POLF_columns[j]
        if print_results:
            print(f'We are minimizing for {col_name}:')

        # Choose the test POLF
        POLF_test = df_POLF[col_name]

        scale_factors = []
        diff_array = []
        
        # Minimize the sum of squares
        for i in range(len(a_max_f_dist_micron)):
            result = minimize(
                dust_sum_of_sq,
                x0=0.01,
                args=(P_times_omega[i, :], POLF_test),
                bounds=[(0, None)]
            )

            sf_opt = result.x[0]
            err = result.fun
            
            # Save the data
            scale_factors.append(sf_opt)
            diff_array.append(err)

        # Add columns to main df
        df[f'sf_{col_name}'] = scale_factors
        df[f'error_{col_name}'] = diff_array
        
        
        best_index_array = []
        best_error_array = []
        best_a_max_f_array = []
        best_sf_array = []

    for col_name in POLF_columns:
        if print_results:
            print(f'We are now looking at the results for {col_name}:')

        best_index = np.argmin(df[f'error_{col_name}'])
        best_error = df[f'error_{col_name}'][best_index]
        best_a_max_f = a_max_f_dist_micron[best_index]
        best_sf = df[f'sf_{col_name}'][best_index]

        best_index_array.append(best_index)
        best_error_array.append(best_error)
        best_a_max_f_array.append(best_a_max_f)
        best_sf_array.append(best_sf)

        if print_results:
            print(f"   The index with the lowest error from minimize is    : {best_index}")
            print(f"   At this index, the corresponding error is           : {best_error:.3f}")
            print(f"   At this index, the corresponding a_max * f is       : {best_a_max_f:.1f} micron")
            print(f"   The scale factor at that a_max * f is               : {best_sf:.5f}")
            print(' ')

    if print_results:    
        print(' ')
        print(' ')

    best_POLF_index = np.argmin(best_error_array)
    best_POLF_column = POLF_columns[best_POLF_index]
    best_POLF = df_POLF[best_POLF_column]
    best_index = best_index_array[best_POLF_index]
    best_error = best_error_array[best_POLF_index]
    best_a_max_f = best_a_max_f_array[best_POLF_index]
    best_sf = best_sf_array[best_POLF_index]


    if print_results:
        print(f'Out of the {len(df_POLF)} POLF options:')
        print(f"   Lowest error occurs at POLF index of                 : {best_POLF_index}")
        print(f"   This POLF index corresponds to                       : {best_POLF_column}")   
        print(f"   The POLFs at this best POLF column are               : {best_POLF[0]:.3f},  {best_POLF[1]:.3f},  {best_POLF[2]:.3f}")
        print(f"   The index with the lowest error from minimize is     : {best_index}")
        print(f"   At this index, the corresponding error is            : {best_error:.3f}")
        print(f"   At this index, the corresponding a_max * f is        : {best_a_max_f:.1f} micron")
        print(f"   The scale factor at that a_max * f is                : {best_sf:.5f}")


    
    
    
    
    return best_a_max_f, best_POLF, best_sf
# --------------------------------------------------------------------------------------






# # --------------------------------------------------------------------------------------
# def find_sf_v2(a_max_f_dist_micron, P_times_omega, df_POLF, print_results = False):
    
#     # Start a df and test
#     df = pd.DataFrame({'a_max micron': a_max_f_dist_micron})
    
    
#     # We are only going to be using POLF_mean
#     POLF_values = 
   
    
#     # Loop over all columns of POLF
#     for j in range(len(POLF_columns)):
#         # Extract the column name
#         col_name = POLF_columns[j]
#         if print_results:
#             print(f'We are minimizing for {col_name}:')

#         # Choose the test POLF
#         POLF_test = df_POLF[col_name]

#         scale_factors = []
#         diff_array = []
        
#         # Minimize the sum of squares
#         for i in range(len(a_max_f_dist_micron)):
#             result = minimize(
#                 dust_sum_of_sq,
#                 x0=0.01,
#                 args=(P_times_omega[i, :], POLF_test),
#                 bounds=[(0, None)]
#             )

#             sf_opt = result.x[0]
#             err = result.fun
            
#             # Save the data
#             scale_factors.append(sf_opt)
#             diff_array.append(err)

#         # Add columns to main df
#         df[f'sf_{col_name}'] = scale_factors
#         df[f'error_{col_name}'] = diff_array
        
        
#         best_index_array = []
#         best_error_array = []
#         best_a_max_f_array = []
#         best_sf_array = []

#     for col_name in POLF_columns:
#         if print_results:
#             print(f'We are now looking at the results for {col_name}:')

#         best_index = np.argmin(df[f'error_{col_name}'])
#         best_error = df[f'error_{col_name}'][best_index]
#         best_a_max_f = a_max_f_dist_micron[best_index]
#         best_sf = df[f'sf_{col_name}'][best_index]

#         best_index_array.append(best_index)
#         best_error_array.append(best_error)
#         best_a_max_f_array.append(best_a_max_f)
#         best_sf_array.append(best_sf)

#         if print_results:
#             print(f"   The index with the lowest error from minimize is    : {best_index}")
#             print(f"   At this index, the corresponding error is           : {best_error:.3f}")
#             print(f"   At this index, the corresponding a_max * f is       : {best_a_max_f:.1f} micron")
#             print(f"   The scale factor at that a_max * f is               : {best_sf:.5f}")
#             print(' ')

#     if print_results:    
#         print(' ')
#         print(' ')

#     best_POLF_index = np.argmin(best_error_array)
#     best_POLF_column = POLF_columns[best_POLF_index]
#     best_POLF = df_POLF[best_POLF_column]
#     best_index = best_index_array[best_POLF_index]
#     best_error = best_error_array[best_POLF_index]
#     best_a_max_f = best_a_max_f_array[best_POLF_index]
#     best_sf = best_sf_array[best_POLF_index]


#     if print_results:
#         print(f'Out of the {len(df_POLF)} POLF options:')
#         print(f"   Lowest error occurs at POLF index of                 : {best_POLF_index}")
#         print(f"   This POLF index corresponds to                       : {best_POLF_column}")   
#         print(f"   The POLFs at this best POLF column are               : {best_POLF[0]:.3f},  {best_POLF[1]:.3f},  {best_POLF[2]:.3f}")
#         print(f"   The index with the lowest error from minimize is     : {best_index}")
#         print(f"   At this index, the corresponding error is            : {best_error:.3f}")
#         print(f"   At this index, the corresponding a_max * f is        : {best_a_max_f:.1f} micron")
#         print(f"   The scale factor at that a_max * f is                : {best_sf:.5f}")


    
    
    
    
#     return best_a_max_f, best_POLF, best_sf
# # --------------------------------------------------------------------------------------
