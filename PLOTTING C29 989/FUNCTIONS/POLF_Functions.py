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
def find_POLF_avg(band, POLF, StokesI_mJy, sigma_cutoff = 5):

    band_paths = {
        "Band 4": (constants.band4_data_folder_path, "constants_BAND4.csv"),
        "Band 4 nterms2": (constants.band4_nterms2_data_folder_path, "constants_BAND4_nterms2.csv"),
        "Band 4 nterms2 robust -1": (constants.band4_nterms2_robust_minus1_data_folder_path, "constants_BAND4_nterms2_robust_minus1.csv"),
        
        "Band 4 nterms2 smooth": (constants.band4_nterms2_smooth_data_folder_path, "constants_BAND4_nterms2_smooth.csv"),
        "Band 4 nterms2 smooth B6": (constants.band4_nterms2_smooth_B6_data_folder_path, "constants_BAND4_nterms2_smooth_B6.csv"),
        "Band 4 nterms2 smooth B6 B7": (constants.band4_nterms2_smooth_B6_B7_data_folder_path, "constants_BAND4_nterms2_smooth_B6_B7.csv"),

        "Band 5": (constants.band5_data_folder_path, "constants_BAND5.csv"),
        "Band 5 v0": (constants.band5_v0_data_folder_path, "constants_BAND5_v0.csv"),
        "Band 5 robust -2": (constants.band5_robust_minus2_data_folder_path, "constants_BAND5_robust_minus2.csv"),
        "Band 5 robust -1": (constants.band5_robust_minus1_data_folder_path, "constants_BAND5_robust_minus1.csv"),
        "Band 5 nterms2": (constants.band5_nterms2_data_folder_path, "constants_BAND5_nterms2.csv"),

        "Band 6": (constants.band6_data_folder_path, "constants_BAND6.csv"),
        "Band 6 smooth": (constants.band6_smooth_data_folder_path, "constants_BAND6_smooth.csv"),
        "Band 6 smooth B7": (constants.band6_smooth_B7_data_folder_path, "constants_BAND6_smooth_B7.csv"),

        "Band 7 nterms2": (constants.band7_nterms2_data_folder_path, "constants_BAND7_nterms2.csv"),
        "Band 7 nterms2 smooth": (constants.band7_nterms2_smooth_data_folder_path, "constants_BAND7_nterms2_smooth.csv"),
        "Band 7 nterms2 smooth B6": (constants.band7_nterms2_smooth_B6_data_folder_path, "constants_BAND7_nterms2_smooth_B6.csv"),
    }

    if band not in band_paths:
        raise ValueError(f"Unsupported band: {band}")

    folder, filename = band_paths[band]
    out_path = Path(folder) / filename


    
    threshold = sigma_cutoff * np.nanstd(StokesI_mJy)   # e.g. 5-sigma cut
    mask = StokesI_mJy > threshold

    POLF_region = POLF[mask]
    POLF_avg = np.nanmean(POLF_region)


    # Your computed value
    POLF_avg = float(POLF_avg)

    # Load existing file
    df = pd.read_csv(out_path)

    # Add / overwrite POLF_avg column
    df["POLF_avg_sigma"] = [POLF_avg]   # assumes 1-row CSV, like your others

    # Save back to disk
    df.to_csv(out_path, index=False)
# --------------------------------------------------------------------------------------



# Function to find 
# --------------------------------------------------------------------------------------
def find_gaussian_POLF(band, UniformRatios, POLF, POLF_err, tolerance=0.001):
    # Convert inputs to numpy arrays in case they aren't already
    UniformRatios = np.array(UniformRatios)
    POLF = np.array(POLF)
    
    # Create a mask where the ratio is within ±tolerance of 1
    mask = np.abs(UniformRatios - 1) <= tolerance
    
    # Use the mask to filter POLF_mJy values
    matching_POLF = POLF[mask]
    matching_POLF_err = POLF_err[mask]
    
    if matching_POLF.size == 0:
        return np.nan, mask  # or raise an exception if preferred
    
    # Compute and return the average
    POLF_average = np.mean(matching_POLF)
    POLF_err_average = np.mean(matching_POLF_err)
    
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
    

    
    return POLF_average, POLF_err_average
# --------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------
def find_POLF_at_max_POLI(StokesI_mJy, POLI_mJy, POLF, POLF_err, print_statements = False):
    # Find index of max POLI, ignoring NaNs
    y_max, x_max = np.unravel_index(np.nanargmax(POLI_mJy), POLI_mJy.shape)

    # Values at that location
    max_POLI = POLI_mJy[y_max, x_max]
    corresponding_I = StokesI_mJy[y_max, x_max]
    calc_avg_POLF = max_POLI / corresponding_I
    POLF_maxPOLI = POLF[y_max, x_max]
    POLF_err_maxPOLI = POLF_err[y_max, x_max]

    if print_statements:
        print(rf'The maximum POLI value is {max_POLI:.3f}, and it happens at the index ({y_max}, {x_max})')
        print(rf'At that index, the Stokes I value is {corresponding_I:.3f}')
        print(rf'At that index, the calculated POLF value is {calc_avg_POLF:.4f} or {calc_avg_POLF * 100:.2f}%')
        print(rf'At that index, the map POLF value is {POLF_maxPOLI:.4f}')
    
    return POLF_maxPOLI, POLF_err_maxPOLI
# --------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------
def find_POLF_at_max_StokesI(StokesI_mJy, POLI_mJy, POLF, POLF_err, print_statements = False):

    
# Find index of max Stokes I, ignoring NaNs
    y_max, x_max = np.unravel_index(np.nanargmax(StokesI_mJy), StokesI_mJy.shape)

    # Extract values at that index
    max_StokesI = StokesI_mJy[y_max, x_max]
    corresponding_POLI = POLI_mJy[y_max, x_max]
    calc_avg_POLF = corresponding_POLI / max_StokesI
    POLF_maxStokesI = POLF[y_max, x_max]
    POLF_err_maxStokesI = POLF_err[y_max, x_max]

    if print_statements:
        print(rf'The maximum Stokes I value is {max_StokesI:.3f}, at index ({y_max}, {x_max})')
        print(rf'At that index, the POLI value is {corresponding_POLI:.3f}')
        print(rf'Calculated POLF = {calc_avg_POLF:.4f} or {calc_avg_POLF * 100:.2f}%')
        print(rf'Map POLF value = {POLF_maxStokesI:.4f}')

    return POLF_maxStokesI, POLF_err_maxStokesI
# --------------------------------------------------------------------------------------





# --------------------------------------------------------------------------------------
def get_all_POLF_and_save(StokesI_mJy, POLI_mJy, POLF, POLF_err, UniformRatios, band, gaussian_tolerance = 0.001, sigma_cutoff = 5, print_statements = False, bias = 'debiased'):

    band_paths = {
    'Band 4': (constants.band4_data_folder_path, "constants_BAND4.csv"),
    'Band 4 nterms2': (constants.band4_nterms2_data_folder_path, "constants_BAND4_nterms2.csv"),
    'Band 4 nterms2 robust -1': (constants.band4_nterms2_robust_minus1_data_folder_path, 
                                 "constants_BAND4_nterms2_robust_minus1.csv"),
    'Band 4 nterms2 smooth': (constants.band4_nterms2_smooth_data_folder_path, "constants_BAND4_nterms2_smooth.csv"),
    'Band 4 nterms2 smooth B6': (constants.band4_nterms2_smooth_B6_data_folder_path, "constants_BAND4_nterms2_smooth_B6.csv"),
    'Band 4 nterms2 smooth B6 B7': (constants.band4_nterms2_smooth_B6_B7_data_folder_path,
                                    "constants_BAND4_nterms2_smooth_B6_B7.csv"),
        
    'Band 5': (constants.band5_data_folder_path, "constants_BAND5.csv"),
    'Band 5 v0': (constants.band5_v0_data_folder_path, "constants_BAND5_v0.csv"),
    'Band 5 robust -2': (constants.band5_robust_minus2_data_folder_path, "constants_BAND5_robust_minus2.csv"),
    'Band 5 robust -1': (constants.band5_robust_minus1_data_folder_path, "constants_BAND5_robust_minus1.csv"),
    'Band 5 nterms2': (constants.band5_nterms2_data_folder_path, "constants_BAND5_nterms2.csv"),

    'Band 6': (constants.band6_data_folder_path, "constants_BAND6.csv"),
    'Band 6 smooth': (constants.band6_smooth_data_folder_path, "constants_BAND6_smooth.csv"),
    'Band 6 smooth B7': (constants.band6_smooth_B7_data_folder_path, "constants_BAND6_smooth_B7.csv"),

    'Band 7 nterms2': (constants.band7_nterms2_data_folder_path, "constants_BAND7_nterms2.csv"),
    'Band 7 nterms2 smooth': (constants.band7_nterms2_smooth_data_folder_path, "constants_BAND7_nterms2_smooth.csv"),
    'Band 7 nterms2 smooth B6': (constants.band7_nterms2_smooth_B6_data_folder_path, "constants_BAND7_nterms2_smooth_B6.csv"),
    }

    if band not in band_paths:
        raise ValueError(
            "Band not supported"
        )

    folder, filename = band_paths[band]
    
    name, ext = filename.rsplit('.', 1)
    filename = f"{name}_{bias}.{ext}"


    path = folder + filename
    
    POLF_Gaussian, POLF_err_Gaussian = find_gaussian_POLF(band, UniformRatios, POLF, POLF_err, gaussian_tolerance)
    
    
    POLF_maxPOLI, POLF_err_maxPOLI = find_POLF_at_max_POLI(StokesI_mJy, POLI_mJy, POLF, POLF_err, print_statements)
    
    POLF_maxStokesI, POLF_err_maxStokesI = find_POLF_at_max_StokesI(StokesI_mJy, POLI_mJy, POLF, POLF_err, print_statements)
    
#     POLF_avg_sigma = find_POLF_avg(band, StokesI_mJy, POLF, sigma_cutoff)
    
    print(rf'POLF_Gaussian = {POLF_Gaussian * 100:.3f}')
    print(rf'POLF_maxPOLI = {POLF_maxPOLI * 100:.3f}')
    print(rf'POLF_maxStokesI = {POLF_maxStokesI * 100:.3f}')
#     print(rf'POLF_avg_sigma = {POLF_avg_sigma}')
    
    data_dict = {
        "POLF_Gaussian": [POLF_Gaussian],
        "POLF_err_Gaussian": [POLF_err_Gaussian],
        "POLF_maxPOLI": [POLF_maxPOLI],
        "POLF_err_maxPOLI": [POLF_err_maxPOLI],
        "POLF_maxStokesI": [POLF_maxStokesI],
        "POLF_err_maxStokesI": [POLF_err_maxStokesI],
#         "POLF_avg_sigma": [POLF_avg_sigma],
    }

    # Convert to DataFrame
    df = pd.DataFrame(data_dict)

    # Save to CSV
    df.to_csv(path, index=False)


# --------------------------------------------------------------------------------------





# Load POLF
# --------------------------------------------------------------------------------------
def load_POLF(bands, print_things = False, bias = 'old'):
    # Initialize storage
    results = {
        'Band': [],
        'POLF_Gaussian': [],
        'POLF_err_Gaussian': [],
        'POLF_maxPOLI': [],
        'POLF_err_maxPOLI': [],
        'POLF_maxStokesI': [],
        'POLF_err_maxStokesI': [],
#         'POLF_avg_sigma': [],
    }

    band_config = {
        "Band 4": ("band4_data_folder_path", "constants_BAND4.csv"),
        "Band 4 nterms2": ("band4_nterms2_data_folder_path", "constants_BAND4_nterms2.csv"),
        "Band 4 nterms2 robust -1": ("band4_nterms2_robust_minus1_data_folder_path", "constants_BAND4_nterms2_robust_minus1.csv"),
        
        "Band 4 nterms2 smooth": ("band4_nterms2_smooth_data_folder_path", "constants_BAND4_nterms2_smooth.csv"),
        "Band 4 nterms2 smooth B6": ("band4_nterms2_smooth_B6_data_folder_path", "constants_BAND4_nterms2_smooth_B6.csv"),
        "Band 4 nterms2 smooth B6 B7": ("band4_nterms2_smooth_B6_B7_data_folder_path",
                                        "constants_BAND4_nterms2_smooth_B6_B7.csv"),

        "Band 5": ("band5_data_folder_path", "constants_BAND5.csv"),
        "Band 5 v0": ("band5_v0_data_folder_path", "constants_BAND5_v0.csv"),
        "Band 5 robust -2": ("band5_robust_minus2_data_folder_path", "constants_BAND5_robust_minus2.csv"),
        "Band 5 robust -1": ("band5_robust_minus1_data_folder_path", "constants_BAND5_robust_minus1.csv"),
        "Band 5 nterms2": ("band5_nterms2_data_folder_path", "constants_BAND5_nterms2.csv"),

        "Band 6": ("band6_data_folder_path", "constants_BAND6.csv"),
        "Band 6 smooth": ("band6_smooth_data_folder_path", "constants_BAND6_smooth.csv"),
        "Band 6 smooth B7": ("band6_smooth_B7_data_folder_path", "constants_BAND6_smooth_B7.csv"),

        "Band 7": ("band7_data_folder_path", "constants_BAND7.csv"),
        "Band 7 nterms2": ("band7_nterms2_data_folder_path", "constants_BAND7_nterms2.csv"),
        "Band 7 nterms2 smooth": ("band7_nterms2_smooth_data_folder_path", "constants_BAND7_nterms2_smooth.csv"),
        "Band 7 nterms2 smooth B6": ("band7_nterms2_smooth_B6_data_folder_path", "constants_BAND7_nterms2_smooth_B6.csv"),
    }

    for band in bands:

        if print_things:
            print(f"Band = {band}")
            
         

        if band not in band_config:
            raise ValueError(f"Unsupported band: {band}")

        folder_var, filename = band_config[band]
        
#         if band != 'Band 6':
        if bias != 'old':
            name, ext = filename.rsplit('.', 1)
            filename = f"{name}_{bias}.{ext}"

        path = Path(getattr(constants, folder_var)) / filename

        df = pd.read_csv(path)
        
        print(rf'filename = {filename}')
        print(df.columns.tolist())

        results['Band'].append(band)
        results['POLF_Gaussian'].append(df.at[0, "POLF_Gaussian"] * 100)
        results['POLF_err_Gaussian'].append(df.at[0, "POLF_err_Gaussian"] * 100)
        results['POLF_maxPOLI'].append(df.at[0, "POLF_maxPOLI"] * 100)
        results['POLF_err_maxPOLI'].append(df.at[0, "POLF_err_maxPOLI"] * 100)
        results['POLF_maxStokesI'].append(df.at[0, "POLF_maxStokesI"] * 100)
        results['POLF_err_maxStokesI'].append(df.at[0, "POLF_err_maxStokesI"] * 100)
#         results['POLF_avg_simga'].append(df.at[0, "POLF_avg_simga"] * 100)
    

    print(rf' in load_POLF, results = {results}')
    # Convert to DataFrame
    df_POLF = pd.DataFrame(results)

    df_POLF["POLF_mean"] = df_POLF[["POLF_Gaussian", "POLF_err_Gaussian", "POLF_maxPOLI", "POLF_err_maxPOLI", "POLF_maxStokesI", "POLF_err_maxStokesI"]].mean(axis=1)

    
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
