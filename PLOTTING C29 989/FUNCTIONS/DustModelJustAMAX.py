import numpy as np 
from scipy.interpolate import interp1d
import pandas as pd


# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
def calculate_scale_factors(a_max_vals,
                            Pw,
                            POLF_obs_all,
                            fit_indices):
    """
    Calculate the median scale factor and its scatter for every a_max.

    Parameters
    ----------
    a_max_vals : ndarray
        Array of a_max values (only used for length).
    Pw : ndarray
        Shape (N_a_max, N_bands). Model P*omega values.
    POLF_obs_all : ndarray
        Observed polarization fractions for all bands.
    fit_indices : list
        Indices of the bands to include in the fit.

    Returns
    -------
    sf_medians : ndarray
        Median scale factor at each a_max.
    sf_stds : ndarray
        Standard deviation of the scale factor at each a_max.
    """

    sf_medians = []
    sf_stds = []
    
    print(f'fit_indicies: {fit_indices}')

    for i in range(len(a_max_vals)):

        # Model values for one a_max
        Pw_model = Pw[i, :]

        # Only include requested bands
        Pw_fit = Pw_model[fit_indices]
        POLF_fit = POLF_obs_all[fit_indices]

        # Avoid divide-by-zero
        valid = Pw_fit > 0

        sf_i = POLF_fit[valid] / Pw_fit[valid]

        sf_medians.append(np.median(sf_i))
        sf_stds.append(np.std(sf_i))

    return np.array(sf_medians), np.array(sf_stds)
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
def calculate_chi_squared_for_sf(Pw,
                          best_idx,
                          best_sf,
                          POLF_obs_all,
                          POLF_err_obs_all,
                          bands):
    """
    Calculate the chi-squared value for the best-fitting model.

    Parameters
    ----------
    Pw : ndarray
        Shape (N_a_max, N_bands). Model P*omega values.
    best_idx : int
        Index of the best a_max.
    best_sf : float
        Best scale factor.
    POLF_obs_all : ndarray
        Observed polarization fractions.
    bands : list
        List of band names (e.g. ["Band 4", "Band 5", ...]).

    Returns
    -------
    chi_sq : float
    """

    chi_sq = 0

    for band_idx, band in enumerate(bands):
#         print(band)

        if band == "Band 5":
            continue

        POLF_obs     = POLF_obs_all[band_idx]
        POLF_obs_err = POLF_err_obs_all[band_idx]
        POLF_model   = Pw[best_idx, band_idx] * best_sf

        chi_sq += ((POLF_obs - POLF_model) / POLF_obs_err) ** 2

    return chi_sq
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
def load_Pomega_from_csv(fname,
                         folder,
                         bands_naming):

    df = pd.read_csv(folder + fname)

    a_max_vals = df["a_max_micron"].values

    Pw = np.column_stack([
        df[f"P_times_omega_{b}"].values
        for b in bands_naming
    ])

    return a_max_vals, Pw
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# I am updating this so we can say what POLF we want to use 
def find_sf_v4(bands,
               bands_naming,
               bands_included_in_fit,
               csv_files,
               P_omega_data_folder_path,
               df_POLF,
               POLF_index,
               print_results=False,
               trouble = False,
               trouble_f = 0.75,
               choose_trouble_index = False,
               trouble_index = 2):
    
    polf_columns = {
        'gaussian': 'POLF_Gaussian',
        'max Stokes I': 'POLF_maxStokesI',
#         'max Stokes I err': 'POLF_err_maxStokesI',
        'POLI': 'POLF_maxPOLI',
#         'POLI_err': 'POLF_err_maxPOLI',
        'mean': 'POLF_mean'
    }
    
    polf_err_columns = {
    #'gaussian': 'POLF_err_Gaussian',
    'max Stokes I': 'POLF_err_maxStokesI',
    'POLI': 'POLF_err_maxPOLI',
    #'mean': 'POLF_err_mean'   # if you have this column
    }
    
    
    if POLF_index not in polf_columns:
        raise ValueError(
            "POLF_index must be 'gaussian', 'poli', 'avg', 'mean', or 'max Stokes I'"
        )
        
        
    
    print(rf'The Bands we are looking at are: {bands}')
    print(rf"The Bands included in the fit are: {bands_included_in_fit}")
 
    # Band names → indices (0, 1, 2, etc)
    band_to_index = {band: i for i, band in enumerate(bands)}
    fit_indices = [band_to_index[b] for b in bands_included_in_fit]
    
    
    # This stores the data for each f value
    a_max_dist_micron_by_f = {}
    P_times_omega_by_f = {}
    a_max_best_by_f = {}
    POLF_obs_by_f = {}
    best_sf_by_f = {}
    best_idx_by_f = {}
    chi_sq_by_f = {}



    POLF_obs_all = df_POLF[polf_columns[POLF_index]].values
    POLF_err_obs_all = df_POLF[polf_err_columns[POLF_index]].values
   

    

    # LOOP OVER FILES (f values)
    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    for f, fname in csv_files.items():

        a_max_vals, Pw = load_Pomega_from_csv(
            fname,
            P_omega_data_folder_path,
            bands_naming,
        )
        

        sf_medians, sf_stds = calculate_scale_factors(
                                a_max_vals,
                                Pw,
                                POLF_obs_all,
                                fit_indices,
                            )
            

        
        
        # Choose best a_max where scatter is minimized
        best_idx = np.argmin(sf_stds)
        


        # Fitting Band 4 and 7 only at f = 0.75 gave me a weird result
        # Want to look at it closer
        # ----------------------------------------
        if trouble:
            if f == trouble_f:
                sorted_indices = np.argsort(sf_stds)

                
                # Print results 
                print(rf'The best a_max is: {a_max_vals[best_idx]}')
                print(rf'At this value, the best sf is : {sf_medians[best_idx]}')
                print(' ')
                
                # Get the next two lwest
                second_best_idx = sorted_indices[1]
                third_best_idx = sorted_indices[2]

                
                # Print results 
                print(rf'The second best a_max is: {a_max_vals[second_best_idx]}')
                print(rf'At this value, the best sf is : {sf_medians[second_best_idx]}')
                print(' ')
                
                # Print results 
                print(rf'The third best a_max is: {a_max_vals[third_best_idx]}')
                print(rf'At this value, the best sf is : {sf_medians[third_best_idx]}')
                print(' ')
                
                
                if choose_trouble_index:
                    best_idx = sorted_indices[trouble_index]
                    
        best_sf = sf_medians[best_idx]
                    
        # ----------------------------------------
        # Calculate chi^2 for each f_value
        # ----------------------------------------
        chi_sq = calculate_chi_squared_for_sf(
            Pw,
            best_idx,
            best_sf,
            POLF_obs_all,
            POLF_err_obs_all,
            bands,
        )
        # ---------------------------------------------
        # ---------------------------------------------
        # ----------------------------------------
                    
        # ----------------------------------------
        
        # Save the information to each by f array
        a_max_dist_micron_by_f[f] = a_max_vals
        P_times_omega_by_f[f] = Pw
        a_max_best_by_f[f] = a_max_vals[best_idx]
        best_sf_by_f[f] = best_sf
        best_idx_by_f[f] = best_idx
        POLF_obs_by_f[f] = np.array([POLF_obs_all[band_to_index[b]] for b in bands])
        chi_sq_by_f[f] = chi_sq
    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    # End of looping over all f values
        

    if print_results:
        print("Best idx:", best_idx_by_f)
        print("Best a_max:", a_max_best_by_f)
        print("Best SF:", best_sf_by_f)
        print("Best chi^2:", chi_sq_by_f)

    # return values
    return a_max_dist_micron_by_f, P_times_omega_by_f, a_max_best_by_f, POLF_obs_by_f, best_sf_by_f, best_idx_by_f, chi_sq_by_f
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
def load_Pomega_from_data(data,
                          wavelengths,
                          f):
    """
    Extract the a_max values and P*omega matrix from the processed data.

    Parameters
    ----------
    data : dict
        Dictionary indexed as data[wavelength][f].
    wavelengths : list
        Wavelengths in the order corresponding to the ALMA bands.
    f : float
        Filling factor.

    Returns
    -------
    a_max_vals : ndarray
        Particle sizes in microns.
    Pw : ndarray
        Shape (N_a_max, N_bands). Model P*omega values.
    """

    # Same a array for every wavelength
    a_max_vals = data[wavelengths[0]][str(f)]["a"] * 1e4 # 'a' are in cm

    Pw = np.column_stack([
        data[wav][str(f)]["P_omega"]
        for wav in wavelengths
    ])

    return a_max_vals, Pw

# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# I am updating this so we can say what POLF we want to use 
def find_sf_v5(bands,
               bands_naming,
               bands_included_in_fit,
               data,
               f_values,
               wavelengths,
               df_POLF,
               POLF_index,
               print_results=False,
               trouble = False,
               trouble_f = 0.75,
               choose_trouble_index = False,
               trouble_index = 2):
    
    polf_columns = {
        'gaussian': 'POLF_Gaussian',
        'max Stokes I': 'POLF_maxStokesI',
        'POLI': 'POLF_maxPOLI',
        'mean': 'POLF_mean'
    }
    
    polf_err_columns = {
    'gaussian': 'POLF_err_maxStokesI',
    'max Stokes I': 'POLF_err_maxStokesI',
    'POLI': 'POLF_err_maxPOLI',
    #'mean': 'POLF_err_mean'   # if you have this column
    }
    
    
    if POLF_index not in polf_columns:
        raise ValueError(
            "POLF_index must be 'gaussian', 'poli', 'avg', 'mean', or 'max Stokes I'"
        )
        
        
    
    print(rf'The Bands we are looking at are: {bands}')
    print(rf"The Bands included in the fit are: {bands_included_in_fit}")
 
    # Band names → indices (0, 1, 2, etc)
    band_to_index = {band: i for i, band in enumerate(bands)}
    fit_indices = [band_to_index[b] for b in bands_included_in_fit]
    
    
    # This stores the data for each f value
    a_max_dist_micron_by_f = {}
    P_times_omega_by_f = {}
    a_max_best_by_f = {}
    POLF_obs_by_f = {}
    best_sf_by_f = {}
    best_idx_by_f = {}
    chi_sq_by_f = {}



    POLF_obs_all = df_POLF[polf_columns[POLF_index]].values
    POLF_err_obs_all = df_POLF[polf_err_columns[POLF_index]].values

    

    # LOOP OVER FILES (f values)
    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    for f in f_values:
        a_max_vals, Pw = load_Pomega_from_data(
            data,
            wavelengths,
            f,
        )
        

        sf_medians, sf_stds = calculate_scale_factors(
                                a_max_vals,
                                Pw,
                                POLF_obs_all,
                                fit_indices,
                            )
            

        
        
        # Choose best a_max where scatter is minimized
        best_idx = np.argmin(sf_stds)
        


        # Fitting Band 4 and 7 only at f = 0.75 gave me a weird result
        # Want to look at it closer
        # ----------------------------------------
        if trouble:
            if f in trouble_f:
                sorted_indices = np.argsort(sf_stds)

                print(f'trouble f = {trouble_f}')
                # Print results 
                print(rf'The best a_max is: {a_max_vals[best_idx]}')
                print(rf'At this value, the best sf is : {sf_medians[best_idx]}')
                print(' ')
                
                # Get the next two lwest
                second_best_idx = sorted_indices[1]
                third_best_idx = sorted_indices[2]
                fourth_best_idx = sorted_indices[3]

                
                # Print results 
                print(rf'The second best a_max is: {a_max_vals[second_best_idx]}')
                print(rf'At this value, the best sf is : {sf_medians[second_best_idx]}')
                print(' ')
                
                # Print results 
                print(rf'The third best a_max is: {a_max_vals[third_best_idx]}')
                print(rf'At this value, the best sf is : {sf_medians[third_best_idx]}')
                print(' ')
                
                # Print results 
                print(rf'The fourth best a_max is: {a_max_vals[fourth_best_idx]}')
                print(rf'At this value, the best sf is : {sf_medians[fourth_best_idx]}')
                print(' ')
                print(' ')
                
                
                if choose_trouble_index:
                    best_idx = sorted_indices[trouble_index]
                    
            
                    
        best_sf = sf_medians[best_idx]
                    
        # ----------------------------------------
        # Calculate chi^2 for each f_value
        # ----------------------------------------
        chi_sq = calculate_chi_squared_for_sf(
            Pw,
            best_idx,
            best_sf,
            POLF_obs_all,
            POLF_err_obs_all,
            bands,
        )
        # ---------------------------------------------
        # ---------------------------------------------
        # ----------------------------------------
                    
        # ----------------------------------------
        
        # Save the information to each by f array
        a_max_dist_micron_by_f[f] = a_max_vals
        P_times_omega_by_f[f] = Pw
        a_max_best_by_f[f] = a_max_vals[best_idx]
        best_sf_by_f[f] = best_sf
        best_idx_by_f[f] = best_idx
        POLF_obs_by_f[f] = np.array([POLF_obs_all[band_to_index[b]] for b in bands])
        chi_sq_by_f[f] = chi_sq
    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    # End of looping over all f values
        

    df_best = pd.DataFrame({
        'f': list(best_idx_by_f.keys()),
        'Best idx': list(best_idx_by_f.values()),
        'Best a_max': list(a_max_best_by_f.values()),
        'Best SF': list(best_sf_by_f.values()),
        'Best chi^2': list(chi_sq_by_f.values())
    })

    df_best['Best a_max'] = df_best['Best a_max'].round(0)
    df_best['Best SF'] = df_best['Best SF'].round(2)
    df_best['Best chi^2'] = df_best['Best chi^2'].round(2)

#         df_best

    # return values
    return a_max_dist_micron_by_f, P_times_omega_by_f, a_max_best_by_f, POLF_obs_by_f, best_sf_by_f, best_idx_by_f, chi_sq_by_f, df_best
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------








# # I am updating this so we can say what POLF we want to use 
# def find_sf_v4(bands,
#                bands_naming,
#                bands_included_in_fit,
#                csv_files,
#                P_omega_data_folder_path,
#                df_POLF,
#                POLF_index,
#                print_results=False,
#                trouble = False,
#                trouble_f = 0.75,
#                choose_trouble_index = False,
#                trouble_index = 2):
    
#     polf_columns = {
#         'gaussian': 'POLF_Gaussian',
#         'max Stokes I': 'POLF_maxStokesI',
#         'POLI': 'POLF_maxPOLI',
#         'mean': 'POLF_mean'
#     }
    
    
#     if POLF_index not in polf_columns:
#         raise ValueError(
#             "POLF_index must be 'gaussian', 'poli', 'avg', 'mean', or 'max Stokes I'"
#         )
        
        
    
#     print(rf'The Bands we are looking at are: {bands}')
#     print(rf"The Bands included in the fit are: {bands_included_in_fit}")
 
#     # Band names → indices (0, 1, 2, etc)
#     band_to_index = {band: i for i, band in enumerate(bands)}
#     fit_indices = [band_to_index[b] for b in bands_included_in_fit]
    
    
#     # This stores the data for each f value
#     a_max_dist_micron_by_f = {}
#     P_times_omega_by_f = {}
#     a_max_best_by_f = {}
#     POLF_obs_by_f = {}
#     best_sf_by_f = {}
#     best_idx_by_f = {}
#     chi_sq_by_f = {}



#     POLF_obs_all = df_POLF[polf_columns[POLF_index]].values

    

#     # LOOP OVER FILES (f values)
#     # ------------------------------------------------------------------------
#     # ------------------------------------------------------------------------
#     for f, fname in csv_files.items():

#         # Read in the df
#         df = pd.read_csv(P_omega_data_folder_path + fname)

#         # Extract the x-axis values
#         a_max_vals = df["a_max_micron"].values

#         # ALL bands
#         Pw = np.column_stack([
#             df[f"P_times_omega_{b}"].values for b in bands_naming
#         ])
        

#         # Make arrays to store the median and standard deviation values 
#         sf_medians = []
#         sf_stds = []
    
    

#         # Loop over each a_max grid point
#         for i in range(len(a_max_vals)):

#             Pw_model = Pw[i, :]  # all bands

#              # ONLY use selected bands for fit
#             Pw_fit = Pw_model[fit_indices]
#             POLF_fit = POLF_obs_all[fit_indices]

#             # Avoid divide by zero
#             valid = Pw_fit > 0

#             sf_i = POLF_fit[valid] / Pw_fit[valid]

#             sf_medians.append(np.median(sf_i))
#             sf_stds.append(np.std(sf_i))
            
# #             print(sf_i)

#         sf_medians = np.array(sf_medians)
#         sf_stds = np.array(sf_stds)

        
        
#         # Choose best a_max where scatter is minimized
#         best_idx = np.argmin(sf_stds)
        


#         # Fitting Band 4 and 7 only at f = 0.75 gave me a weird result
#         # Want to look at it closer
#         # ----------------------------------------
#         if trouble:
#             if f == trouble_f:
#                 sorted_indices = np.argsort(sf_stds)

                
#                 # Print results 
#                 print(rf'The best a_max is: {a_max_vals[best_idx]}')
#                 print(rf'At this value, the best sf is : {sf_medians[best_idx]}')
#                 print(' ')
                
#                 # Get the next two lwest
#                 second_best_idx = sorted_indices[1]
#                 third_best_idx = sorted_indices[2]

                
#                 # Print results 
#                 print(rf'The second best a_max is: {a_max_vals[second_best_idx]}')
#                 print(rf'At this value, the best sf is : {sf_medians[second_best_idx]}')
#                 print(' ')
                
#                 # Print results 
#                 print(rf'The third best a_max is: {a_max_vals[third_best_idx]}')
#                 print(rf'At this value, the best sf is : {sf_medians[third_best_idx]}')
#                 print(' ')
                
                
#                 if choose_trouble_index:
#                     best_idx = sorted_indices[trouble_index]
                    
#         best_sf = sf_medians[best_idx]
                    
#         # ----------------------------------------
#         # Calculate chi^2 for each f_value
#         # ----------------------------------------
#         chi_sq_running_sum = 0

#         for j, b in enumerate(bands):

#             if b == "Band 5":
#                 continue

#             band_idx = bands.index(b)

#             POLF_obs_i = POLF_obs_all[band_idx]
#             POLF_model_i = Pw[best_idx, band_idx] * best_sf

#             chi_sq_running_sum += (POLF_obs_i - POLF_model_i)**2

#         chi_sq = chi_sq_running_sum
#         # ---------------------------------------------
#         # ---------------------------------------------
#         # ----------------------------------------
                    
#         # ----------------------------------------
        
#         # Save the information to each by f array
#         a_max_dist_micron_by_f[f] = a_max_vals
#         P_times_omega_by_f[f] = Pw
#         a_max_best_by_f[f] = a_max_vals[best_idx]
#         best_sf_by_f[f] = best_sf
#         best_idx_by_f[f] = best_idx
#         POLF_obs_by_f[f] = np.array([POLF_obs_all[band_to_index[b]] for b in bands])
#         chi_sq_by_f[f] = chi_sq
#     # ------------------------------------------------------------------------
#     # ------------------------------------------------------------------------
#     # End of looping over all f values
        

#     if print_results:
#         print("Best idx:", best_idx_by_f)
#         print("Best a_max:", a_max_best_by_f)
#         print("Best SF:", best_sf_by_f)
#         print("Best chi^2:", chi_sq_by_f)

#     # return values
#     return a_max_dist_micron_by_f, P_times_omega_by_f, a_max_best_by_f, POLF_obs_by_f, best_sf_by_f, best_idx_by_f, chi_sq_by_f