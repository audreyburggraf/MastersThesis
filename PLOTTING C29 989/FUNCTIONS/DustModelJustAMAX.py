import numpy as np 
from scipy.interpolate import interp1d
import pandas as pd

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
        'POLI': 'POLF_maxPOLI',
        'mean': 'POLF_mean'
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

    

    # LOOP OVER FILES (f values)
    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    for f, fname in csv_files.items():

        # Read in the df
        df = pd.read_csv(P_omega_data_folder_path + fname)

        # Extract the x-axis values
        a_max_vals = df["a_max_micron"].values

        # ALL bands
        Pw = np.column_stack([
            df[f"P_times_omega_{b}"].values for b in bands_naming
        ])
        

        # Make arrays to store the median and standard deviation values 
        sf_medians = []
        sf_stds = []
    
    

        # Loop over each a_max grid point
        for i in range(len(a_max_vals)):

            Pw_model = Pw[i, :]  # all bands

             # ONLY use selected bands for fit
            Pw_fit = Pw_model[fit_indices]
            POLF_fit = POLF_obs_all[fit_indices]

            # Avoid divide by zero
            valid = Pw_fit > 0

            sf_i = POLF_fit[valid] / Pw_fit[valid]

            sf_medians.append(np.median(sf_i))
            sf_stds.append(np.std(sf_i))
            
#             print(sf_i)

        sf_medians = np.array(sf_medians)
        sf_stds = np.array(sf_stds)

        
        
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
        chi_sq_running_sum = 0

        for j, b in enumerate(bands):

            if b == "Band 5":
                continue

            band_idx = bands.index(b)

            POLF_obs_i = POLF_obs_all[band_idx]
            POLF_model_i = Pw[best_idx, band_idx] * best_sf

            chi_sq_running_sum += (POLF_obs_i - POLF_model_i)**2

        chi_sq = chi_sq_running_sum
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