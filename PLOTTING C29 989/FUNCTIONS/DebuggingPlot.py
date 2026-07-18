import numpy as np 
from scipy.interpolate import interp1d
import pandas as pd


from UnitConversion import * 
from DSHARP_Functions import * 



# Import necessary packages
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns
import math
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Import Functions
from FITS_Image_Functions import *
from PolarizationFunctions import *
from DataAnalysisFunctions import *
from GaussianFunctions import *
from SlicesFunctions import * 

# Import the constants
# -----------------------------------------------------------------------------------------
import sys

sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/") 

import constants

normalized_cbar_ticks = constants.normalized_cbar_ticks

title_fs = constants.title_fs
axis_label_fs = constants.axis_label_fs
axis_num_fs = constants.axis_num_fs
legend_title_fs = constants.legend_title_fs
legend_text_fs = constants.legend_text_fs
cbar_fs = constants.cbar_fs
text_fs = constants.text_fs

alma_band_colors = constants.alma_band_colors
alma_band_ms = constants.alma_band_ms





def MakeCurves(DebuggingVersion, f_values, min_bounds_micron, max_bounds_micron, N_grains, lambda_bands_cm, bands_naming, zhang = False):
    
    
    # ---------------------------------------------------------------------
    # a_max values for this plot
    a_max_test_micron = [1, 100]
    a_max_test_cm = micron_to_cm(a_max_test_micron)
    # ---------------------------------------------------------------------
    # Set up wavelength array
    logwave_vals = np.linspace(0.1, 4, 10000)

    lambda_dist_micron = 10**logwave_vals

    lambda_dist_cm = micron_to_cm(lambda_dist_micron)
    # ---------------------------------------------------------------------
    
    base_path =  "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/DUST MODEL NOTEBOOKS/P_omega_Data/"
    
    if DebuggingVersion == 'v1':
        path = base_path + "Debugging_v1/"
    
    elif DebuggingVersion == 'v2':
        path = base_path + "Debugging_v2/"     
    
    elif DebuggingVersion == 'v3':
        path = base_path + "Debugging_v3/"
            
    elif DebuggingVersion == 'v4':
        path = base_path + "Debugging_v4/"
            
    elif DebuggingVersion == 'v5':
        path = base_path + "Debugging_v5/"
    elif DebuggingVersion == 'v5 poster':
        path = base_path + "Poster/"

    # Loop over the f_values
    # ------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------
    for f in f_values:
#         debug_label = None
#         if f == 0.5:
#             debug_label = f"{DebuggingVersion}, f={f}"
    
#         print(rf"Now working on f = {f}")
        
        
        # Choose how the curves are made based on the version we are debugging
        # ------------------------------
        # In this version (I think) the data is incorrectly made to be a/f, and then fed into DSHARP
        if DebuggingVersion == 'v1':
            
            Xdata_cm = micron_to_cm((1/f) * np.linspace(min_bounds_micron[f],
                                                        max_bounds_micron[f],
                                                        N_grains[f]))
            
            a_min_micron = cm_to_micron(min(Xdata_cm))
            a_max_micron = cm_to_micron(max(Xdata_cm))
            print(rf'f = {f}, a_min/f = {a_min_micron:.3f}, a_max/f = {a_max_micron:.3f}')


            print(rf'f = {f}, a_min = {a_min_micron * f:.3f}, a_max = {a_max_micron * f:.3f}')
            print(' ')
            
            # Run DSHARP for this f
            P, omega, _, P_times_omega, _ = run_DSHARP(f, a_max_test_cm, Xdata_cm, lambda_bands_cm, lambda_dist_cm)
        # ------------------------------
        # In this version (I think) the data is incorrectly made to be a/f, and then fed into DSHARP correctly
        elif DebuggingVersion == 'v2':
            
            Xdata_cm = micron_to_cm((1/f) * np.linspace(min_bounds_micron[f],
                                                        max_bounds_micron[f],
                                                        N_grains[f]))
            
            a_min_micron = cm_to_micron(min(Xdata_cm))
            a_max_micron = cm_to_micron(max(Xdata_cm))
            print(rf'f = {f}, a_min/f = {a_min_micron:.3f}, a_max/f = {a_max_micron:.3f}')


            print(rf'f = {f}, a_min = {a_min_micron * f:.3f}, a_max = {a_max_micron * f:.3f}')
            print(' ')
            
            # Run DSHARP for this f
            P, omega, _, P_times_omega, _ = run_DSHARP(f, a_max_test_cm, Xdata_cm * f, lambda_bands_cm, lambda_dist_cm)
        # ------------------------------
        if DebuggingVersion == 'v3':
            
            Xdata_cm = micron_to_cm(f * np.linspace(min_bounds_micron[f],
                                                        max_bounds_micron[f],
                                                        N_grains[f]))
            
            
    
            # Run DSHARP for this f
            P, omega, _, P_times_omega, _ = run_DSHARP(f, a_max_test_cm, Xdata_cm, lambda_bands_cm, lambda_dist_cm)
        # ------------------------------
        # In this version (I think) the data is incorrectly made to be a/f, and then fed into DSHARP correctly
        elif DebuggingVersion == 'v4':
            
            Xdata_cm = micron_to_cm(f * np.linspace(min_bounds_micron[f],
                                                        max_bounds_micron[f],
                                                        N_grains[f]))
            
            print(f"f = {f}")
            print("DSHARP min a_max =", cm_to_micron(min(Xdata_cm / f)))
            print("DSHARP max a_max =", cm_to_micron(max(Xdata_cm / f)))
            print()
        
            # Run DSHARP for this f
            P, omega, _, P_times_omega, _ = run_DSHARP(f, a_max_test_cm, Xdata_cm /f, lambda_bands_cm, lambda_dist_cm)
        # ------------------------------
        # In this version (I think) the data is incorrectly made to be a/f, and then fed into DSHARP correctly
        elif DebuggingVersion in ['v5', 'v5 poster']:
            
            Xdata_cm = micron_to_cm(np.linspace(min_bounds_micron[f],
                                                max_bounds_micron[f],
                                                N_grains[f]))
          
            
            # Run DSHARP for this f
            P, omega, _, P_times_omega, _ = run_DSHARP(f, a_max_test_cm, Xdata_cm, lambda_bands_cm, lambda_dist_cm)
        # ------------------------------


        # Run DSHARP for this f
        #P, omega, P_times_omega = run_DSHARP(f, a_max_test_cm, Xdata_cm * f, lambda_bands_cm, lambda_dist_cm)

        # Prepare data dictionary
        data = {
            "Xdata": Xdata_cm * 1e4
        }

        # Add per-band columns
        for i, b in enumerate(bands_naming):  # bands = [4, 5, 6, ...]
            data[f"P_{b}"] = P[:, i]
            data[f"omega_{b}"] = omega[:, i]
            data[f"P_times_omega_{b}"] = P_times_omega[:, i]

        # Convert to DataFrame
        df = pd.DataFrame(data)

        # Save to CSV with f in the filename
        f_str = str(f).replace('.', '_')  # e.g., 0.3 → '0_3'

        if DebuggingVersion == 'v5 poster':
            file_name = f"data_v5_poster_{f_str}_JUST_AMAX_POSTER.csv"

        elif zhang:
            file_name = f"data_debugging_{DebuggingVersion}_{f_str}_ZHANG.csv"

        else:
            file_name = f"data_debugging_{DebuggingVersion}_{f_str}.csv"
            

        df.to_csv(path + file_name, index=False)

        
        print(f"Saved {file_name} to {path + file_name}")
   # ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------





def MakeCurvesZHANG(DebuggingVersion, f_values, min_bounds_micron, max_bounds_micron, N_grains, lambda_bands_cm, bands_naming, zhang = False):
    
    
    # ---------------------------------------------------------------------
    # a_max values for this plot
    a_max_test_micron = [1, 100]
    a_max_test_cm = micron_to_cm(a_max_test_micron)
    # ---------------------------------------------------------------------
    # Set up wavelength array
    logwave_vals = np.linspace(0.1, 4, 10000)

    lambda_dist_micron = 10**logwave_vals

    lambda_dist_cm = micron_to_cm(lambda_dist_micron)
    # ---------------------------------------------------------------------
    
    base_path =  "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/DUST MODEL NOTEBOOKS/P_omega_Data/"
    
    if DebuggingVersion == 'v1':
        path = base_path + "Debugging_v1/"
    
    elif DebuggingVersion == 'v2':
        path = base_path + "Debugging_v2/"     
    
    elif DebuggingVersion == 'v3':
        path = base_path + "Debugging_v3/"
            
    elif DebuggingVersion == 'v4':
        path = base_path + "Debugging_v4/"
            
    elif DebuggingVersion == 'v5':
        path = base_path + "Debugging_v5/"

    # Loop over the f_values
    # ------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------
    for f in f_values:
#         debug_label = None
#         if f == 0.5:
#             debug_label = f"{DebuggingVersion}, f={f}"
    
#         print(rf"Now working on f = {f}")
        
        
        # Choose how the curves are made based on the version we are debugging
        # ------------------------------
        # In this version (I think) the data is incorrectly made to be a/f, and then fed into DSHARP
        if DebuggingVersion == 'v1':
            
            Xdata_cm = micron_to_cm((1/f) * np.logspace(min_bounds_micron[f],
                                                        max_bounds_micron[f],
                                                        N_grains[f]))
            
            a_min_micron = cm_to_micron(min(Xdata_cm))
            a_max_micron = cm_to_micron(max(Xdata_cm))
            print(rf'f = {f}, a_min/f = {a_min_micron:.3f}, a_max/f = {a_max_micron:.3f}')


            print(rf'f = {f}, a_min = {a_min_micron * f:.3f}, a_max = {a_max_micron * f:.3f}')
            print(' ')
            
            # Run DSHARP for this f
            P, omega, _, P_times_omega, _ = run_DSHARP(f, a_max_test_cm, Xdata_cm, lambda_bands_cm, lambda_dist_cm)
        # ------------------------------
        # In this version (I think) the data is incorrectly made to be a/f, and then fed into DSHARP correctly
        elif DebuggingVersion == 'v2':
            
            Xdata_cm = micron_to_cm((1/f) * np.logspace(min_bounds_micron[f],
                                                        max_bounds_micron[f],
                                                        N_grains[f]))
            
            a_min_micron = cm_to_micron(min(Xdata_cm))
            a_max_micron = cm_to_micron(max(Xdata_cm))
            print(rf'f = {f}, a_min/f = {a_min_micron:.3f}, a_max/f = {a_max_micron:.3f}')


            print(rf'f = {f}, a_min = {a_min_micron * f:.3f}, a_max = {a_max_micron * f:.3f}')
            print(' ')
            
            # Run DSHARP for this f
            P, omega, _, P_times_omega, _ = run_DSHARP(f, a_max_test_cm, Xdata_cm * f, lambda_bands_cm, lambda_dist_cm)
        # ------------------------------
        if DebuggingVersion == 'v3':
            
            Xdata_cm = micron_to_cm(f * np.logspace(min_bounds_micron[f],
                                                        max_bounds_micron[f],
                                                        N_grains[f]))
            
            
            # Run DSHARP for this f
            P, omega, _, P_times_omega, _ = run_DSHARP(f, a_max_test_cm, Xdata_cm, lambda_bands_cm, lambda_dist_cm)
        # ------------------------------
        # In this version (I think) the data is incorrectly made to be a/f, and then fed into DSHARP correctly
        elif DebuggingVersion == 'v4':
            
            Xdata_cm = micron_to_cm(f * np.logspace(min_bounds_micron[f],
                                                        max_bounds_micron[f],
                                                        N_grains[f]))
            
            # Run DSHARP for this f
            P, omega, omega_eff, P_times_omega, P_times_omega_eff= run_DSHARP(f, a_max_test_cm, Xdata_cm /f, lambda_bands_cm, lambda_dist_cm)
        # ------------------------------
        # In this version (I think) the data is incorrectly made to be a/f, and then fed into DSHARP correctly
        elif DebuggingVersion == 'v5':
            
            Xdata_cm = micron_to_cm(np.logspace(min_bounds_micron[f],
                                                max_bounds_micron[f],
                                                N_grains[f]))
          
            
            # Run DSHARP for this f
            P, omega, _, P_times_omega, _ = run_DSHARP(f, a_max_test_cm, Xdata_cm, lambda_bands_cm, lambda_dist_cm)
        # ------------------------------


        # Run DSHARP for this f
        #P, omega, P_times_omega = run_DSHARP(f, a_max_test_cm, Xdata_cm * f, lambda_bands_cm, lambda_dist_cm)

        # Prepare data dictionary
        data = {
            "Xdata": Xdata_cm * 1e4
        }

        # Add per-band columns
        for i, b in enumerate(bands_naming):  # bands = [4, 5, 6, ...]
            data[f"P_{b}"] = P[:, i]
            data[f"omega_{b}"] = omega[:, i]
            data[f"P_times_omega_{b}"] = P_times_omega[:, i]
            data[f"omega_eff_{b}"] = omega_eff[:, i]
            data[f"P_times_omega_eff_{b}"] = P_times_omega_eff[:, i]
            

        # Convert to DataFrame
        df = pd.DataFrame(data)

        # Save to CSV with f in the filename
        f_str = str(f).replace('.', '_')  # e.g., 0.3 → '0_3'

        if zhang:
            file_name = f"data_debugging_{DebuggingVersion}_{f_str}_ZHANG.csv"
        else:
            file_name = f"data_debugging_{DebuggingVersion}_{f_str}.csv"

        df.to_csv(path + file_name, index=False)

        
        print(f"Saved {file_name} to {path + file_name}")
   # ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------








def STDEV_Calculations(Xdata, Pw, fit_indices, POLF_obs_all):
    
    # Make arrays to store the median and standard deviation values 
    sf_medians = []
    sf_stds = []

    # Loop over each a_max grid point
    for i in range(len(Xdata)):

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

    return sf_medians, best_idx










# --------------------------------------------------------------------------------------
def DebuggingPlot(DebuggingVersion, 
                              bands,        # These are the bands we are looking at
                              bands_labels, 
                              bands_included_in_fit, # These are the values that were fit
                              f_values,     # These are the porosity values we are testing
                              Xdata, 
                              P_times_omega,
                              sf,           # This is the values the data is being scaled by 
                              best_idx, 
                              Xdata_best,        # These are the best a_max values for each f
                              POLF_markers,    # These are the values for each band we are plotting the marker at
                              x_axis, 
                              plot_sf = 1.5, # This controls how much the text labels are scaled by
                              ymin = -0.1, 
                              ymax = 1.7, 
#                               hlines = False,
                              custom_lw = None,
                              custom_text_x = None,
                             chi_sq_precision = 3,
                            
                             for_poster = False):
    
    
    
    band_ls = [
    '-' if b in bands_included_in_fit else '--'
    for b in bands
    ]
    
    if custom_lw is None:
        # default behavior
        bands_lw = [
            5 if b in bands_included_in_fit else 3
            for b in bands
        ]
    else:
        # override where specified
        bands_lw = [
            custom_lw.get(b, 5 if b in bands_included_in_fit else 3)
            for b in bands
        ]
    
    # Set the colors and marker size for each band 
    band_colors = [alma_band_colors[b] for b in bands]
    ms = [constants.alma_band_ms[b] for b in bands]
    
    
    
    # Make the figure
    # Later I will make this so you can customize how many figures you need
    fig, ax = plt.subplots(1, 4, figsize=(30, 8))
    
    
    
    
    # Loop over the f values you want to plot
    for i, f in enumerate(f_values): 
        
        
        x_pos = 0.05 
        
        if custom_text_x is not None:
            x_pos = x_pos + custom_text_x[i]
#         print(rf'f = {f}')
        
       

        if i == 0:
            ax[i].set_ylabel('P $\omega$', fontsize=axis_label_fs * plot_sf)
        else:
            ax[i].set_yticklabels([])
            ax[i].set_ylabel('')
            ax[i].tick_params(left=False)

        # Make the x-axis log scale 
        ax[i].set_xscale("log")


        # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        if DebuggingVersion in ['v1', 'v2']:
            
            
           
                    # Add the a_max to the text
            ax[i].text(x_pos, 0.05, rf'$a_{{\mathrm{{max}}}}$ = {Xdata_best[f] * f:.0f} $\mu$m', 
                       transform=ax[i].transAxes, 
                       fontsize = 30)
            ax[i].text(x_pos, 0.90, rf'$a_{{\mathrm{{max}}}}/f$ = {Xdata_best[f]:.0f} $\mu$m', 
           transform=ax[i].transAxes, 
           fontsize = 30)
             # ---------------------------------------------------------------------------------   
            if x_axis == 'a/f':
                xloc_marker = Xdata_best[f] 
                ax[i].set_xlabel(r'$a_{\mathrm{max}} / f$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = Xdata[f] # this is a/f
                
            elif x_axis == 'a':
                xloc_marker = Xdata_best[f] 
                ax[i].set_xlabel(r'$a_{\mathrm{max}}$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = Xdata[f] * f
              # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        if DebuggingVersion in ['v3', 'v4']:
            
            
            
                    # Add the a_max to the text
            ax[i].text(x_pos, 0.05, rf'$a_{{\mathrm{{max}}}}$ = {Xdata_best[f] / f:.0f} $\mu$m', 
                       transform=ax[i].transAxes, 
                       fontsize = 30)
            ax[i].text(x_pos, 0.90, rf'$a_{{\mathrm{{max}}}}*f$ = {Xdata_best[f]:.0f} $\mu$m', 
           transform=ax[i].transAxes, 
           fontsize = 30)
             # ---------------------------------------------------------------------------------   
            if x_axis == 'a*f':
                xloc_marker = Xdata_best[f]
                ax[i].set_xlabel(r'$a_{\mathrm{max}} * f$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = Xdata[f] # this is a/f
                
            elif x_axis == 'a':
                xloc_marker = Xdata_best[f] /f
                ax[i].set_xlabel(r'$a_{\mathrm{max}}$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = Xdata[f] /f
              # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        if DebuggingVersion == 'v5':
            
            
            xloc_marker = Xdata_best[f] 
                    # Add the a_max to the text
            ax[i].text(x_pos, 0.05, rf'$a_{{\mathrm{{max}}}}$ = {Xdata_best[f]:.0f} $\mu$m', 
                       transform=ax[i].transAxes, 
                       fontsize = 30)
             # ---------------------------------------------------------------------------------   

                
            if x_axis == 'a':
                ax[i].set_xlabel(r'$a_{\mathrm{max}}$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = Xdata[f] 
              # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------


    
    
    
    

    
        # Set ticks
        ax[i].minorticks_on()
        ax[i].tick_params(axis="x", which="minor", direction="in", left=True, right=True, length=4)
        ax[i].tick_params(axis="y", which="minor", direction="in", left=True, right=True, length=4)
        ax[i].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs*plot_sf)
        ax[i].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs*plot_sf)

        
        
        # ---------------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------------

        # f labels on plot
        # ---------------------------------------------------------------------------------------
        if for_poster == False:
            ax[i].set_title(f'$f$ = {f:.2f}', fontsize=50)
        
        if for_poster == True:
            ax[i].text(x_pos, 0.15,  f'$\chi^2$ = {chi_sq:.{chi_sq_precision}f}', transform=ax[i].transAxes, fontsize = 30)
        # ---------------------------------------------------------------------------------------


        # In the loop over each plot
        # ------------------------------------------
#         if hlines:
#             ax[i].axhline(1.367706, color = alma_band_colors[4], ls = '--')
#             ax[i].axhline(0.990896, color = alma_band_colors[5], ls = '--')
#             ax[i].axhline(1.463922, color = alma_band_colors[6], ls = '--')
#             ax[i].axhline(1.042140, color = alma_band_colors[7], ls = '--')
        
        ax[i].set_ylim(ymin, ymax)
        # ------------------------------------------
        
        
        marker_handles = []
        # plot each band for this f
        for j, b in enumerate(bands):
            
            idx = bands.index(b)
            
            
#             print(rf'At band {b}, the POLF_obs value is: {POLF_obs[f][idx]}')
#             print(' ')

            fc = band_colors[j]
            ec = band_colors[j]
            marker_lw = 0 
    
            if b in ["Band 5"]:
                fc = 'none'
                marker_lw = 5

            # ---------------------------------------------------------------------
            # Plot the POLF markers
            ax[i].scatter(xloc_marker, 
                          POLF_markers[b],
                          marker=ms[j], 
                          s=400,
                          facecolors=fc,
                          edgecolors=ec,
                          linewidths=marker_lw)
            # ---------------------------------------------------------------------

        
            marker_handles.append(
            Line2D([0], [0],
                   marker=ms[j],
                   color='none',
                   markerfacecolor=fc,
                   markeredgecolor=ec,
                   markeredgewidth=marker_lw,
                   markersize=15,
                   linestyle='None'))
            
            
            if b not in ["Band 5 robust -1", "Band 5 robust -2"]:

                ax[i].plot(x_data_plot,
                           P_times_omega[f][:, idx] * sf[f],
                           color=band_colors[j],
                           label=f'{bands_labels[j]}',
                           lw=bands_lw[j],
                           ls=band_ls[j])

                

        # Calculate chi_sq
        # ---------------------------------------------
        # ---------------------------------------------
        chi_sq_running_sum = 0
        chi_sq = 0
        
        #print('In chi squared calculation:')
        #print(rf'f = {f}')
        #print(rf'best_idx = {best_idx[f]}')
#         print(rf'The best_index is {best_idx[f]}')
        # Loop over each band 
        for j, b in enumerate(bands):
            
            if b == "Band 5":
                continue
            
            band_idx = bands.index(b)
            
            
            POLF_obs = POLF_markers[b]
            POLF_model = P_times_omega[f][best_idx[f], band_idx] * sf[f]
            
            chi_sq_running_sum = chi_sq_running_sum + (POLF_obs - POLF_model)**2
            
            
            #print(rf'At band {bf}, or index {band_idx}:')
            #print(rf'   the observed POLF is {POLF_obs:.2f}, and the measured POLF is {POLF_model:.2f} and the difference = {(POLF_obs - POLF_model):.2f}')
        chi_sq = chi_sq_running_sum
        #print(' ')
        # ---------------------------------------------
        # ---------------------------------------------
        
        
        
                
       
        
        
        
        ax[i].text(x_pos, 0.15,  f'$\chi^2$ = {chi_sq:.{chi_sq_precision}f}', transform=ax[i].transAxes, fontsize = 30)

        # Only add sf if for_poster = False
        if for_poster == False:
            ax[i].text(x_pos, 0.25, f'sf = {sf[f]:.2f}', transform=ax[i].transAxes, fontsize = 30)


        
        
        
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.01, hspace=0)
    
    
    # Add coloured text for each Band in the far right plot
    # ----------------------------------------------------
    # ----------------------------------------------------
    start = 0.9 # This is the starting vertical position
    c = 0
    marker_legend_handles = []
    
    # Loop over each band 
    for j, b_label in enumerate(bands_labels):
        
        if "robust" in b_label:
            label = b_label.replace(" robust ", "\nrob. ")
            c = c + 0.075 # This controls the vertical spacing 
        else:
            label = b_label
            
            
        handle = Line2D(
            [0], [0],
            marker = ms[j],  
            color='none',
            markerfacecolor=band_colors[j],
            markeredgecolor=band_colors[j],
            markersize=20,
            linestyle='None',
            label=label)

        marker_legend_handles.append(handle)

#         ax[3].text(0.75, # This controls the horizontal position
#                    start - c, 
#                    label,
#                    transform=ax[3].transAxes,
#                    fontsize=30,
#                    color=band_colors[j])
        c = c + 0.1
        
        
        
    legend_elements = [
        Line2D([0], [0], color='black', linestyle='-', lw = 5, label='Included'),
        Line2D([0], [0], color='black', linestyle='--', lw = 5, label='Excluded')
    ]

    # Line style
#     ax[2].legend(handles=legend_elements, fontsize = 30, frameon=False, title="Fit", title_fontsize = 30)   
    
    # ----------------------------------------------------
    # ----------------------------------------------------

    
    
    # Markers
#     ax[3].legend(handles=marker_handles, 
#                  fontsize = 30, 
#                  frameon=False, 
#                  # title="Obs. POLF", 
#                  title_fontsize = 30,
#                  loc = 'upper right',
#                  handletextpad=3.5,
#                  handlelength=1,  
#             )
    legend = ax[3].legend(
            handles=marker_legend_handles,
            fontsize=30,
            frameon=False,
            loc='upper right',
            handletextpad=0.5, # Cotrols distance between marker and text 
            handlelength=1, # Cotrols distance between marker and text 
            labelspacing=0.6 # Controls vertical gap between the markers
        )
    
    for text, color in zip(legend.get_texts(), band_colors):
        text.set_color(color)
#     # ----------------------------------------------------
#     # ----------------------------------------------------



# 

    return fig, ax
    

# --------------------------------------------------------------------------------------\










# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
def DebuggingPlotZHANG(DebuggingVersion, 
                       bands,        # These are the bands we are looking at
                       bands_labels, 
                             # bands_included_in_fit, # These are the values that were fit
                              f_values,     # These are the porosity values we are testing
                              Xdata, 
                              P_times_omega, P, omega, omega_eff, P_times_omega_eff, 
                              x_axis, 
                              plot_sf = 1.5, # This controls how much the text labels are scaled by
                              ymin = -0.1, 
                              ymax = 1.7, 
#                               hlines = False,
                              custom_lw = None,
                              custom_text_x = None,
                            fig_x = 30,
                       fig_y = 8,
                             for_poster = False,
                      band_colors = ['red', 'blue', 'black', 'magenta'],
                      bands_lw = 5,
                      sigma = 0.1, alpha = 0.75,
                      plot_omega = True, omega_type = 'omega_eff',
                      plot_P = True, P_type = 'P',
                      plot_Pw = True, 
                      Pw_type = 'P_omega_eff',
                      add_legend = True):
    
    
    
#     band_ls = [
#     '-' if b in bands_included_in_fit else '--'
#     for b in bands
#     ]
    
#     if custom_lw is None:
#         # default behavior
#         bands_lw = [
#             5 if b in bands_included_in_fit else 3
#             for b in bands
#         ]
#     else:
#         # override where specified
#         bands_lw = [
#             custom_lw.get(b, 5 if b in bands_included_in_fit else 3)
#             for b in bands
#         ]
    
    # Set the colors and marker size for each band 
    #band_colors = [alma_band_colors[b] for b in bands]
    #ms = [constants.alma_band_ms[b] for b in bands]
    
    
    
    # Make the figure
    # Later I will make this so you can customize how many figures you need
    fig, ax = plt.subplots(1, 4, figsize=(fig_x, fig_y))
    
    
    
    
    # Loop over the f values you want to plot
    for i, f in enumerate(f_values): 
        
        
        x_pos = 0.05 
        
        if custom_text_x is not None:
            x_pos = x_pos + custom_text_x[i]
#         print(rf'f = {f}')
        
       

        if i == 0:
            ax[i].set_ylabel('P $\omega$', fontsize=axis_label_fs * plot_sf)
        else:
            ax[i].set_yticklabels([])
            ax[i].set_ylabel('')
            ax[i].tick_params(left=False)

        # Make the x-axis log scale 
        if x_axis != 'log a*f':
            ax[i].set_xscale("log")


        # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        if DebuggingVersion in ['v1', 'v2']:
            
             # ---------------------------------------------------------------------------------   
            if x_axis == 'a/f':
                ax[i].set_xlabel(r'$a_{\mathrm{max}} / f um$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = Xdata[f] # this is a/f
                
            elif x_axis == 'a':
                ax[i].set_xlabel(r'$a_{\mathrm{max}} um$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = Xdata[f] * f
              # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        if DebuggingVersion in ['v3', 'v4']:
             # ---------------------------------------------------------------------------------   
            if x_axis == 'a*f':
                ax[i].set_xlabel(r'$a_{\mathrm{max}} * f/cm$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = micron_to_cm(Xdata[f]) # this is a/f
                
            elif x_axis == 'a':
                ax[i].set_xlabel(r'$a_{\mathrm{max}} um$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = Xdata[f] /f 
                
            if x_axis == 'log a*f':
                ax[i].set_xlabel(r'$\log_{10}(a_{\mathrm{max}}f/cm)$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = np.log10(micron_to_cm(Xdata[f])) 
              # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        if DebuggingVersion == 'v5':
             # ---------------------------------------------------------------------------------   

                
            if x_axis == 'a':
                ax[i].set_xlabel(r'$a_{\mathrm{max}} um$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = Xdata[f] 
            elif x_axis == 'a*f':
                ax[i].set_xlabel(r'$a_{\mathrm{max}} * f/cm$', fontsize=axis_label_fs*plot_sf)
                x_data_plot = micron_to_cm(Xdata[f] * f)
              # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------


    
    
    
    

    
        # Set ticks
        ax[i].minorticks_on()
        ax[i].tick_params(axis="x", which="minor", direction="in", left=True, right=True, length=4)
        ax[i].tick_params(axis="y", which="minor", direction="in", left=True, right=True, length=4)
        ax[i].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs*plot_sf)
        ax[i].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs*plot_sf)

        
        
        # ---------------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------------

        # f labels on plot
        # ---------------------------------------------------------------------------------------
        if for_poster == False:
            ax[i].set_title(f'$f$ = {f:.2f}', fontsize=50 * plot_sf)
        
        if for_poster == True:
            ax[i].text(x_pos, 0.15,  f'$\chi^2$ = {chi_sq:.{chi_sq_precision}f}', transform=ax[i].transAxes, fontsize = 30)
        # ---------------------------------------------------------------------------------------

        
        ax[i].set_ylim(ymin, ymax)
        # ------------------------------------------
        
        
        marker_handles = []
        # plot each band for this f
        for j, b in enumerate(bands):
            
            idx = bands.index(b)
            
            
#             print(rf'At band {b}, the POLF_obs value is: {POLF_obs[f][idx]}')
#             print(' ')

#             fc = band_colors[j]
#             ec = band_colors[j]
#             marker_lw = 0 
    
#             if b in ["Band 5"]:
#                 fc = 'none'
#                 marker_lw = 5

#             # ---------------------------------------------------------------------
#             # Plot the POLF markers
#             ax[i].scatter(xloc_marker, 
#                           POLF_markers[b],
#                           marker=ms[j], 
#                           s=400,
#                           facecolors=fc,
#                           edgecolors=ec,
#                           linewidths=marker_lw)
#             # ---------------------------------------------------------------------

        
#             marker_handles.append(
#             Line2D([0], [0],
#                    marker=ms[j],
#                    color='none',
#                    markerfacecolor=fc,
#                    markeredgecolor=ec,
#                    markeredgewidth=marker_lw,
#                    markersize=15,
#                    linestyle='None'))
            
    
            # Start the Pw, P, w legend
            style_legend = []
            
            # Plot P * omega
            # --------------------------------------
                    
            if plot_Pw:     
                
                if Pw_type == 'P_omega':
                    Pw_plot = P_times_omega[f][:, idx]
                    Pw_label = r'$P\omega$'
                elif Pw_type == 'P_omega_eff':
                    Pw_plot = P_times_omega_eff[f][:, idx]
                    Pw_label = r'$P\omega_{eff}$'
                else:
                    raise ValueError(f"Unknown Pw_type: {Pw_type}")
                
                style_legend.append(Line2D([0], [0], color='black', lw=2, ls='-', label=Pw_label))
                    
                ax[i].plot(x_data_plot,
                           Pw_plot,
                           color=band_colors[j],
                           label=f'{bands_labels[j]}',
                           lw=bands_lw,
                           ls= '-')
            # --------------------------------------
            # Plot P
            if plot_P:
               
  
                    
                style_legend.append(Line2D([0], [0], color='black', lw=2, ls='-.', label= '$P$'))
                    
                ax[i].plot(x_data_plot,
                           P[f][:, idx],
                           color=band_colors[j],
                           lw=bands_lw,
                           ls= '-.',
                           alpha = alpha)
            # --------------------------------------
            # Plot omega
            if plot_omega:

                if omega_type == 'omega_eff':
                    omega_plot = omega_eff[f][:, idx]
                    omega_label = r'$\omega_{eff}$'
                elif omega_type == 'omega':
                    omega_plot = omega[f][:, idx]
                    omega_label = r'$\omega$'
                else:
                    raise ValueError(f"Unknown omega_type: {omega_type}")

                style_legend.append(Line2D([0], [0], color='black', lw=2, ls=':', label=omega_label))
                        
                ax[i].plot(x_data_plot,
                           omega_plot,
                           color=band_colors[j],
                           lw=bands_lw,
                           ls=':',
                           alpha=alpha)
            # --------------------------------------
    # Add the Pw, P, w legend
    ax[2].legend(handles=style_legend)
                

                
       
        
        
        
        
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.01, hspace=0)
    
    
#     Add coloured text for each Band in the far right plot
#     ----------------------------------------------------
#     ----------------------------------------------------
    start = 0.9 # This is the starting vertical position
    c = 0
    marker_legend_handles = []
    
    # Loop over each band 
    for j, b_label in enumerate(bands_labels):
        
        label = rf'$\lambda$ = {b_label}'
        
#         if "robust" in b_label:
#             label = b_label.replace(" robust ", "\nrob. ")
#             c = c + 0.075 # This controls the vertical spacing 
#         else:
#             label = b_label
            
            
#         handle = Line2D(
#             [0], [0],
#             marker = ms[j],  
#             color='none',
#             markerfacecolor=band_colors[j],
#             markeredgecolor=band_colors[j],
#             markersize=20,
#             linestyle='None',
#             label=label)

#         marker_legend_handles.append(handle)

        if add_legend:
            ax[3].text(0.6, # This controls the horizontal position
                       start - c, 
                       label,
                       transform=ax[3].transAxes,
                       fontsize=30 * plot_sf,
                       color=band_colors[j])
        c = c + sigma
        
        
        

# 

    return fig, ax
    

# --------------------------------------------------------------------------------------






