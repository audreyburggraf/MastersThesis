import numpy as np
import pandas as pd
import math 


from FITS_Image_Functions import * 
from C29_functions import * 


# Import necessary packages
import matplotlib.pyplot as plt

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
# -----------------------------------------------------------------------------------------



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
        # These values and equation are from Carniani, 2015 (dN/dS)
        phi = 1800 # deg^-2
        S0 = 1.7 # mJy
        alpha = -2.08 # unitless
        
        dN_dS = phi * (S/S0)**(alpha) * np.exp(-S/S0)
        
    elif band == 'Band 7':
        # Values from Statch, 2018 (https://arxiv.org/pdf/1805.05362), page 6, below equation 2
        N0    = 2.2 # mJy ^-1 deg^-2
        S0    = 3.9  # mJy
        alpha = -1.7  # unitless
        
        dN_dS = 1.15e3 * (N0 / S0) * (S / S0)**alpha * np.exp(-S/S0)
            
    else:
        return print('Currently function only supports Band 6 and Band 7')
    
    

    return dN_dS
# ---------------------------------------------------



# ---------------------------------------------------
# This function is used to find N from dN/dS using the trapezoidal rule

# I got ChatGPT to help me with log spacing i am eating lunch very busy need this done will clean it up later
# ---------------------------------------------------
def find_N_from_dN_dS(S_min, band, spacing):
    
    # Band 6 uses a linear spaced S grid for now 
    if spacing == 'linspace':  
        S_grid = np.logspace(-2, 2, 1000)  # 0.01 to 100 mJy
            
        dN_dS_grid = np.array([galaxy_probability(s, band) for s in S_grid])

        # Select fluxes above threshold
        mask = S_grid >= S_min
        
        N_per_area = np.trapz(dN_dS_grid[mask], S_grid[mask])

    elif spacing == 'logspace':
        
        # Use log-spaced bins and integrate in log-space
        S_grid = np.logspace(-2, 3, 2000)   # 0.01 to 1000 mJy
        logS_grid = np.log10(S_grid)
        
        # Compute dN/dS
        dN_dS_grid = np.array([galaxy_probability(s, band) for s in S_grid])

        
        # Convert to dN/dlog10S
        dN_dlogS_grid = S_grid * np.log(10) * dN_dS_grid
        
        # Select fluxes above threshold
        mask = S_grid >= S_min
        
        # Integrate in log-space
        N_per_area = np.trapz(dN_dlogS_grid[mask], logS_grid[mask])

        
    else:
        
        return print("Currently function only supports 'linspace' or 'logspace'")


    return N_per_area
# ---------------------------------------------------




# -----------------------------------------------------------------------------------------
def run_galaxy_contamination(band, sn_array, spacing = 'linspace', dr_arcsec = 0.25, print_things = True):
    if band == "Band 6":
        pb_path = constants.band6_data_folder_path + "c2d_989_pbeam_233GHz.fits"
        
        if print_things:
            print(rf"The Schecter function form that is being used is: dN_dS = phi * (S/S0)**(alpha) * np.exp(-S/S0)")
                
            print(rf"The Schecter best-fit parameters being used are: phi = 1800 deg^-2, S0 = 1.7 mJy, and alpha =  -2.08")
        
            print(rf"These values are from Carniani et al. (2015")

    elif band == "Band 7":
        pb_path = constants.band7_data_folder_path + "IRS63_BAND7_pb.fits"
        
        if print_things:
            print(rf"The Schecter function form that is being used is: dN_dS = 1.15e3 * (N0 / S0) * (S / S0)**alpha * np.exp(-S/S0)")
            
            print(rf"The Schecter best-fit parameters being used are: N0 = 2.2 mJy^-1 deg^-2, S0 = 3.9 mJy, and alpha = -1.7")
         
            
            print(rf"These values are from Chen et al. (2022)")
        
    else:
        return print('Currently function only supports Band 5 and Band 6')
    
    print(" ")
    print(rf"The dr value is: {dr_arcsec} arcsec")
    print(" ")
    
    # Primary Beam
    # ------------------------------------------------------------------------------------------------
    header, _, pb_map, wcs = read_in_file(pb_path)

    ny, nx = pb_map.shape
    
    flat_index = np.nanargmax(pb_map)   # ignores NaNs, finds max
    y_cen, x_cen = np.unravel_index(flat_index, pb_map.shape)

    if print_things:
        print(rf'The centre position is: (y, x) = ({y_cen}, {x_cen}), and at this position pb = {pb_map[y_cen, x_cen]}')
        print(" ")
    # ------------------------------------------------------------------------------------------------
    

    # Size of annuli
    # ------------------------------------------------------------------------------------------------
    dr_deg = arcsec_to_degrees(dr_arcsec)

    dr_pix = arcsec_to_pixels(header, dr_arcsec)
    # ------------------------------------------------------------------------------------------------


    # Radius
    # ------------------------------------------------------------------------------------------------
    r_pix = np.arange(0, nx - x_cen, dr_pix)
    r_arcsec = pixels_to_arcsec(header, r_pix)


    df = pd.DataFrame({
        "r_pix": r_pix,
        "r_arcsec": r_arcsec,
        "r_deg": arcsec_to_degrees(r_arcsec)
    })
    # ------------------------------------------------------------------------------------------------

    
    # sigma, area
    # ------------------------------------------------------------------------------------------------
    x_pos = x_cen + r_pix.astype(int)

    pb_val_arr = pb_map[y_cen, x_pos].astype(np.float32).newbyteorder('=')  # native endian
    
    # to do: find sentiviity for BAND 7 
    sigma_r_arr = (0.035 / pb_val_arr).astype(np.float32).newbyteorder('=')

    df["pb_val"] = pb_val_arr
    df["sigma_r"] = sigma_r_arr

    df["area_annulus"] = 2 * np.pi * df["r_deg"] * dr_deg
    # ------------------------------------------------------------------------------------------------

    
    
    # Dictionary to hold each df
    dfs_sn = {}
    dfs_list = [] 

    for sn in sn_array:
        # Compute flux threshold per annulus
        S = sn * df["sigma_r"].values

        # Compute dN/dS
        dN_dS = np.array([galaxy_probability(s, band) for s in S])

        N_per_area_list = []

        # Loop over each annulus
        for S_min in S:
            # Use your function to integrate dN/dS above threshold
            N_per_area = find_N_from_dN_dS(S_min, band, spacing)
            N_per_area_list.append(N_per_area)

        N = np.array(N_per_area_list)

        # Multiply by annulus area to get expected number per annulus
        N_expected = N * df["area_annulus"].values

        prob = 100*(1 - np.exp(-N_expected))

        N_cumulative = np.cumsum(N_expected)
        prob_cumulative = 100 * (1 - np.exp(-N_cumulative))

        df_sn = pd.DataFrame({
            "r_arcsec": df["r_arcsec"].values,
            "S_mJy": S,
            "dN_dS": dN_dS,
            "N": N,
            "N_expected": N_expected,
            "prob": prob,
            "prob_cumulative": prob_cumulative
        })

#         # Save in dictionary
#         dfs_sn[f"df_{sn}sigma"] = df_sn
        
        
        dfs_list.append(df_sn)





    return y_cen, x_cen, df, dfs_list
# -----------------------------------------------------------------------------------------






# Function to find probability
# -----------------------------------------------------------------------------------------
def GC_probability(x_cen, y_cen, df, dfs, StokesI_mJy, sn_array, band, distance = "x", print_things = True):
    
    # Note: CARTA and the fits maps are different
    # If on CARTA the source is at Image: (50, 25)
    # I need to index my StokesI_mJy map as [25, 50]
    if band == "Band 6":
        S_source_mJy = 1 * 1000 # From CARTA
        x_source = 1 # From CARTA
        y_source = 1 # From CARTA
        
        x_dist = 1 # From CARTA
        y_dist = 1 # From CARTA
        hypot_dist = 1 # From CARTA

        err_mJy = constants.StokesI_err_mJy_band6
        
    elif band == "Band 7":
        S_source_mJy = 4.008985706605e-4 * 1000 # From CARTA
        x_source = 773 # From CARTA
        y_source = 472 # From CARTA
        
        x_dist = 4.752000 # From CARTA
        y_dist = 0.684000 # From CARTA
        hypot_dist = 4.800975 # From CARTA
        
        err_mJy = constants.StokesI_err_mJy_band7 / (0.827)
        
        print("S/N = ", S_source_mJy/err_mJy)
       
        
    else:
        return print('Currently function only supports Band 5 and Band 6')
    
    
    # Check that the inputted S_mJy at (x_source, y_source) matches my StokesI_mJy array
    if math.isclose(StokesI_mJy[y_source, x_source], S_source_mJy, rel_tol = 0.001):
        if print_things: print("     Checkpoint: Stokes I inputted and mapped values match")
    else: return("The inputted values of Stokes I and the indexed ones do not match within the tolerence")          
    
         
    # Find the S/N ratio of the source
    S_N_source = S_source_mJy / err_mJy
    
    # Now check that the hypotenuse value from CARTA is correct based on x_dist and y_dist 
    calculated_hypot = np.sqrt(x_dist**2 + y_dist**2)
                  
    if math.isclose(np.sqrt(x_dist**2 + y_dist**2), hypot_dist, rel_tol = 0.001):
        if print_things: print("     Checkpoint: The hypotenuse values match up")
    else: return("The hypotenuse values do not match within the tolerence")
                  
    # Choose what distance value we are looking at based on the input  
    if distance == "x":
        d_source = x_dist
    elif distance == "y": 
        d_source = y_dist                  
    elif distance == "h" :
        d_source = hypot_dist 
    else:
        return print("Currently function only supports Distance of 'x', 'y', or 'h'")    
                  
    # Print the distance from the centre to the object 
    if print_things: print(" ")     
    if print_things: print(rf"The distance to the source in {distance} is  : {d_source} arcsec")
                  
    # Find the closest r value to the distance from dr, and its corresponding index
    closest_r = df["r_arcsec"].iloc[(df["r_arcsec"] - d_source).abs().argmin()]
    closest_idx = (df["r_arcsec"] - x_dist).abs().argmin()
                  
    if print_things: print(rf"The distance to the soruce in df is : {closest_r } arcsec")
    if print_things: print(" ")            
                  
    # Find the probability of an object of the S/N valye at its corresponding distance
    for i in range(len(sn_array)):

        print(rf"If the source S/N > {sn_array[i]}:")
        print(rf"     The cumulative prob is: {dfs[i]['prob_cumulative'][closest_idx]:.2f}%")
        print(" ")              
                  

# -----------------------------------------------------------------------------------------





# I will also make some functions to make plotting easier




# -----------------------------------------------------------------------------------------
# This function makes a plot to replicate Figure 11 and 12 from Sadavoy 2019 
# (https://iopscience.iop.org/article/10.3847/1538-4365/ab4257)
# -----------------------------------------------------------------------------------------
galaxy_ls = ['-', '--', '-.']


def recreate_figs_11_12_sadavoy2019(dfs, sn_array, 
                                    fig_size_x = 20, fig_size_y = 7):
    
    # sn_array will usually be [5, 10, 20] but I will leave it open ended 
    
        # make the figures
    fig, ax = plt.subplots(1, 2, figsize=(fig_size_x, fig_size_y))
    
    for i in range(len(sn_array)):

        # Plot data on left plot
        ax[0].plot(dfs[i]["r_arcsec"], dfs[i]["prob"], color = 'black', ls = galaxy_ls[i], label=f'S/N > {sn_array[i]}')

        # Cumulative probability
        ax[1].plot(dfs[i]["r_arcsec"], dfs[i]["prob_cumulative"], color = 'black', ls = galaxy_ls[i], label=f'S/N > {sn_array[i]}')

    
    
    # Set y-axis labels 
    ax[0].set_ylabel('Probability of Galaxy in Annulus (%)', fontsize = axis_label_fs)
    ax[1].set_ylabel('Cumulative Probability of Galaxy(%)',  fontsize = axis_label_fs)
    
    
    # Add an x-label, proper ticks on axes and add a legend
    for i in range(2):
    
        # Add axis labels
        ax[i].set_xlabel(r'Radius (arcsec)', fontsize=axis_label_fs)


        # Adjust ticks
        ax[i].minorticks_on()
        ax[i].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
        ax[i].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)
    
        ax[i].legend(fontsize = legend_text_fs)
    
    
    return fig, ax
# -----------------------------------------------------------------------------------------





# -----------------------------------------------------------------------------------------
# This plot is so we can see how pb, sigma and area behave with change in radius
# -----------------------------------------------------------------------------------------
def plot_galaxy_contamination_models(df,
                                     fig_size_x = 18, fig_size_y = 10):
    
    # Define the grid and size 
    fig, ax = plt.subplots(2, 2, figsize=(fig_size_x, fig_size_y))
    
    
    # Top left 
    ax[0, 0].plot(df["r_arcsec"], df["pb_val"], color = 'tab:orange')
    ax[0, 0].set_ylabel('pb val', fontsize=axis_label_fs)


    # Top right 
    ax[0, 1].plot(df["r_arcsec"], df["sigma_r"], color = 'tab:blue')
    ax[0, 1].set_ylabel('$\sigma(r)$ = 0.03/pb val', fontsize=axis_label_fs)
    
    # Bottom right 
    ax[1, 0].plot(df["r_arcsec"], df["area_annulus"], color = 'tab:red')
    ax[1, 0].set_ylabel('area [deg$^2$]', fontsize=axis_label_fs)
    
    
    for i in range(2):
        for j in range(2):
            ax[i, j].minorticks_on()
            ax[i, j].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
            ax[i, j].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)
            ax[i, j].tick_params(which='minor', length=4, direction="in")
            
            ax[i, j].set_xlabel('Radius (arcsec)', fontsize=axis_label_fs)  
    
    
    
    return fig, ax
# -----------------------------------------------------------------------------------------





# -----------------------------------------------------------------------------------------
# Here we will plot S, dN/dS, N, etc
# -----------------------------------------------------------------------------------------
def plot_probability_models(dfs, sn_array,
                            fig_size_x = 18, fig_size_y = 10):
    
    # Define the grid and size 
    fig, ax = plt.subplots(2, 2, figsize=(fig_size_x, fig_size_y))
    
    for i in range(len(sn_array)):
        # Top left 
        ax[0, 0].plot(dfs[i]["r_arcsec"],  dfs[i]["S_mJy"], color = 'black', ls = galaxy_ls[i], label=f'S/N > {sn_array[i]}')
        ax[0, 0].set_ylabel('S_mJy', fontsize=axis_label_fs)
        ax[0, 0].legend(fontsize = legend_text_fs)

        # TopRight
        ax[0, 1].plot(dfs[i]["r_arcsec"], dfs[i]["dN_dS"], color = 'black', ls = galaxy_ls[i], label=f'S/N > {sn_array[i]}')
        ax[0, 1].set_ylabel('dN/dS', fontsize=axis_label_fs)
        ax[0, 1].legend(fontsize = legend_text_fs)


        # Bottom Left
        ax[1, 0].plot(dfs[i]["r_arcsec"], dfs[i]["N"], color = 'black', ls = galaxy_ls[i], label=f'S/N > {sn_array[i]}')
        ax[1, 0].set_ylabel('N (from dN/dS using trap. rule)', fontsize=axis_label_fs)
        ax[1, 0].legend(fontsize = legend_text_fs)

        # Bottom Right
        ax[1, 1].plot(dfs[i]["r_arcsec"], dfs[i]["N_expected"], color = 'black', ls = galaxy_ls[i], label=f'S/N > {sn_array[i]}')
        ax[1, 1].set_ylabel('N * area (#, not %)', fontsize=axis_label_fs)
        ax[1, 1].legend(fontsize = legend_text_fs)

    
    
    for i in range(2):
        for j in range(2):
            ax[i, j].minorticks_on()
            ax[i, j].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
            ax[i, j].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)
            ax[i, j].tick_params(which='minor', length=4, direction="in")
            
            if i == 1:
                ax[i, j].set_xlabel('Radius (arcsec)', fontsize=axis_label_fs)  
    
    
    fig.subplots_adjust(
    left=0.08,   # space on left
    right=0.95,  # space on right
    bottom=0.10, # space on bottom
    top=0.92,    # space on top
    wspace=0.25, # horizontal space between subplots
    hspace=0.20  # vertical space between subplots
)
    
    return fig, ax
# -----------------------------------------------------------------------------------------










