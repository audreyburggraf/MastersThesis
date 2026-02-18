import pandas as pd
import numpy as np
from scipy.optimize import *
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import seaborn as sns

# Import the constants
# -----------------------------------------------------------------------------------------
import sys

sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/") 

import constants

title_fs = constants.title_fs
axis_label_fs = constants.axis_label_fs
axis_num_fs = constants.axis_num_fs
legend_title_fs = constants.legend_title_fs
legend_text_fs = constants.legend_text_fs
cbar_fs = constants.cbar_fs
text_fs = constants.text_fs
# -----------------------------------------------------------------------------------------



# Data
# -----------------------------------------------------------------------------------------------------
# Make a df to store the data from Chen Table 2 
df_cumulative_Chen = pd.DataFrame({
    "S"              : np.array([0.19,  0.32,  0.52, 0.86, 1.43, 2.36, 3.91, 6.47]), # [mJy]
    "N"              : np.array([92,    84,    69,   53,   39,   25,   13,   4]),    # [unitless]
    "N(>S)"          : np.array([21900, 12950, 7390, 4120, 2660, 1540, 770,  230]),  # [deg^-2]
    "N(>S) err lower": np.array([3250,  1640,  930,  560,  420,  300,  210,  110]),  # [deg^-2] 
    "N(>S) err upper": np.array([3250,  1640,  930,  560,  420,  370,  270,  180])   # [deg^-2] 
})

cumulative_data = [df_cumulative_Chen["S"], 
                   df_cumulative_Chen["N(>S)"],
                   df_cumulative_Chen["N(>S) err lower"], 
                   df_cumulative_Chen["N(>S) err upper"]]
# df_cumulative_Chen[""] = ["|" for S in df_cumulative_Chen["S"]]
# -----------------------------------------------------------------------------------------------------
# Make a df to store the differential data from Chen Table 2
df_differential_Chen = pd.DataFrame({
    "S"              : np.array([0.25,  0.42,  0.69,  1.15, 1.90, 3.14, 5.19]), # [mJy]
    "N"              : np.array([8,     15,    16,    14,   14,   12,   9]),    # [unitless]
    "dN/dS"          : np.array([73860, 28600, 10230, 2780, 1280, 530,  220]),  # [mJy^-1 deg^-2]
    "dN/dS err lower": np.array([29940, 8210,  2670,  740,  340,  150,  70]),   # [mJy^-1 deg^-2] 
    "dN/dS err upper": np.array([39660, 10180, 3350,  960,  440,  200,  100])   # [mJy^-1 deg^-2] 
})

differential_data = [df_differential_Chen["S"], 
                     df_differential_Chen["dN/dS"],
                     df_differential_Chen["dN/dS err lower"], 
                     df_differential_Chen["dN/dS err upper"]]

# df_diff_Chen[""] = ["|" for S in df_diff_Chen["S"]]
# -----------------------------------------------------------------------------------------------------
# Make a df to store the data from Table 3
BestFitValues = {
    # Schecter: Cumulative
    "schechter":{
        "cumulative":{
            "alpha"          : -0.6,  # [unitless]
            "S0"             : 2.2,   # [mJy]
            "N0"             : 9.6,   # [deg^-2]   
            "alpha_err_upper": 0.1,   # [unitless]
            "alpha_err_lower": -0.1,  # [unitless]
            "S0_err_upper"   : 0.1,   # [mJy]
            "S0_err_lower"   : -0.1,  # [mJy]
            "N0_err_upper"   : 0.9,   # [deg^-2]  
            "N0_err_lower"   : -0.9,  # [deg^-2]   
            },
      # Schecter: Differential
        "differential":{
            "alpha"          : -1.7,  # [unitless]
            "S0"             : 3.9,   # [mJy]
            "N0"             : 2.2,   # [mJy^-1 deg^-2]   
            "alpha_err_upper": 0.2,   # [unitless]
            "alpha_err_lower": -0.2,  # [unitless]
            "S0_err_upper"   : 0.7,   # [mJy]
            "S0_err_lower"   : -0.7,  # [mJy]
            "N0_err_upper"   : 1.4,   # [mJy^-1 deg^-2]  
            "N0_err_lower"   : -0.9,  # [mJy^-1 deg^-2]   
            },
        "joint_fit":{
            "alpha_alpha1"          : -1.7,  # [unitless]
            "S0"                    : 3.9,   # [mJy]
            "N0"                    : 2.1,   # [mJy^-1 deg^-2]   
            "alpha_alpha1_err_upper": 0.1,   # [unitless]
            "alpha_alpha1_err_lower": -0.1,  # [unitless]
            "S0_err_upper"          : 0.5,   # [mJy]
            "S0_err_lower"          : -0.5,  # [mJy]
            "N0_err_upper"          : 0.7,   # [mJy^-1 deg^-2] 
            "N0_err_lower"          : -0.5,  # [mJy^-1 deg^-2]   
            }
      },
    # Double power law: Cumulative
    "dpl":{
        "cumulative":{ 
            "alpha1"           : 1.0,   # [unitless]
            "alpha2"           : 4.7,   # [mJy]
            "S0"               : 4.5,   # [mJy]
            "N0"               : 3.6,   # [deg^-2]   
            "alpha1_err_upper" : 0.1,   # [unitless]
            "alpha1_err_lower" : -0.1,  # [unitless]
            "alpha2_err_upper" : 0.2,   # [unitless]
            "alpha2_err_lower" : -0.2,  # [unitless]
            "S0_err_upper"     : 0.4,   # [mJy]
            "S0_err_lower"     : -0.4,  # [mJy]
            "N0_err_upper"     : 0.7,   # [deg^-2]  
            "N0_err_lower"     : -0.6,  # [deg^-2]   
            },
    # Double power law: Differential
        "differential":{ 
            "alpha1"           : 2.0,   # [unitless]
            "alpha2"           : 4.7,   # [mJy]
            "S0"               : 5.6,   # [mJy]
            "N0"               : 0.8,   # [mJy^-1 deg^-2]   
            "alpha1_err_upper" : 0.2,   # [unitless]
            "alpha1_err_lower" : -0.3,  # [unitless]
            "alpha2_err_upper" : 0.9,   # [unitless]
            "alpha2_err_lower" : -0.5,  # [unitless]
            "S0_err_upper"     : 2.2,   # [mJy]
            "S0_err_lower"     : -1.5,  # [mJy]
            "N0_err_upper"     : 1.6,   # [mJy^-1 deg^-2]  
            "N0_err_lower"     : -0.8,  # [mJy^-1 deg^-2]   
            },
        # double power law: joint fit
        "joint_fit":{
            "alpha1"                : 1.9,   # [unitless]
            "alpha2"                : 4.2,   # [unitless]
            "S0"                    : 4.3,   # [mJy]
            "N0"                    : 1.3,   # [mJy^-1 deg^-2]   
            "alpha_alpha1_err_upper": 0.1,   # [unitless]
            "alpha_alpha1_err_lower": -0.1,  # [unitless]
            "alpha_alpha2_err_upper": 0.3,   # [unitless]
            "alpha_alpha2_err_lower": -0.3,  # [unitless]
            "S0_err_upper"          : 0.7,   # [mJy]
            "S0_err_lower"          : -0.6,  # [mJy]
            "N0_err_upper"          : 0.6,   # [mJy^-1 deg^-2]  
            "N0_err_lower"          : -0.5,  # [mJy^-1 deg^-2]   
            }
    }
}
# You can index each value like this: B7["differential"]["schechter"]["S0"]
# -----------------------------------------------------------------------------------------------------


# Functions
# -----------------------------------------------------------------------------------------------------
# This is the schecter function from the code the Chen paper gave us 
def schechter_chencode(s, a=0, S0=1, N0=1):
    return N0*3600/S0*(s/S0)**(a+1)*np.exp(-s/S0)
# -----------------------------------------------------------------------------------------------------



# This is code directly from the paper that they provided us
# -----------------------------------------------------------------------------------------------------
from astropy.table import Table


def flux_scaling(wavelength, target_wavelength, flux=1, SED_template='data/uds_composite_sed_Lir_ld.dat', z=2.0,
                 debug=False):
    """interpolate the flux at a specified wavelength
    
    Params:
        wavelength : in um
        flux : any units, self consistent
    """
    SED = np.loadtxt(SED_template)
    f_SED = interpolate.interp1d(SED[:,0]*(1.0+z), SED[:,1]/(1+z)**4)
    scale_factor = flux / f_SED(wavelength)
    if debug:
        print("flux at {}: {}".format(wavelength, f_SED(wavelength)))
        print("flux at {}: {}".format(target_wavelength, f_SED(target_wavelength)))
        print("scale_factor: {}".format(scale_factor))
        print("flux difference: {}".format(f_SED(wavelength)/f_SED(target_wavelength)))
    return scale_factor * f_SED(target_wavelength)



S_alma_scale = {}
S_alma = {'B3': 3027.9403678720932, 'B4': 2085.4643523587724, 'B5': 1595.2429758332012, 'B6': 1259.8894674554122, 
          'B7': 886.746475901551, 'B8': 673.4616464220171, 'B9': 444.86026069217337, 'B10': 341.799906987015}
S_ref = {'B3': 3000., 'B4': 2000., 'B5': 1500., 'B6': 1200., 'B7': 870., 'B8': 650., 'B9': 440., 'B10': 340.}
for band in ['B3','B4','B5','B6','B7','B8','B9','B10']:
    S_alma_scale[band] = flux_scaling(S_alma[band], S_ref[band], z=2.0)
# print(S_alma_scale)



S_almacal_B7toS870 = 1.07
# -----------------------------------------------------------------------------------------------------





# Comparison plot
# -----------------------------------------------------------------------------------------------------

def plot_Chen2022():
    # Set the size and layout of the figures
    fig, ax = plt.subplots(1, 2, figsize=(20, 7))
    
    
    # Set colors 
    # ---------------------------------------------
    my_line_color = 'tab:red'
    chen_line_color = 'tab:cyan'

    my_data_point_color = 'green'
    chen_data_point_color = 'darkblue'
    # ---------------------------------------------
    
    
    # Set y-axis labels 
    # ----------------------------------------------------------------------------------------------
    ax[0].set_ylabel('dN/dS x $S^{2.5}$ [mJy$^{-1.5}$ deg$^{-2}$]', fontsize = axis_label_fs)
    ax[1].set_ylabel('N(>S) [deg$^{-2}$]',  fontsize = axis_label_fs)
    # ----------------------------------------------------------------------------------------------
    
    
   
    
    
    # Add the green markers using the data from Table 2 to the left plot
    # ----------------------------------------------------------------------------------------------
    y_error_lower_differential = df_differential_Chen["dN/dS err lower"] * df_differential_Chen["S"]**2.5
    y_error_upper_differential = df_differential_Chen["dN/dS err upper"] * df_differential_Chen["S"]**2.5
    
    y_error_differential = np.array([y_error_lower_differential, y_error_upper_differential])
    
    
    ax[0].errorbar(
        df_differential_Chen["S"],
        df_differential_Chen["dN/dS"] * df_differential_Chen["S"]**2.5,
        yerr = y_error_differential,
        fmt        = 'D',                  # diamond marker
        color      = my_data_point_color,  # Set the colour
        markersize = 8,
        capsize    = 5,                    # little horizontal bar at the ends
        elinewidth = 1.5,
        label      = 'Table 2 Data Points')
    # ----------------------------------------------------------------------------------------------
    
    
    
    # Add the green markers using the data from Table 2 to the right plot
    # ----------------------------------------------------------------------------------------------
    y_error_lower_cumulative = df_cumulative_Chen["N(>S) err lower"]
    y_error_upper_cumulative = df_cumulative_Chen["N(>S) err upper"]
    
    y_error_cumulative = np.array([y_error_lower_cumulative, y_error_upper_cumulative])
    
    
    ax[1].errorbar(
        df_cumulative_Chen["S"],
        df_cumulative_Chen["N(>S)"],
        yerr = y_error_cumulative,
        fmt        = 'D',                  # diamond marker
        color      = my_data_point_color,  # Set the colour
        markersize = 8,
        capsize    = 5,                    # little horizontal bar at the ends
        elinewidth = 1.5,
        label      = 'Table 2 Data Points')
    # ----------------------------------------------------------------------------------------------
    
    
    
    
    
    # Add an x-label, proper ticks on axes and add a legend
    # ----------------------------------------------------------------------------------------------
    for i in range(2):

        # Add axis labels
        ax[i].set_xlabel(r'Flux density [mJy]', fontsize=axis_label_fs)

        # Adjust ticks
        ax[i].minorticks_on()
        ax[i].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
        ax[i].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)

        ax[i].tick_params(axis="x", which="minor", direction="in", bottom=True, top=True, length=4, labelsize=axis_num_fs)
        ax[i].tick_params(axis="y", which="minor", direction="in", left=True, right=True, length=4, labelsize=axis_num_fs)

#         ax[i].legend(fontsize = legend_text_fs)

        ax[i].set_xscale("log")
        ax[i].set_yscale("log")

        ax[i].xaxis.set_major_formatter(ScalarFormatter())
        ax[i].ticklabel_format(axis="x", style="plain")
    # ----------------------------------------------------------------------------------------------

    # Set the bounds
    # ------------------------------------------
    # Bounds of left plot from the paper code
    ax[0].set_xlim(0.06, 30)
    ax[0].set_ylim(5., 1e5)
    # ------------------------------------------
    # Bounds of right plot from the paper code
    ax[1].set_xlim(0.03, 50)
    ax[1].set_ylim(0.1, 10**5)
    # ------------------------------------------
    
    return fig, ax
# -----------------------------------------------------------------------------------------------------