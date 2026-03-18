# Import necessary packages
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns
import math
from matplotlib.ticker import ScalarFormatter
import sys

# Import Functions
from DustModelFunctions import *


# Add the directory where constants.py is located to sys.path
sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/")  # Replace with the actual path

# Now you can import constants.py
import constants


# Set the colours for each band 
# band1_color  = "#f63fe3"  # pink from Kataoka et al. 2015
# band3_color  = "#2321f6"  # blue from Kataoka et al. 2015
# band4_color  = "#d4aa00"  # warm gold/mustard yellow that ChatGPT reccomended
# band6_color  =  "#00bfc4" # bright cyan that ChatGPT reccomended
# band7_color  = "#086608"  # green from Kataoka et al. 2015
# band10_color = "#fb0d0d" # red from Kataoka et al. 2015

# Set the colours for each band 
# alma_band_colors = {
#     1: "#f63fe3",  # pink from Kataoka et al. 2015
#     3: "#2321f6",  # blue from Kataoka et al. 2015
   
#     #--------
#     # Band 4
#     #--------
#     4: "red", 
#     "Band 4": "red",
#     "Band 4 nterms2": 'red',
#     #--------
#     # Band 5
#     #--------
#     5: "orange",  
#     "Band 5": "orange",
#     #--------
#     # Band 6
#     #--------
#     6: "#00bfc4", # Bright cyan
#     "Band 6": "#00bfc4",
#     #--------
#     # Band 7
#     #--------
#     7: "darkviolet", 
#     "Band 7":  "darkviolet",
#     "Band 7 nterms2": "darkviolet",
#     10:"#fb0d0d", # red from Kataoka et al. 2015
# }

alma_band_colors = constants.alma_band_colors

# alma_band_ms = {
#     1: '.', 
#     3: 'P', 
    
#     4: 'o',  
#     'Band 4': 'o',  
#     'Band 4 nterms2': 'o',  
#     'Band 4 nterms2 smooth': 'o', 
#     'Band 4 nterms2 smooth B6': 'o', 
#     'Band 4 nterms2 smooth B6 B7': 'o', 
    
#     5: '^',  
#     50: '^', 
#     52: '^', 
#     'Band 5': '^',  
#     'Band 5 v0': '^', 
    
#     # Band 5 robust values
#     'Band 5 robust -2': 'd',
#     'Band 5 robust -1': 'P',
    
#     6: 's', 
#     'Band 6': 's', 
#     'Band 6 smooth': 's', 
#     'Band 6 smooth B7': 's', 
    
#     7: '*',  
#     'Band 7': '*', 
#     'Band 7 nterms2': '*', 
#     'Band 7 nterms2 smooth': '*', 
#     'Band 7 nterms2 smooth B6': '*', 
    
#     10:'p', 
# }


kataoka_fig2_colors = fig2_colors = ['red', 'orange', 'blue', 'purple']







def plot_a_vs_POLF(a_max_dist_micron, P_times_omega,
                   sf, bands, y_plot, legend, 
                   text_fs = constants.text_fs , 
                   axis_label_fs = constants.axis_label_fs, 
                   axis_num_fs = constants.axis_num_fs, 
                   legend_text_fs = constants.legend_text_fs,
                   fig_size_x = 10,
                   fig_size_y = 6,
                   
   ):
    
    # Get bands colors, and marker style
    band_colors = [alma_band_colors[b] for b in bands]
    ms = [alma_band_ms[b] for b in bands]
    
    
    
    # Make list for labels
    # -----------------------------------
    labels = []

    for band in bands:

        if band == 'Band 4' or band == 'Band 4 nterms2':
            label = 'Band 4'
        elif band == 'Band 7' or band == 'Band 7 nterms2':
            label = 'Band 7'
        else:
            label = band

        labels.append(label)

    labels = np.array(labels)
    # -----------------------------------

    # mm

    # Create a figure with the WCS projection
    fig, ax = plt.subplots(figsize=(fig_size_x, fig_size_y))


    for i in range(len(bands)):
        if legend == 'w sf':
            label = f'{labels[i]} $\cdot$ sf'
        elif legend == 'wo sf':
            label = f'{labels[i]}'
        else:
            label = None  # in case legend is neither 'with sf' nor 'wo sf'
        
        
        ax.plot(a_max_dist_micron,
                P_times_omega[:, i] * sf,
                color = band_colors[i],
                ls = '-',
                label = label, 
                lw = 3)


    # Add axis labels
    if y_plot == 'a':
        ax.set_xlabel('Maximum grain size [$\mu$m]', fontsize=axis_label_fs)
    elif y_plot == 'af':
        ax.set_xlabel(rf"$a_{{\rm max}} \cdot f$", fontsize=axis_label_fs)
        
    ax.set_ylabel('POLF', fontsize=axis_label_fs)


    # Adjust ticks
    ax.minorticks_on()

    ax.tick_params(axis="x", which="minor", direction="in", left=True, right=True, length=4)
    ax.tick_params(axis="y", which="minor", direction="in", left=True, right=True, length=4)

    ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True,   length=7, labelsize=axis_num_fs)
    ax.tick_params(axis="y", which="major", direction="in", left=True,   right=True, length=7, labelsize=axis_num_fs)



    # Make axes log if needed
    ax.set_xscale('log')
    # ax.set_yscale('log')

    # Remove scientific notation from x-axis
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
    ax.ticklabel_format(style='plain', axis='x')

    # ---------------------------------------------------------------------------------
#     ax.legend(fontsize = legend_text_fs, title = 'Model', title_fontsize = legend_text_fs)

    
    return fig, ax