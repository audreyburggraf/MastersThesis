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
from PlottingWithFunction import * 

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
# -----------------------------------------------------------------------------------------




def MakeAllBandGrid(bands,
                   fig_x_size = 20, 
                    fig_y_size = 6.5,
                    label = 'wavelength',
                    cmap = soft_colormap_v2,
                   vector_lw = constants.writeup_grid_lw, # 2.2
                   axis_num_fs = 16,
                axis_label_fs = 22, 
                   scalebar_fs = 16,
                   beam_lw = 1.5,
                   cbar_label_fs = 20,
                    # This is for the 100 AU and label text
                    x_pos = 0.06, 
                    y_pos = 0.85,
                    AU100_lw = 2.5,
                    AU100_fs = 20, 
                    text_fs = 20, 
                    wspace = 0.12,
                   for_poster = False,
                    for_writeup = False, 
                   ps = 1,
                   cbar_shrink = 1):
    
    
    #poster_text1_fs
    # Change sizes of things if for_poster or not
    # ---------------------------
    if for_poster:
        text_fs = poster_text1_fs * ps
        AU100_fs = constants.poster_text2_fs * ps
        spine_width = poster_spine_width 
        cbar_label_fs = constants.poster_text2_fs * ps
    # ---------------------------
    else: 
        text_fs = text_fs
        spine_width = 1
    # ---------------------------
    
    

    #bands = [B4, B5, B6, B7]
    B4, B5, B6, B7 = bands
    bands_idx = [4, 5, 6, 7]
    
    # These things should be the same for all bands
    StokesI_wcs = B4['StokesI_wcs']
    
    
    
    # Make the figure
    fig = plt.figure(figsize=(fig_x_size, fig_y_size))
    axes = []

    for i, b in enumerate(bands):
        
        
        
        ax = fig.add_subplot(1, 4, i + 1, projection=b['StokesI_wcs'])
        axes.append(ax)
        
        
        ax.coords.frame.set_linewidth(spine_width)

        
        # Plot POLI
        im = ax.imshow(b['POLI_debiased_mJy'], cmap = cmap)
        
        # Add a cbar
        
        # Make colorbar axis above plot
        cax = inset_axes(
            ax,
            width="100%",
            height="6%",
            loc="upper center",
            bbox_to_anchor=(0, 0.075, 1, 1),
            bbox_transform=ax.transAxes,
            borderpad=0
        )

        # Add colorbar
        cbar1 = plt.colorbar(im, cax=cax, orientation="horizontal")

        # Put ticks on top
        cbar1.ax.xaxis.set_ticks_position("top")

        # Format colorbar
        cbar1.ax.tick_params(
            labelsize=axis_num_fs,
            which="major",
            length=4,
            direction="in",
            width=spine_width
        )

        cbar1.outline.set_linewidth(spine_width)
        
        #cbar1 = plt.colorbar(im, cax = cax, orientation='horizontal', shrink=cbar_shrink, location = 'top')
        # cbar1.set_label('POLI (mJy beam$^{-1}$)', fontsize=cbar_fs)
#         cbar1.ax.tick_params(labelsize=axis_num_fs, which='major', length=4, direction="in", width=spine_width,)
#         cbar1.outline.set_linewidth(spine_width)
        #cbar1.set_label("Polarized Intensity [mJy beam$^{-1}$]",fontsize=cbar_label_fs)
        
        # Add the vectors
        for j, row in enumerate(b['VectorsActual']):
            if j == 0:
                ax.plot([row[0], row[1]], [row[2], row[3]], color='black', lw = vector_lw, label = 'Observations')
            else:
                ax.plot([row[0], row[1]], [row[2], row[3]], color='black', lw = vector_lw)
            
            
        
        
        xmin = b['xmin'] 
        xmax = b['xmax'] 
        ymin = b['ymin']
        ymax = b['ymax'] 
        
        print(rf'({xmin}, {xmax}), ({ymin}, {ymax})')
        
#         print(f'Band = {i}, and ymax = {ymax}')
    

        # Beam location
        beam_x_pos = xmin - 0.1 * (xmin - xmax)
        beam_y_pos = ymin - 0.1 * (ymin - ymax)

        
        # Add y-axis label to the left
        if i == 0:
            ax.set_ylabel("Declination", fontsize=axis_label_fs)
        
        # Add an x-axis label to all plots
        ax.set_xlabel(" ")
     
        # Set the x and y lims
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        
        # Add the line to each plot
        # ---------------------------------------------------------
        # Scalebar line
        line_x_pos = xmax - 0.04 * (xmax - xmin)
        line_y_pos = ymin - 0.14 * (ymin - ymax)  -10

        x0 = line_x_pos - b['reference_length_pix']
        x1 = line_x_pos

        ax.plot(
            [x0, x1],
            [line_y_pos, line_y_pos],
            color='black',
            linewidth=AU100_lw
        )

        # ---------------------------------------------------------
        # Convert midpoint from data coords -> axes coords
        x_mid_data = (x0 + x1) / 2
        y_mid_data = line_y_pos

        display_coords = ax.transData.transform((x_mid_data, y_mid_data))
        axes_coords = ax.transAxes.inverted().transform(display_coords)

        x_axes, y_axes = axes_coords

        # Text directly above middle of line
        ax.text(
            x_axes,
            y_axes + 0.015,
            f"{b['reference_length_AU']} AU",
            transform=ax.transAxes,
            fontsize=AU100_fs,
            ha='center',
            va='bottom'
        )
        
        add_band_label(ax, bands_idx[i], label, constants, x_pos = x_pos, y_pos = y_pos, 
                       label_fs = text_fs, va = 'bottom')
        # ---------------------------------------------------------
    
    
        # Add the beam
        beam = Ellipse((beam_x_pos, beam_y_pos), width=b['BMAJ_pix'], height=b['BMIN_pix'],  
                   angle=b['BPA_deg_cartesian'], edgecolor='black', facecolor='none', alpha=1, lw=beam_lw)
        
        ax.add_patch(beam)
        
        # Axis ticks and labels 
        # ---------------------------------------------------------
        # Only show y-axis numbers on far left plot
        if i != 0:
            ax.coords[1].set_ticklabel_visible(False)
            ax.coords[1].set_axislabel('')
        # ---------------------------------------------------------
        # Tick locations and style
        ax.coords[0].set_ticks(number=3, direction='in', width=spine_width)
        ax.coords[1].set_ticks(number=5, direction='in', width=spine_width)

        # Tick label sizes
        ax.coords[0].set_ticklabel(size=axis_num_fs)
        ax.coords[1].set_ticklabel(size=axis_num_fs)
#         # --------------------------------------------------------


        
    
    #fig.subplots_adjust(wspace=wspace)
    
    fig.supxlabel(
    "Right Ascension",
    fontsize=axis_label_fs,
    y=0.04
)
    if for_poster:
        fig.text(
            0.5,
            0.03,
            "Right Ascension",
            ha='center',
            fontsize=axis_label_fs
    )
    
#     fig.text(
#         0.025,
#         0.08,
#         r'Colour scale: ' '\n' 'Polarized Intensity' '\n' r'(mJy beam$^{-1}$)',
#         ha='left',
#         fontsize=15
#     )

    fig.text(
        0.5,          # centered
        0.995,         # adjust as needed
        r'Polarized Intensity [mJy beam$^{-1}$]',
        ha='center',
        va='top',
        fontsize=cbar_label_fs
    )
    
    
    return fig, axes