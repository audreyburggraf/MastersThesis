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
from UnitConversion import * 

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
lambda_mm = constants.lambda_mm


writeup_grid_axis_label_fs = constants.writeup_grid_axis_label_fs 
writeup_grid_axis_num_fs = constants.writeup_grid_axis_num_fs 
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
            f"{b['reference_length_AU']} au",
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












from astropy.coordinates import SkyCoord, SkyOffsetFrame
# from matplotlib.ticker import FuncFormatter

# def ra_offset_formatter(x, pos, centre_pix):
#     return f'{(x - centre_pix) * 0.018:.1f}'




# ChatGPT wrote this for me
# ------------------
def format_writeup_grid(
    fig, cbar_label,
    stokes_i_y=0.97,
    ra_y=0.05,
    dec_x=0.05,
    wspace=0.085
):

    fig.text(
        0.5,
        stokes_i_y,
        cbar_label,
        ha='center',
        va='top',
        fontsize= constants.writeup_grid_axis_label_fs
    )

    fig.supxlabel(
        r"$\Delta$ RA [arcsec]",
        fontsize=constants.writeup_grid_axis_label_fs,
        y=ra_y
    )

    fig.supylabel(
        r"$\Delta$ Dec [arcsec]",
        fontsize=constants.writeup_grid_axis_label_fs,
        x=dec_x
    )

    fig.subplots_adjust(wspace=wspace)
# ------------------    

def MakeAllBandGridDelta(bands,
                         centre_strings,
                         plotting, 
                   fig_x_size = 20, 
                    fig_y_size = 6.5,
                    label = 'wavelength',
                    cmap = soft_colormap_v2,
                   vector_lw = 2.2,
                   axis_num_fs = 16,
                axis_label_fs = 22, 
                   scalebar_fs = 16,
                   beam_lw = 1.5,
                   cbar_label_fs = 20,
                    # This is for the 100 AU and label text
                    x_pos = 0.06, 
                    y_pos = 0.85,
                    AU100_lw = 2.5,
                    text_fs = 20, 
                    wspace = 0.12,
                   for_poster = False,
                    for_writeup = False, 
                   ps = 1,
                   cbar_shrink = 1,
                    show_vectors = True):
    
    
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
        AU100_fs = 20
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
        
        
        
        
        # Get the centre position of the disk
        centre_string = centre_strings[i]

        centre_string = centre_strings[i]

        ra_part = centre_string[1:11]
        dec_part = centre_string[11:]


        ra_string = (
        ra_part[:2] + 'h' +
        ra_part[2:4] + 'm' +
        ra_part[4:] + 's'
        )

        dec_string = (
            dec_part[:3] + 'd' +
            dec_part[3:5] + 'm' +
            dec_part[5:] + 's'
        )

#         print(rf'ra_string = {ra_string}')
#         print(rf'dec_string = {dec_string}')

        centre = SkyCoord(
            ra_string,
            dec_string,
            frame='icrs'
        )
        

        
        
        
        ax = fig.add_subplot(1, 4, i + 1, projection=b['StokesI_wcs'])
        axes.append(ax)
        
        offset_frame = SkyOffsetFrame(origin=centre)

        overlay = ax.get_coords_overlay(offset_frame)
        
        ax.coords.frame.set_linewidth(spine_width)
        
    



        # Plot based on parameter
        # ---------------------------------------------------------
        if plotting == 'Stokes I':
            cbar_label = r'Stokes I [mJy beam$^{-1}$]'
            
            im = ax.imshow(b['StokesI_mJy'], cmap = cmap)
            #im = ax.imshow(b['StokesI_stretched_mJy'], cmap = cmap)
        # ---------------------------------------------------------
        elif plotting == 'Stokes I stretched':
            cbar_label = r'Stokes I [mJy beam$^{-1}$]'
            
            # im = ax.imshow(b['StokesI_mJy'], cmap = cmap)
            im = ax.imshow(b['StokesI_stretched_mJy_grid'], cmap = cmap)
            
        # ---------------------------------------------------------
        elif plotting == 'POLI':
            cbar_label = r'Polarized Intensity [mJy beam$^{-1}$]'
            
            im = ax.imshow(b['POLI_debiased_mJy'], cmap = cmap)
        # ---------------------------------------------------------
        elif plotting == 'Gaussian Map':
            cbar_label = r"$\mathrm{Gaussian\ Uniform\ Weights}\ [W_{\mathrm{Uniform}}]$"
            
            im = ax.imshow(b['ratios'], cmap = soft_colormap_no_red)
        # ---------------------------------------------------------
        
        else:
            return("parameter 'plotting' needs to equal 'Stokes I' or 'POLI'")
        # --------------------------------------------------------
        
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
        
        if plotting == 'Stokes I stretched':
       
            
            cbar1.set_ticks(constants.normalized_cbar_ticks_grid)
            cbar1.set_ticklabels([f"{val:.1f}" for val in b['StokesI_unstretched_cbar_ticks_grid']])
            
#             if i == 0:
#                 cbar1.set_ticklabels(
#                     [f"{val:.2f}" for val in b['StokesI_unstretched_cbar_ticks_grid']]
#                 )
#             else:
#                 cbar1.set_ticklabels(
#                     [f"{val:.1f}" for val in b['StokesI_unstretched_cbar_ticks_grid']])
                
                
            cbar1.ax.get_xticklabels()[0].set_ha('left')
            cbar1.ax.get_xticklabels()[-1].set_ha('right')

        # Put ticks on top
        cbar1.ax.xaxis.set_ticks_position("top")

        # Format colorbar
        cbar1.ax.tick_params(
            labelsize=writeup_grid_axis_num_fs,
            which="major",
            length=4,
            direction="in",
            width=spine_width
        )
        
#         if plotting == 'Stokes I stretched':
#             if i != 0:
#                 cbar1.ax.get_xticklabels()[0].set_visible(False)
        
        if plotting == 'POLI':
            if i == 2:
                for tick, label_cbar1 in zip(cbar1.get_ticks(), cbar1.ax.get_xticklabels()):
                    if tick == 0:
                        label_cbar1.set_visible(False)
                    
                    
        if plotting == 'Gaussian Map':
            if i != 3:
                for tick, label_cbar1 in zip(cbar1.get_ticks(), cbar1.ax.get_xticklabels()):
                    if tick == 1.0:
                        label_cbar1.set_visible(False)

        cbar1.outline.set_linewidth(spine_width)
        
        
        #cbar1 = plt.colorbar(im, cax = cax, orientation='horizontal', shrink=cbar_shrink, location = 'top')
        # cbar1.set_label('POLI (mJy beam$^{-1}$)', fontsize=cbar_fs)
#         cbar1.ax.tick_params(labelsize=axis_num_fs, which='major', length=4, direction="in", width=spine_width,)
#         cbar1.outline.set_linewidth(spine_width)
        #cbar1.set_label("Polarized Intensity [mJy beam$^{-1}$]",fontsize=cbar_label_fs)
        
        # Add the vectors
        j = 0 
        if show_vectors:
            for j, row in enumerate(b['VectorsActual']):
                if j == 0:
                    ax.plot([row[0], row[1]], [row[2], row[3]], color='black', lw = vector_lw, label = 'Observations')
                else:
                    ax.plot([row[0], row[1]], [row[2], row[3]], color='black', lw = vector_lw)

        #centre_pix_x, centre_pix_y = ax.world_to_pixel(centre)

        half_size_pix = 0.9 / 0.018

        xmin = b['RA_centre_pix'] - half_size_pix
        xmax = b['RA_centre_pix'] + half_size_pix
        ymin = b['Dec_centre_pix'] - half_size_pix
        ymax = b['Dec_centre_pix'] + half_size_pix
        
#         print(f'Band = {i}, and ymax = {ymax}')
    

        # Beam location
        beam_x_pos = xmin - 0.1 * (xmin - xmax)
        beam_y_pos = ymin - 0.1 * (ymin - ymax)

        
        # Add y-axis label to the left
#         if for_poster:
#             if i == 0:
#                 ax.set_ylabel('$\Delta$ RA (arcsec)', fontsize=axis_label_fs)
        
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

        if plotting != "Gaussian Map":
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
        if plotting != "Gaussian Map":
            ax.text(
                x_axes,
                y_axes + 0.015,
                f"{b['reference_length_AU']} au",
                transform=ax.transAxes,
                fontsize=AU100_fs,
                ha='center',
                va='bottom'
            )
            
        print(rf'AU 100 text location: x_axes = {x_axes} and y_axes + 0.015 = {y_axes + 0.015}')
        
        add_band_label(ax, bands_idx[i], label, constants, x_pos = x_pos, y_pos = y_pos, 
                       label_fs = 25, va = 'bottom')
        # ---------------------------------------------------------
    
    
        # Add the beam
        beam = Ellipse((beam_x_pos, beam_y_pos), width=b['BMAJ_pix'], height=b['BMIN_pix'],  
                   angle=b['BPA_deg_cartesian'], edgecolor='black', facecolor='none', alpha=1, lw=beam_lw)
        
        if plotting != "Gaussian Map":
            ax.add_patch(beam)
        
        # Hide the original RA/Dec WCS coordinates
        # ---------------------------------------------------------
        ax.coords[0].set_ticks_visible(False)
        ax.coords[0].set_ticklabel_visible(False)
        ax.coords[0].set_axislabel('')

        ax.coords[1].set_ticks_visible(False)
        ax.coords[1].set_ticklabel_visible(False)
        ax.coords[1].set_axislabel('')
        # --------------------------------------------------------


        # Format the offset coordinates
        # --------------------------------------------------------
        # Put Delta RA ticks/labels on the bottom
        overlay[0].set_ticks_position('b')
        overlay[0].set_ticklabel_position('b')
        overlay[0].set_axislabel_position('b')

        # Put Delta Dec ticks/labels on the left
        overlay[1].set_ticks_position('l')
        overlay[1].set_ticklabel_position('l')
        overlay[1].set_axislabel_position('l')

        overlay[0].set_format_unit(u.arcsec)
        overlay[1].set_format_unit(u.arcsec)

        overlay[0].set_axislabel(r' ', fontsize=0)
        overlay[1].set_axislabel(r' ', fontsize=0)

        overlay[0].set_ticks(spacing=0.8* u.arcsec, direction='in', width=spine_width)
        overlay[1].set_ticks(spacing=0.8 * u.arcsec, direction='in', width=spine_width)

        overlay[0].set_ticklabel(size=writeup_grid_axis_num_fs)
        overlay[1].set_ticklabel(size=writeup_grid_axis_num_fs)

        if i != 0:
            overlay[1].set_ticklabel_visible(False)
        

        # Put Delta RA ticks/labels on the bottom
#         x_ticks = np.arange(-0.6, 0, 0.6) * u.arcsec
#         y_ticks = np.arange(-0.4, 0, 0.4) * u.arcsec

#         overlay[0].set_ticks(values=x_ticks)
#         overlay[1].set_ticks(values=y_ticks)
        # --------------------------------------------------------
        
    
#     axes[0].set_ylabel('$\Delta$ Dec (arcsec)', fontsize=axis_label_fs)
    
    #fig.subplots_adjust(wspace=wspace)
   
    format_writeup_grid(fig, cbar_label)
    
    
    return fig, axes














# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------
def scale_factor_plot_for_writeup(bands,                  # These are the bands we are looking at
                                        bands_labels, 
                                        bands_included_in_fit, # These are the values that were fit
                                        f_values,            # These are the porosity values we are testing
                                        a_max_dist_micron_by_f,
                                        P_times_omega_by_f,
                                        sf_by_f,             # This is the values the data is being scaled by 
                                        best_idx_by_f, 
                                        best_a_max_by_f,     # These are the best a_max values for each f
                                        POLF_markers,        # These are the values for each band we are plotting the marker
                                        chi_sq_by_f, 
                                        ymin = -0.1, 
                                        ymax = 2, 
                                        xmin = 15, xmax = 1e3 + 20, 
                                        custom_lw = None,
                                        custom_text_x = None,
                                        chi_sq_precision = 2,
                                        fig_x = 30, fig_y = 10,
                                        f_fs = 45, 
                                        f_x = 0.05, f_y = 0.86,
                                        spine_width = 1,
                                        included_lw = 5,
                                        marker_border = False,
                                        BandLegend_bbox_to_anchor=(0.5, 0.85),
                                        BandLegend_handletextpad=0.5,
                                        BandLegend_handlelength=1,
                                        BandLegend_labelspacing=0.6,
                                        plot_markers = True,
                                        plot_band_legend = True,
                                        bands_lw = 5, 
                                        num_fs = 35,
                                        xy_axis_fs = 40,
                                        stats_fs = 35,
                                        legend_marker_size = 20,
                                        band_legend_fs = 32.5,
                                        marker_size = 600,
                                        plot_sf = False,
                                        plot_chi_sq = False):
    
    
    # Set the colors and marker size for each band 
    band_colors = []
    ms = []

    for b in bands:
        if b == "Band 5 robust -1":
            band_colors.append(alma_band_colors["Band 5"])
            ms.append(constants.alma_band_ms["Band 5"])
        else:
            band_colors.append(alma_band_colors[b])
            ms.append(constants.alma_band_ms[b])
    
    
    
    # Make the figure
    # Later I will make this so you can customize how many figures you need
    fig, ax = plt.subplots(1, 4, figsize=(fig_x, fig_y))
    
    
    
    
    # Loop over the f values you want to plot
    for i, f in enumerate(f_values): 
        
        for spine in ax[i].spines.values():
            spine.set_linewidth(spine_width)
        
       

        if i == 0:
            ax[i].set_ylabel(r'$P\omega_{\mathrm{eff}}$', fontsize = xy_axis_fs)
        else:
            ax[i].set_yticklabels([])
            ax[i].set_ylabel('')
            ax[i].tick_params(left=False)

        # Make the x-axis log scale 
        ax[i].set_xscale("log")
        ax[i].set_xlim(xmin, xmax)


            
        fig.supxlabel('Maximum grain size [$\mu$m]', fontsize = xy_axis_fs)
        ax[i].set_title(f'$f$ = {f:.2f}', fontsize=45)

            
    
        # Set ticks
        ax[i].minorticks_on()
        ax[i].tick_params(axis="x", which="minor", direction="in", left=True, right=True, length=4, width=spine_width)
        ax[i].tick_params(axis="y", which="minor", direction="in", left=True, right=True, length=4, width=spine_width)
        
        
        # Control the font size of the numbers on the axes
        # ---------------------------------------------------------------------------------------
            
        ax[i].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize = num_fs, width=spine_width)
        ax[i].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize = num_fs, width=spine_width)
        # ---------------------------------------------------------------------------------------
        
        
        # ---------------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------------


           

            
        
#         if for_poster == True:
#             ax[i].text(x_pos, 0.15,  f'$\chi^2$ = {chi_sq:.{chi_sq_precision}f}', transform=ax[i].transAxes, fontsize = 30)
        # ---------------------------------------------------------------------------------------

        ax[i].set_ylim(ymin, ymax)
        # ------------------------------------------
        
        
        marker_handles = []
        # plot each band for this f
        for j, b in enumerate(bands):

            idx = bands.index(b)


#             print(rf'At band {b}, the POLF_obs value is: {POLF_obs[f][idx]}')
#             print(' ')

            fc = band_colors[j]

            if marker_border:
                ec = 'black'
                marker_lw = 1.5
            else:
                ec = band_colors[j]
                marker_lw = 0
            #ec = band_colors[j]
#             marker_lw = 0 

#             if for_poster and b == "Band 5 robust -1":
#                 fc = alma_band_colors["Band 5"]
#                 ec = alma_band_colors["Band 5"]

            if b == "Band 5" and not for_poster:
                fc = 'none'
                marker_lw = 5

            # Plot the POLF markers
            plot_marker = True

            if b == "Band 5":
                plot_marker = False

            if plot_markers:
                ax[i].scatter(
                    best_a_max_by_f[f],
                    POLF_markers[b],
                    marker=ms[j], 
                    s=marker_size,
                    facecolors=fc,
                    edgecolors=ec,
                    linewidths=marker_lw,
                    zorder = 5
                )
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

       # if b not in ["Band 5 robust -1", "Band 5 robust -2"]:
            ax[i].plot(a_max_dist_micron_by_f[f],
                           P_times_omega_by_f[f][:, idx] * sf_by_f[f],
                           color=band_colors[j],
                           label=b,
                           lw=bands_lw,
                           ls='-' if b in bands_included_in_fit else '--',
                           zorder = 1)
        
                
        x_pos = 0.05 
        
        if custom_text_x is not None:
            x_pos = x_pos + custom_text_x[i]
        
        if plot_chi_sq:
            ax[i].text(x_pos, 0.8,  f'$\chi^2$ = {chi_sq_by_f[f]:.{chi_sq_precision}f}', transform=ax[i].transAxes, fontsize = stats_fs)

        if plot_sf:
            ax[i].text(x_pos, 0.7, f'sf = {sf_by_f[f]:.2f}', transform=ax[i].transAxes, fontsize = stats_fs)


        # Add the a_max to the text
        ax[i].text(x_pos, 0.9, rf'$a_{{\mathrm{{max}}}}$ = {best_a_max_by_f[f]:.0f} $\mu$m', 
                   transform=ax[i].transAxes, 
                   fontsize = stats_fs)
        
        
        
    #plt.tight_layout()
    
    plt.subplots_adjust(
        left=0.05,
        right=0.98,
        bottom=0.15,
        top=0.95,
        wspace=0.02
    )
    
    
    # Add coloured text for each Band in the far right plot
    # ----------------------------------------------------
    # ----------------------------------------------------
    #start = start # This is the starting vertical position
    c = 0
    marker_legend_handles = []
    
    # Choose legend entries
    legend_bands = ["Band 4", "Band 5 robust -1", "Band 6", "Band 7"]

    for b_label in legend_bands:

        if b_label == "Band 5 robust -1":
            label = lambda_mm["Band 5"]
            color = alma_band_colors["Band 5"]
            marker = constants.alma_band_ms["Band 5"]
        else:
            label = lambda_mm[b_label]
            color = alma_band_colors[b_label]
            marker = constants.alma_band_ms[b_label]


        handle = Line2D(
            [0], [0],
            marker=marker,
            color='none',
            markerfacecolor=color,
            markeredgecolor=color,
            markersize=legend_marker_size,
            linestyle='None',
            label=f"{label} mm"
        )

        marker_legend_handles.append(handle)
        
        
        
    legend_elements = [
        Line2D([0], [0], color='black', linestyle='-', lw = 5, label='Included'),
        Line2D([0], [0], color='black', linestyle='--', lw = 5, label='Excluded')
    ]
    # ----------------------------------------------------
    # ----------------------------------------------------
    #if plot_markers: 
    if plot_band_legend:
        legend = ax[3].legend(
                handles=marker_legend_handles,
                fontsize= band_legend_fs,
                frameon=False,
                loc='upper right',
                bbox_to_anchor=BandLegend_bbox_to_anchor,
                handletextpad=BandLegend_handletextpad, # Cotrols distance between marker and text 
                handlelength=BandLegend_handlelength, # Cotrols distance between marker and text 
                labelspacing=BandLegend_labelspacing # Controls vertical gap between the markers
            )

        legend_colors = [
            alma_band_colors["Band 4"],
            alma_band_colors["Band 5"],
            alma_band_colors["Band 6"],
            alma_band_colors["Band 7"]
        ]

        for text, color in zip(legend.get_texts(), legend_colors):
            text.set_color(color)
#     # ----------------------------------------------------
#     # ----------------------------------------------------



# 

    return fig, ax
# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------





def scale_factor_plot_for_writeup_glitterin(bands,                  # These are the bands we are looking at
                                        bands_labels, 
                                            results_array, 
                                        bands_included_in_fit_arr, # These are the values that were fit
                                        df_POLF,        # These are the values for each band we are plotting the marker
                                        albedo = 'w',
                                        #chi_sq_by_f, 
                                        ymin = -0.1, 
                                        ymax = 2, 
                                        xmin = 15, xmax = 1e3 + 20, 
                                        custom_lw = None,
                                        custom_text_x = None,
                                        chi_sq_precision = 2,
                                        fig_x = 30, fig_y = 8,
                                        f_fs = 45, 
                                        f_x = 0.05, f_y = 0.86,
                                        spine_width = 1,
                                        included_lw = 5,
                                        marker_border = False,
                                        BandLegend_bbox_to_anchor=(0.425, 0.9),
                                        BandLegend_handletextpad=0.5,
                                        BandLegend_handlelength=1,
                                        BandLegend_labelspacing=0.6,
                                        plot_markers = True,
                                        plot_band_legend = True,
                                        bands_lw = 5, 
                                        num_fs = 35,
                                        xy_axis_fs = 40,
                                        stats_fs = 35,
                                        legend_marker_size = 20,
                                        band_legend_fs = 32.5,
                                        marker_size = 600,
                                        plot_sf = False,
                                        plot_chi_sq = False,
                                           legend_loc = 3,
                                           titles = ['All wavelengths', 'No 1.5 mm', 'No 1.3 mm', 'Just 2.1 mm and 0.87 mm'],
                                           x_axis_pad = 10):
    
    POLF_markers = dict(zip(df_POLF["Band"], df_POLF["POLF_maxStokesI"]))
    
    # Set the colors and marker size for each band 
    band_colors = []
    ms = []

    for b in bands:
        if b == "Band 5 robust -1":
            band_colors.append(alma_band_colors["Band 5"])
            ms.append(constants.alma_band_ms["Band 5"])
        else:
            band_colors.append(alma_band_colors[b])
            ms.append(constants.alma_band_ms[b])
    
    
    num_fits = len(results_array)
    
    # Make the figure
    # Later I will make this so you can customize how many figures you need
    fig, ax = plt.subplots(1, num_fits, figsize=(fig_x, fig_y))
    
    
    
    
    # Loop over the f values you want to plot
#     titles = ['All wavelengths', 'No 1.5 mm', 'No 1.3 mm', 'Just 2.1 mm and 0.87 mm']
    
    
    
    for i in range(num_fits):
        
        bands_included_in_fit = bands_included_in_fit_arr[i]
        
        res = results_array[i]
        a_max_dist_micron = cm_to_micron(res['rvol_max_cm'])
        
        if albedo == 'w_eff':
            P_times_omega = res['Pw_eff']
            y_label = '$P\omega_{\mathrm{eff}}$'
        elif albedo == 'w':
            P_times_omega = res['Pw']
            y_label = '$P\omega$'
        else:
            return("albedo should equal w or w_eff")
        
        
        sf = res['best_sf']
        best_a_max = res['best_rvol_max_micron']
        best_idx = res['best_idx_arr']
        chi_sq = res['chi_sq']
        
        for spine in ax[i].spines.values():
            spine.set_linewidth(spine_width)
        
       

        if i == 0:
            ax[i].set_ylabel(y_label, fontsize = xy_axis_fs)
        else:
            ax[i].set_yticklabels([])
            ax[i].set_ylabel('')
            ax[i].tick_params(left=False)

        # Make the x-axis log scale 
        ax[i].set_xscale("log")
        ax[i].set_xlim(xmin, xmax)


            
        fig.supxlabel('Maximum grain size [$\mu$m]', fontsize = xy_axis_fs)
        ax[i].set_title(f'{titles[i]}', fontsize=40, pad = 15)

            
    
        # Set ticks
        ax[i].minorticks_on()
        ax[i].tick_params(axis="x", which="minor", direction="in", left=True, right=True, length=4, width=spine_width)
        ax[i].tick_params(axis="y", which="minor", direction="in", left=True, right=True, length=4, width=spine_width)
        
        
        # Control the font size of the numbers on the axes
        # ---------------------------------------------------------------------------------------
            
        ax[i].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize = num_fs, width=spine_width, pad = x_axis_pad)
        ax[i].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize = num_fs, width=spine_width)
        
        
        # Hide the 10^3 label on the left and middle panels
#         if i in [0, 1]:
#             for tick, label in zip(ax[i].get_xticks(), ax[i].get_xticklabels()):
#                 if np.isclose(tick, 1e3):
#                     label.set_visible(False)
        # ---------------------------------------------------------------------------------------
        
        
        # ---------------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------------


           

            
        
#         if for_poster == True:
#             ax[i].text(x_pos, 0.15,  f'$\chi^2$ = {chi_sq:.{chi_sq_precision}f}', transform=ax[i].transAxes, fontsize = 30)
#         # ---------------------------------------------------------------------------------------

        ax[i].set_ylim(ymin, ymax)
        # ------------------------------------------
        
        
        marker_handles = []
        # plot each band for this f
        for j, b in enumerate(bands):

            idx = bands.index(b)


#             print(rf'At band {b}, the POLF_obs value is: {POLF_obs[f][idx]}')
#             print(' ')

            fc = band_colors[j]

            if marker_border:
                ec = 'black'
                marker_lw = 1.5
            else:
                ec = band_colors[j]
                marker_lw = 0
            #ec = band_colors[j]
#             marker_lw = 0 

#             if for_poster and b == "Band 5 robust -1":
#                 fc = alma_band_colors["Band 5"]
#                 ec = alma_band_colors["Band 5"]

            if b == "Band 5" and not for_poster:
                fc = 'none'
                marker_lw = 5

            # Plot the POLF markers
            plot_marker = True

            if b == "Band 5":
                plot_marker = False

            if plot_markers:
                ax[i].scatter(
                    best_a_max,
                    POLF_markers[b],
                    marker=ms[j], 
                    s=marker_size,
                    facecolors=fc,
                    edgecolors=ec,
                    linewidths=marker_lw,
                    zorder = 5
                )
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

       # if b not in ["Band 5 robust -1", "Band 5 robust -2"]:
            print(rf'i = {i} and P_times_omega[:, j] * sf max is {max(P_times_omega[:, j] * sf) }')
            ax[i].plot(a_max_dist_micron,
                           P_times_omega[:, j] * sf,
                           color=band_colors[j],
                           label=b,
                           lw=bands_lw,
                           ls='-' if b in bands_included_in_fit else '--',
                           zorder = 1)
        
                
        x_pos = 0.05 
        
        if custom_text_x is not None:
            x_pos = x_pos + custom_text_x[i]
        
        if plot_chi_sq:
            ax[i].text(x_pos, 0.8,  f'$\chi^2$ = {chi_sq:.{chi_sq_precision}f}', transform=ax[i].transAxes, fontsize = stats_fs)

        if plot_sf:
            ax[i].text(x_pos, 0.7, f'sf = {sf:.2f}', transform=ax[i].transAxes, fontsize = stats_fs)


        # Add the a_max to the text
        ax[i].text(x_pos, 0.9, rf'$a_{{\mathrm{{max}}}}$ = {best_a_max:.0f} $\mu$m', 
                   transform=ax[i].transAxes, 
                   fontsize = stats_fs)
        
        
        
    #plt.tight_layout()
    
    plt.subplots_adjust(
        left=0.05,
        right=0.98,
        bottom=0.15,
        top=0.95,
        wspace=0.02
    )
    
    
    # Add coloured text for each Band in the far right plot
    # ----------------------------------------------------
    # ----------------------------------------------------
    #start = start # This is the starting vertical position
    c = 0
    marker_legend_handles = []
    
    # Choose legend entries
    legend_bands = ["Band 4", "Band 5 robust -1", "Band 6", "Band 7"]

    for b_label in legend_bands:

        if b_label == "Band 5 robust -1":
            label = constants.lambda_mm["Band 5"]
            color = alma_band_colors["Band 5"]
            marker = constants.alma_band_ms["Band 5"]
        else:
            label = constants.lambda_mm[b_label]
            color = alma_band_colors[b_label]
            marker = constants.alma_band_ms[b_label]


        handle = Line2D(
            [0], [0],
            marker=marker,
            color='none',
            markerfacecolor=color,
            markeredgecolor=color,
            markersize=legend_marker_size,
            linestyle='None',
            label=f"{label} mm"
        )

        marker_legend_handles.append(handle)
        
        
        
    legend_elements = [
        Line2D([0], [0], color='black', linestyle='-', lw = 5, label='Included'),
        Line2D([0], [0], color='black', linestyle='--', lw = 5, label='Excluded')
    ]
    # ----------------------------------------------------
    # ----------------------------------------------------
    #if plot_markers: 
    if plot_band_legend:
        legend = ax[legend_loc].legend(
                handles=marker_legend_handles,
                fontsize= band_legend_fs,
                frameon=False,
                loc='upper right',
                bbox_to_anchor=BandLegend_bbox_to_anchor,
                handletextpad=BandLegend_handletextpad, # Cotrols distance between marker and text 
                handlelength=BandLegend_handlelength, # Cotrols distance between marker and text 
                labelspacing=BandLegend_labelspacing # Controls vertical gap between the markers
            )

        legend_colors = [
            alma_band_colors["Band 4"],
            alma_band_colors["Band 5"],
            alma_band_colors["Band 6"],
            alma_band_colors["Band 7"]
        ]

        for text, color in zip(legend.get_texts(), legend_colors):
            text.set_color(color)
#     # ----------------------------------------------------
#     # ----------------------------------------------------



# 

    return fig, ax
# -----------------------------------------------------------------------------------


def plot_Pw_vs_grainsize(data, bands, lambda_bands_cm, bands_labels = ['B4', 'B5', 'B6', 'B7'],
                         
                         porosity='1', albedo='w',
                        idx = 'wav'):
    
    # Choose which P*omega quantity to plot
    if albedo == 'w':
        Pw_key = 'P_omega'
        ylabel = r'$P\omega$'
        save_name = f'Pw_vs_grainsize_f{porosity}.pdf'
        
    elif albedo == 'w_eff':
        Pw_key = 'P_omega_eff'
        ylabel = r'$P\omega_{\mathrm{eff}}$'
        save_name = f'Pw_eff_vs_grainsize_f{porosity}.pdf'
        
    else:
        raise ValueError("albedo must be 'w' or 'w_eff'")
    
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    wavelengths = [lambda_mm[band] for band in bands]
    nums = [4, 5, 6, 7]
    for i, band in enumerate(bands):
        
        if idx == 'wav':
            lam = wavelengths[i]
        elif idx == 'num':
            lam = nums[i]
        
        ax.plot(
            cm_to_micron(data[lam][porosity]['a']),
            data[lam][porosity][Pw_key],
            color=alma_band_colors[band],
            ls='-',
            lw=2,
            label=f'{cm_to_mm(lambda_bands_cm[i]):.1f} mm'
                  if i != 3
                  else f'{cm_to_mm(lambda_bands_cm[i]):.2f} mm'
        )
        
        # Characteristic grain size
        a_char_micron = cm_to_micron(lambda_bands_cm[i] / (2 * np.pi))
        
        ax.axvline(
            a_char_micron,
            color=alma_band_colors[band],
            alpha=0.25,
            lw=2,
            ls='--'
        )
        
        # Find maximum
        max_idx = np.argmax(data[lam][porosity][Pw_key])
        max_val = data[lam][porosity][Pw_key][max_idx]
        amax_at_max = cm_to_micron(data[lam][porosity]['a'][max_idx])

        print(f'At Band {bands_labels[i]}, the highest Pw value is {max_val}')
        print(f'At this Pw value, the a_max value is {amax_at_max:.2f} micron')
        print(f'a_char_micron = {a_char_micron:.2f} micron')
        print()
    
    
    # Axis labels
    ax.set_xlabel('Maximum grain size [$\mu$m]', fontsize=axis_label_fs)
    ax.set_ylabel(ylabel, fontsize=axis_label_fs)
    
    # Ticks
    ax.minorticks_on()

    ax.tick_params(axis="x", which="minor",
                   direction="in", left=True, right=True, length=4)
    ax.tick_params(axis="y", which="minor",
                   direction="in", left=True, right=True, length=4)

    ax.tick_params(axis="x", which="major",
                   direction="in", bottom=True, top=True,
                   length=7, labelsize=axis_num_fs)
    ax.tick_params(axis="y", which="major",
                   direction="in", left=True, right=True,
                   length=7, labelsize=axis_num_fs)

    # Log x-axis
    ax.set_xscale('log')

    # Remove scientific notation
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
    ax.ticklabel_format(style='plain', axis='x')

    # Characteristic grain size legend entry
    ax.axvline(
        -100,
        color='black',
        alpha=0.25,
        lw=2,
        ls='--',
        label=r'$a_{\mathrm{char}} = \lambda / 2\pi$'
    )

    ax.set_xlim(10, 1000)

    ax.legend(fontsize=legend_text_fs, frameon=False)

#     plt.savefig(
#         writeup_image_folder_path + '/DustModels/' + save_name,
#         dpi=300,
#         facecolor='white',
#         bbox_inches='tight'
#     )

    return fig, ax
# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------