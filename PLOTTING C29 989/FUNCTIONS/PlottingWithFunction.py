# Import necessary packages
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns
import math

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
# -----------------------------------------------------------------------------------------



# Function to make StokesI normalized colorbar ticks 
# -----------------------------------------------------------------------------------------
def normalize_stokesI_for_cmap(StokesI_data_2d_mJy, custom_min=None, custom_max=None, base=100):
    """
    Normalize Stokes I data for visualization and generate colorbar tick values.

    Parameters
    ----------
    StokesI_data_2d_mJy : ndarray
        2D array of Stokes I intensity values in mJy/beam.
    normalized_cbar_ticks : array-like
        Tick positions on the normalized color scale (e.g., [0.2, 0.4, 0.6, 0.8]).
    custom_min : float, optional
        Custom minimum value to use for clipping and stretching. If None, uses np.nanmin.
    custom_max : float, optional
        Custom maximum value to use for clipping and stretching. If None, uses np.nanmax.
    base : float, optional
        Stretching base used in the normalization (default is 100).

    Returns
    -------
    StokesI_stretched : ndarray
        The normalized and stretched version of the input Stokes I data.
    StokesI_unstretched_cbar_ticks : ndarray
        The corresponding physical intensity values (in mJy/beam) for the colorbar ticks.
    """

    # Use custom min/max if provided, else fallback to nan-safe defaults
    StokesI_min = custom_min if custom_min is not None else np.nanmin(StokesI_data_2d_mJy)
    StokesI_max = custom_max if custom_max is not None else np.nanmax(StokesI_data_2d_mJy)

    # Clip the data to the selected range
    StokesI_clipped = np.clip(StokesI_data_2d_mJy, StokesI_min, StokesI_max)

    # Apply stretch transformation
    StokesI_stretched = stretch(StokesI_clipped, base=base, vmin=StokesI_min, vmax=StokesI_max)

    # Map normalized colorbar ticks back to original units
    StokesI_unstretched_cbar_ticks = unstretch(normalized_cbar_ticks, StokesI_min, StokesI_max)

    return StokesI_stretched, StokesI_unstretched_cbar_ticks
# -----------------------------------------------------------------------------------------




# ---------------------------------------------------------------------------------------------------------------
def create_stokes_i_base_plot(StokesI_wcs, StokesI_stretched, cmap, 
                              normalized_cbar_ticks, StokesI_unstretched_cbar_ticks, 
                              xmin, xmax, ymin, ymax, reference_length_pix, reference_length_AU,
                              BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
                              max_length_pix, reference_fraction,
                              text_fs = constants.text_fs , 
                              axis_label_fs = constants.axis_label_fs, 
                              axis_num_fs = constants.axis_num_fs, 
                              cbar_fs = constants.cbar_fs):
    
    # Create a figure with the WCS projection
    fig, ax = plt.subplots(figsize=(14, 12), subplot_kw={'projection': StokesI_wcs})

    # Add data
    im = ax.imshow(StokesI_stretched, cmap=cmap)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r'Stokes I (mJy/beam)', fontsize=cbar_fs)
    cbar.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")
    cbar.ax.tick_params(which='minor', length=4, direction="in")
    cbar.set_ticks(normalized_cbar_ticks)
    cbar.set_ticklabels([f"{val:.2f}" for val in StokesI_unstretched_cbar_ticks])

    # Add axis labels
    ax.set_xlabel('Right Ascension', fontsize=axis_label_fs)
    ax.set_ylabel('Declination', fontsize=axis_label_fs)

    # Set x and y limits
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # Add line and text for 100 AU 
    line_x_pos = xmax - 0.05 * (xmax - xmin) 
    line_y_pos = ymax - 0.1  * (ymax - ymin) 

    ax.plot([(line_x_pos - reference_length_pix), (line_x_pos)], 
            [line_y_pos, line_y_pos],
            color='black', linewidth=3)

    ax.text((line_x_pos - reference_length_pix/2), (line_y_pos - 2), 
            f'{reference_length_AU} AU', fontsize=text_fs, ha='center', va='top') 

    # Create and add the circular beam
    beam_x_pos = xmin - 0.1 * (xmin - xmax)
    beam_y_pos = ymin - 0.1 * (ymin - ymax)

    beam = Ellipse((beam_x_pos, beam_y_pos), width=BMAJ_pix, height=BMIN_pix,  
                   angle=BPA_deg_cartesian, edgecolor='black', facecolor='none', alpha=1, lw=2)

    ax.add_patch(beam)

    # Adjust ticks
    ax.minorticks_on()
    ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
    ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)

    return fig, ax
# ---------------------------------------------------------------------------------------------------------------


def create_base_plot(StokesI_wcs, plotting_data, cbar_label, cmap, 
                     xmin, xmax, ymin, ymax, reference_length_pix, reference_length_AU,
                     BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
                     max_length_pix, reference_fraction,
                     text_fs = constants.text_fs , 
                     axis_label_fs = constants.axis_label_fs, 
                     axis_num_fs = constants.axis_num_fs, 
                     cbar_fs = constants.cbar_fs):
    
    # Create a figure with the WCS projection
    fig, ax = plt.subplots(figsize=(14, 12), subplot_kw={'projection': StokesI_wcs})

    # Add data
    im = ax.imshow(plotting_data, cmap=cmap)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(cbar_label, fontsize=cbar_fs)
    cbar.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")
    cbar.ax.tick_params(which='minor', length=4, direction="in")
    #cbar.set_ticks(normalized_cbar_ticks)
    #cbar.set_ticklabels([f"{val:.2f}" for val in StokesI_unstretched_cbar_ticks])

    # Add axis labels
    ax.set_xlabel('Right Ascension', fontsize=axis_label_fs)
    ax.set_ylabel('Declination', fontsize=axis_label_fs)

    # Set x and y limits
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # Add line and text for 100 AU 
    line_x_pos = xmax - 0.05 * (xmax - xmin) 
    line_y_pos = ymax - 0.1  * (ymax - ymin) 

    ax.plot([(line_x_pos - reference_length_pix), (line_x_pos)], 
            [line_y_pos, line_y_pos],
            color='black', linewidth=3)

    ax.text((line_x_pos - reference_length_pix/2), (line_y_pos - 2), 
            f'{reference_length_AU} AU', fontsize=text_fs, ha='center', va='top') 

    # Create and add the circular beam
    beam_x_pos = xmin - 0.1 * (xmin - xmax)
    beam_y_pos = ymin - 0.1 * (ymin - ymax)

    beam = Ellipse((beam_x_pos, beam_y_pos), width=BMAJ_pix, height=BMIN_pix,  
                   angle=BPA_deg_cartesian, edgecolor='black', facecolor='none', alpha=1, lw=2)

    ax.add_patch(beam)

    # Adjust ticks
    ax.minorticks_on()
    ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
    ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)

    return fig, ax
# ---------------------------------------------------------------------------------------------


def create_blank_grid(i, j, ax,
                      data_plotting, cbar_label, 
                      StokesI_wcs, soft_colormap_v2, 
                      xmin, xmax, ymin, ymax, reference_length_pix, reference_length_AU,
                      BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
                      max_length_pix, reference_fraction,
                      fs_scale = 0.5):
    
    
    # Plot the data on the given ax
    im = ax.imshow(data_plotting, cmap=soft_colormap_v2)

    # Colorbar for every rightmost column (in a 6x2 grid) to avoid redundancy
    if j == 1:
        cbar = plt.colorbar(im, ax=ax, pad=0.02, fraction=0.05)
        cbar.set_label(cbar_label, fontsize= cbar_fs)
        cbar.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")
        cbar.ax.tick_params(which='minor', length=4, direction="in")

    # Set axis labels
    ax.set_xlabel('Right Ascension', fontsize=axis_label_fs * fs_scale)
    ax.set_ylabel('Declination', fontsize=axis_label_fs * fs_scale)

    # Set x and y limits
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # Add reference scale line
    line_x_pos = xmax - 0.05 * (xmax - xmin) 
    line_y_pos = ymax - 0.1  * (ymax - ymin) 

    ax.plot([(line_x_pos - reference_length_pix), (line_x_pos)], 
            [line_y_pos, line_y_pos],
            color='black', linewidth=3)

    ax.text((line_x_pos - reference_length_pix/2), (line_y_pos - 2), 
            f'{reference_length_AU} AU', fontsize=text_fs * fs_scale, ha='center', va='top') 

    # Create and add the circular beam
    beam_x_pos = xmin - 0.1 * (xmin - xmax)
    beam_y_pos = ymin - 0.1 * (ymin - ymax)

    beam = Ellipse((beam_x_pos, beam_y_pos), width=BMAJ_pix, height=BMIN_pix,  
                   angle=BPA_deg_cartesian, edgecolor='black', facecolor='none', alpha=1, lw=2)

    ax.add_patch(beam)

#     # Draw the reference vector line
#     vector_x_pos = xmax - 0.05 * (xmax - xmin)
#     vector_y_pos = ymin - 0.1 * (ymin - ymax)

#     ax.plot([(vector_x_pos - max_length_pix * reference_fraction), vector_x_pos], 
#             [vector_y_pos, vector_y_pos], color='black', linewidth=3)

#     ax.text(vector_x_pos - max_length_pix * reference_fraction / 2,  
#             vector_y_pos - 2,  
#             f'{reference_fraction * 100:.0f}%',  
#             fontsize=text_fs, ha='center', va='top')

    # Adjust ticks
    ax.minorticks_on()
    ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs * fs_scale)
    ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs* fs_scale)

#     # Optional: Title with the row/column info
#     ax.set_title(f"Plot {i+1},{j+1}", fontsize=10)

# ---------------------------------------------------------------------------------------------
  

# ---------------------------------------------------------------------------------------------
def plot_grids(data_list, soft_colormap_v2, StokesI_wcs, axis_label_fs, axis_num_fs, cbar_fs, xmin, xmax, ymin, ymax, BMAJ_pix, BMIN_pix, BPA_deg_cartesian):
    grid_titles = ['100% U', '50% U 50% A', '100% A']
    cb_titles = [
        "Polarization Intensity", "Polarization Intensity", "Polarization Intensity",
        "Stokes Q ", "Stokes Q", "Stokes Q",
        "Stokes U", "Stokes U", "Stokes U"
    ]

    top_cb_pad = 0.1
    middle_cb_pad = 0.1
    bottom_cb_pad = 0.2

    cb_pads = [
        top_cb_pad , top_cb_pad , top_cb_pad ,
        middle_cb_pad, middle_cb_pad, middle_cb_pad,
        bottom_cb_pad, bottom_cb_pad, bottom_cb_pad
    ]

    POLA_cmap = soft_colormap_v2
    StokesQ_cmap = soft_colormap_v2
    StokesU_cmap = soft_colormap_v2

    cb_cmaps = [
        POLA_cmap, POLA_cmap, POLA_cmap, 
        StokesQ_cmap, StokesQ_cmap, StokesQ_cmap, 
        StokesU_cmap, StokesU_cmap, StokesU_cmap   
    ]

    # Create figure with 3x3 grid of subplots
    fig, axes = plt.subplots(3, 3, figsize=(18, 18), subplot_kw={'projection': StokesI_wcs},
                             gridspec_kw={'hspace': 0.2, 'wspace': -0.3})  # Remove spacing

    # Loop through each subplot and plot corresponding data
    for i, ax in enumerate(axes.flat):
        # Plot the data 
        im = ax.imshow(data_list[i], cmap=cb_cmaps[i])

        # Set titles for the top row only
        row, col = divmod(i, 3)  # Calculate the row and column index
        if row == 0:
            ax.set_title(grid_titles[i], fontsize=axis_label_fs)

        # Y-axis 
        # ---------------------------------------------------------------------------------------------
        if col == 0:
            ax.set_ylabel('Dec', fontsize=axis_label_fs)
            ax.tick_params(axis="y", which="both", left=True, labelleft=True)  # Show ticks and labels on left
        else:
            ax.tick_params(axis="y", which="both", left=True, labelleft=False)  # Show ticks but no labels on middle and right
            ax.set_ylabel("")  # No label on middle or right column
        # ---------------------------------------------------------------------------------------------
        
        # X-axis 
        # ---------------------------------------------------------------------------------------------
        if row == 2:
            ax.set_xlabel('RA', fontsize=axis_label_fs)
            ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)
        else:
            ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
            ax.set_xlabel("")  # No label on top or middle row
        # ---------------------------------------------------------------------------------------------

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        
        # Create and add the circular beam
        beam_x_pos = xmin - 0.1 * (xmin - xmax)
        beam_y_pos = ymin - 0.1 * (ymin - ymax)

        beam = Ellipse((beam_x_pos, beam_y_pos), width=BMAJ_pix, height=BMIN_pix,  
                   angle=BPA_deg_cartesian, edgecolor='black', facecolor='none', alpha=1, lw=2)

        ax.add_patch(beam)


        # Adjust ticks
        ax.minorticks_on()
        ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs - 10)
        ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)

        # Add a colorbar to each subplot (horizontal colorbars)
        cbar = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.04, pad=cb_pads[i])
        cbar.set_label(cb_titles[i], fontsize=cbar_fs)  # Assign corresponding label
        cbar.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")
        cbar.ax.tick_params(which='minor', length=4, direction="in")

    
    return axes  # Return axes for further modifications


# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
def create_stokes_i_plus_one_base_plot(StokesI_wcs, StokesI_stretched, 
                                       normalized_cbar_ticks, StokesI_unstretched_cbar_ticks, 
                                       other_plotting, other_plotting_label, 
                                       soft_colormap_v2, 
                                       xmin, xmax, ymin, ymax, reference_length_pix, reference_length_AU,
                                       BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
                                       max_length_pix, reference_fraction,
                                       text_fs = constants.text_fs , 
                                       axis_label_fs = constants.axis_label_fs, 
                                       axis_num_fs = constants.axis_num_fs, 
                                       cbar_fs = constants.cbar_fs):

    fig, ax = plt.subplots(1, 2, figsize=(28, 12),
                           gridspec_kw={'wspace': -0.5},  # Increase wspace for more horizontal space
                           subplot_kw={'projection': StokesI_wcs})


    im0 = ax[0].imshow(StokesI_stretched, cmap=soft_colormap_v2)
    im1 = ax[1].imshow(other_plotting, cmap=soft_colormap_v2)

    # Colorbar for Stokes I
    cbar0 = plt.colorbar(im0, ax=ax[0], orientation='horizontal', shrink=0.4)
    cbar0.set_label(r'Stokes I (mJy/beam)', fontsize=cbar_fs)
    cbar0.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")
    cbar0.ax.tick_params(which='minor', length=4, direction="in")
    cbar0.set_ticks(normalized_cbar_ticks)
    cbar0.set_ticklabels([f"{val:.2f}" for val in StokesI_unstretched_cbar_ticks])
    #cbar0.ax.set_aspect(5)  # Reduce the height of the colorbar by adjusting the aspect ratio
    
    # Colorbar for other plot
    cbar1 = plt.colorbar(im1, ax=ax[1], orientation='horizontal', shrink=0.4)
    cbar1.set_label(other_plotting_label, fontsize=cbar_fs)
    cbar1.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")
    #cbar1.ax.set_aspect(5)  # Reduce the height of the colorbar by adjusting the aspect ratio

    # Add AU scalebar
    line_x_pos = xmax - 0.05 * (xmax - xmin) 
    line_y_pos = ymax - 0.1  * (ymax - ymin) 


    # Add beam ellipse
    beam_x_pos = xmin - 0.1 * (xmin - xmax)
    beam_y_pos = ymin - 0.1 * (ymin - ymax)



    # --- Shared formatting for all axes ---
    for i in range(2):
        # Axis limits and labels
        ax[i].set_xlim(xmin, xmax)
        ax[i].set_ylim(ymin, ymax)
    
    
        ax[i].plot([(line_x_pos - reference_length_pix), (line_x_pos)], 
               [line_y_pos, line_y_pos], color='black', linewidth=3)

        ax[i].text((line_x_pos - reference_length_pix / 2), (line_y_pos - 2), 
                   f'{reference_length_AU} AU', fontsize=text_fs, ha='center', va='top')
    
        
        beam = Ellipse((beam_x_pos, beam_y_pos), width=BMAJ_pix, height=BMIN_pix,  
                   angle=BPA_deg_cartesian, edgecolor='black', facecolor='none', alpha=1, lw=2)
        
        ax[i].add_patch(beam)
        
        ax[i].set_xlabel('Right Ascension', fontsize=axis_label_fs)
        ax[i].set_ylabel('Declination', fontsize=axis_label_fs)
        
        ax[i].minorticks_on()

        ax[i].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
        ax[i].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)

    ax[1].tick_params(axis='y', labelleft=False)
    
    return fig, ax

# ---------------------------------------------------------------------------------------------





# ----------------------------------------------------------------------------------------------------------------------
def create_two_plots(StokesI_wcs,  
                     plot1_data, plot1_label, 
                     plot2_data, plot2_label, 
                     cmap, 
                     xmin, xmax, ymin, ymax, reference_length_pix, reference_length_AU,
                     BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
                     max_length_pix, reference_fraction,
                     text_fs = constants.text_fs , 
                     axis_label_fs = constants.axis_label_fs, 
                     axis_num_fs = constants.axis_num_fs, 
                     cbar_fs = constants.cbar_fs):

    fig, ax = plt.subplots(1, 2, figsize=(28, 12),
                           gridspec_kw={'wspace': -0.5},  # Increase wspace for more horizontal space
                           subplot_kw={'projection': StokesI_wcs})


    im0 = ax[0].imshow(plot1_data, cmap=cmap)
    im1 = ax[1].imshow(plot2_data, cmap=cmap)

    
    # Colorbar for other plot
    cbar1 = plt.colorbar(im0, ax=ax[0], orientation='horizontal', shrink=0.4)
    cbar1.set_label(plot1_label, fontsize=cbar_fs)
    cbar1.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")
    
    cbar2 = plt.colorbar(im1, ax=ax[1], orientation='horizontal', shrink=0.4)
    cbar2.set_label(plot2_label, fontsize=cbar_fs)
    cbar2.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")

    # Add AU scalebar
    line_x_pos = xmax - 0.05 * (xmax - xmin) 
    line_y_pos = ymax - 0.1  * (ymax - ymin) 


    # Add beam ellipse
    beam_x_pos = xmin - 0.1 * (xmin - xmax)
    beam_y_pos = ymin - 0.1 * (ymin - ymax)



    # --- Shared formatting for all axes ---
    for i in range(2):
        # Axis limits and labels
        ax[i].set_xlim(xmin, xmax)
        ax[i].set_ylim(ymin, ymax)
    
    
        ax[i].plot([(line_x_pos - reference_length_pix), (line_x_pos)], 
               [line_y_pos, line_y_pos], color='black', linewidth=3)

        ax[i].text((line_x_pos - reference_length_pix / 2), (line_y_pos - 2), 
                   f'{reference_length_AU} AU', fontsize=text_fs, ha='center', va='top')
    
        
        beam = Ellipse((beam_x_pos, beam_y_pos), width=BMAJ_pix, height=BMIN_pix,  
                   angle=BPA_deg_cartesian, edgecolor='black', facecolor='none', alpha=1, lw=2)
        
        ax[i].add_patch(beam)
        
        ax[i].set_xlabel('Right Ascension', fontsize=axis_label_fs)
        ax[i].set_ylabel('Declination', fontsize=axis_label_fs)
        
        ax[i].minorticks_on()

        ax[i].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
        ax[i].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)

    ax[1].tick_params(axis='y', labelleft=False)
    
    return fig, ax
# ----------------------------------------------------------------------------------------------------------------------








# ----------------------------------------------------------------------------------------------------------------------
def plot_grids_4x3(data_list, soft_colormap_v2, StokesI_wcs, axis_label_fs, axis_num_fs, cbar_fs, xmin, xmax, ymin, ymax, BMAJ_pix, BMIN_pix, BPA_deg_cartesian):
    grid_titles = ['Actual Data', '100% U', '50% U 50% A', '100% A']  # Updated to 4 columns
    cb_titles = [
        "Polarization Intensity", "Polarization Intensity", "Polarization Intensity", "Polarization Intensity",
        "Stokes Q", "Stokes Q", "Stokes Q", "Stokes Q",
        "Stokes U", "Stokes U", "Stokes U", "Stokes U"
    ]

    top_cb_pad = 0.1
    middle_cb_pad = 0.1
    bottom_cb_pad = 0.2

    cb_pads = [
        top_cb_pad , top_cb_pad , top_cb_pad , top_cb_pad,
        middle_cb_pad, middle_cb_pad, middle_cb_pad, middle_cb_pad,
        bottom_cb_pad, bottom_cb_pad, bottom_cb_pad, bottom_cb_pad
    ]

    POLA_cmap = soft_colormap_v2
    StokesQ_cmap = soft_colormap_v2
    StokesU_cmap = soft_colormap_v2

    cb_cmaps = [
        POLA_cmap, POLA_cmap, POLA_cmap, POLA_cmap,
        StokesQ_cmap, StokesQ_cmap, StokesQ_cmap, StokesQ_cmap,
        StokesU_cmap, StokesU_cmap, StokesU_cmap, StokesU_cmap
    ]

    # Create figure with 4x3 grid of subplots
    fig, axes = plt.subplots(3, 4, figsize=(24, 18), subplot_kw={'projection': StokesI_wcs},
                             gridspec_kw={'hspace': 0.2, 'wspace': -0.3})  # Remove spacing

    # Loop through each subplot and plot corresponding data
    for i, ax in enumerate(axes.flat):
        # Plot the data 
        im = ax.imshow(data_list[i], cmap=cb_cmaps[i])

        # Set titles for the top row only
        row, col = divmod(i, 4)  # Calculate the row and column index
        if row == 0:
            ax.set_title(grid_titles[col], fontsize=axis_label_fs)

        # Y-axis 
        # ---------------------------------------------------------------------------------------------
        if col == 0:
            ax.set_ylabel('Dec', fontsize=axis_label_fs)
            ax.tick_params(axis="y", which="both", left=True, labelleft=True)  # Show ticks and labels on left
        else:
            ax.tick_params(axis="y", which="both", left=True, labelleft=False)  # Show ticks but no labels on middle and right
            ax.set_ylabel("")  # No label on middle or right column
        # ---------------------------------------------------------------------------------------------
        
        # X-axis 
        # ---------------------------------------------------------------------------------------------
        if row == 2:
            ax.set_xlabel('RA', fontsize=axis_label_fs)
            ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)
        else:
            ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
            ax.set_xlabel("")  # No label on top or middle row
        # ---------------------------------------------------------------------------------------------

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        # Adjust ticks
        ax.minorticks_on()
        ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs - 10)
        ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)

        # Add a colorbar to each subplot (horizontal colorbars)
        cbar = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.04, pad=cb_pads[i])
        cbar.set_label(cb_titles[i], fontsize=cbar_fs)  # Assign corresponding label
        cbar.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")
        cbar.ax.tick_params(which='minor', length=4, direction="in")
        
        # Create and add the circular beam
        beam_x_pos = xmin - 0.1 * (xmin - xmax)
        beam_y_pos = ymin - 0.1 * (ymin - ymax)

        beam = Ellipse((beam_x_pos, beam_y_pos), width=BMAJ_pix, height=BMIN_pix,  
                   angle=BPA_deg_cartesian, edgecolor='black', facecolor='none', alpha=1, lw=2)

        ax.add_patch(beam)

    return axes  # Return axes for further modifications



def plot_grids_test(data_list, soft_colormap_v2, StokesI_wcs, axis_label_fs, axis_num_fs, cbar_fs, xmin, xmax, ymin, ymax, ncols, nrows, grid_titles, cb_titles, top_cb_pad, middle_cb_pad, bottom_cb_pad):
    
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6 * nrows), subplot_kw={'projection': StokesI_wcs},
                             gridspec_kw={'hspace': 0.2, 'wspace': -0.3})  # Adjust spacing
    
    cb_cmaps = [soft_colormap_v2] * (ncols * nrows)  # Assuming the same colormap for all plots
    
    # Create colorbar padding based on ncols and nrows
    cbar_pads = []
    
    for row in range(nrows):
        if row == 0:
            cbar_pads.extend([top_cb_pad] * ncols)  # Top padding for the first row
        elif row == nrows - 1:
            cbar_pads.extend([bottom_cb_pad] * ncols)  # Bottom padding for the last row
        else:
            cbar_pads.extend([middle_cb_pad] * ncols)  # Middle padding for all other rows
    
    
    for i, ax in enumerate(axes.flat):
        if i >= len(data_list):
            ax.axis("off")  # Hide empty subplots if data_list is shorter than grid size
            continue
        
        im = ax.imshow(data_list[i], cmap=cb_cmaps[i])
        row, col = divmod(i, ncols)
        
        if row == 0 and col < len(grid_titles):
            ax.set_title(grid_titles[col], fontsize=axis_label_fs)
        
        if col == 0:
            ax.set_ylabel('Dec', fontsize=axis_label_fs)
            ax.tick_params(axis="y", which="both", left=True, labelleft=True)
        else:
            ax.tick_params(axis="y", which="both", left=True, labelleft=False)
        
        if row == nrows - 1:
            ax.set_xlabel('RA', fontsize=axis_label_fs)
            ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)
        else:
            ax.tick_params(axis="x", which="both", bottom=True, labelbottom=False)
        
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.minorticks_on()
        ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs - 10)
        ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)
        
        cbar = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.04, pad=cb_pads[i])
        if i < len(cb_titles):
            cbar.set_label(cb_titles[i], fontsize=cbar_fs)
        cbar.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")
        cbar.ax.tick_params(which='minor', length=4, direction="in")
    
    return axes
# ----------------------------------------------------------------------------------------------------------------------









# ----------------------------------------------------------------------------------------------------------------------
# Adding titles to each subplot
ratio_grid_plot_titles = [
    "100% Uniform 0% Azimuthal", "0% Uniform 100% Azimuthal", 
    "90% Uniform 10% Azimuthal", "10% Uniform 90% Azimuthal",
    "80% Uniform 20% Azimuthal", "20% Uniform 80% Azimuthal",
    "70% Uniform 30% Azimuthal", "30% Uniform 70% Azimuthal",
    "60% Uniform 40% Azimuthal", "40% Uniform 60% Azimuthal",
    "50% Uniform 50% Azimuthal", "50% Uniform 50% Azimuthal"
]


def plot_ratio_grid(POLI_mJy, 
                    StokesI_wcs, soft_colormap_v2, 
                    xmin, xmax, ymin, ymax, reference_length_pix, reference_length_AU,
                    BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
                    max_length_pix, reference_fraction,
                    vector_data_actual_cartesian, vector_data_plotting_grid, 
                    fs_scale = 0.5):
    
    # Define the figure and subplots
    nrows, ncols = 6, 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 38), constrained_layout=True, 
                             subplot_kw={'projection': StokesI_wcs},
                             gridspec_kw={'wspace': -1})


    # Loop through subplots
    for i, ax in enumerate(axes.flat):
        if i >= len(ratio_grid_plot_titles):
            ax.axis("off")  # Hide empty subplots if there is extra space
            continue

        # Create blank grid
        row, col = divmod(i, ncols)
        create_blank_grid(row, col, ax, 
                          POLI_mJy, 'Polarized Intensity (mJy/beam)',  
                          StokesI_wcs, soft_colormap_v2, 
                          xmin, xmax, ymin, ymax, reference_length_pix, reference_length_AU,
                          BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
                          max_length_pix, reference_fraction, fs_scale = fs_scale)

        ax.set_title(ratio_grid_plot_titles[i], fontsize=20)

        # Add vector plots
        for row in vector_data_plotting_grid[i]:
            ax.plot([row[0], row[1]], [row[2], row[3]], color='black', lw = 3, label = 'Model')

        # Plot vector data on the first subplot (axes[0, 0])
        for row in vector_data_actual_cartesian:
            ax.plot([row[0], row[1]], [row[2], row[3]], color='red', lw = 3, label = 'Real')

        # Set axis labels and ticks
        if col == 0:
            ax.set_ylabel('Dec', fontsize=axis_label_fs)
            ax.tick_params(axis="y", which="both", left=True, labelleft=True)
        else:
            ax.tick_params(axis="y", which="both", left=True, labelleft=False)

        if row == nrows - 1:
            ax.set_xlabel('RA', fontsize=axis_label_fs)
            ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)
        else:
            ax.tick_params(axis="x", which="both", bottom=True, labelbottom=False)


        # ax.text(xmin + 5, ymax - 5, f"{i}", fontsize=16, color='blue', ha='left', va='top')


        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.minorticks_on()
        ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs - 10)
        ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)

    
    return ax
# ----------------------------------------------------------------------------------------------------------------------


# Plot CARTA slices along major and minor axis
# ----------------------------------------------------------------------------------------------------------------------
def plot_slices_along_axes(
    my_major_data, my_minor_data, 
    carta_major_data, carta_minor_data, 
    my_major_offset, my_minor_offset, 
    carta_major_offset, carta_minor_offset,
    BMAJ_deg, BMIN_deg, 
    y_label, type_of_plot="Line", 
    cb_friendly=False,
    vline = True,
    chi_sq = True,
    fit_slices = True
):
    """
    Plot data along major and minor axes.

    Parameters:
    - *_data: Your and CARTA's data for major/minor axes.
    - *_offset: Corresponding offsets for major/minor axes.
    - y_label: Label for the Y-axis.
    - type_of_plot: "Line" or "Scatter" to determine plot style.
    - cb_friendly: Whether to use a colorblind-friendly palette.
    - vline: Whether to draw a vertical line at x = 0.
    chi_sq : bool, optional
        If True, compute and print chi-squared values between user and CARTA slices.

    """
    
    # Compute chi-squared if requested
    if chi_sq:
        chi_squared_minor_values = calculate_chi_squared(my_minor_data, carta_minor_data)
        chi_squared_major_values = calculate_chi_squared(my_major_data, carta_major_data)
        
        chi_squared_minor_offset = calculate_chi_squared(my_minor_offset, carta_minor_offset)
        chi_squared_major_offset = calculate_chi_squared(my_major_offset, carta_major_offset)

        print("Chi-squared comparison between user and CARTA slices:")
        print(f"  Minor axis values χ²: {chi_squared_minor_values:.3f}")
        print(f"  Major axis values χ²: {chi_squared_major_values:.3f}")
        print(f"  Minor axis offset χ²: {chi_squared_minor_offset:.3f}")
        print(f"  Major axis offset χ²: {chi_squared_major_offset:.3f}")
        
        
    if fit_slices:
        
        # Find the area of the beam
        BMAJ_arcsec = deg_to_arcsec(BMAJ_deg)
        BMIN_arcsec = deg_to_arcsec(BMIN_deg)

        # Find the area of the beam in arcsec^2
        beam_area_arcsec2 = 1.1331 * BMAJ_arcsec * BMIN_arcsec 
    
    
        alpha_major_pos, _ , alpha_major_neg, _ = fit_slices_slope(beam_area_arcsec2, 
                                                                   my_major_data, 
                                                                   my_major_offset)
        
        alpha_minor_pos, _ , alpha_minor_neg, _ = fit_slices_slope(beam_area_arcsec2, 
                                                                   my_minor_data, 
                                                                   my_minor_offset)
        
        alpha_major_pos_carta, _ , alpha_major_neg_carta, _ = fit_slices_slope(beam_area_arcsec2, 
                                                                               carta_major_data, 
                                                                               carta_major_offset)
        
        alpha_minor_pos_carta, _ , alpha_minor_neg_carta, _ = fit_slices_slope(beam_area_arcsec2, 
                                                                               carta_minor_data, 
                                                                               carta_minor_offset)
        
        
        
        print(" ")
        print("Fitting the slope of the slices:")
        print(rf"  Major Axis α:")
        print(rf"    Mine : {alpha_major_pos:.1f} (pos) and {alpha_major_neg:.1f} (neg)")
        print(rf"    CARTA: {alpha_major_pos_carta:.1f} (pos) and {alpha_major_neg_carta:.1f} (neg)")
        print(rf"  Minor Axis α:")
        print(rf"    Mine : {alpha_minor_pos:.1f} (pos) and {alpha_minor_neg:.1f} (neg)")
        print(rf"    CARTA: {alpha_minor_pos_carta:.1f} (pos) and {alpha_minor_neg_carta:.1f} (neg)")
        
        
        
    
    if cb_friendly:
        palette = sns.color_palette("colorblind", 4)
        slices_color_my_major = palette[0]
        slices_color_carta_major = palette[1]
        slices_color_my_minor = palette[2]
        slices_color_carta_minor = palette[3]
        vline_color = 'black'
    else:
        slices_color_my_major = 'blue'
        slices_color_my_minor = 'red'
        slices_color_carta_major = 'darkorange'
        slices_color_carta_minor = 'limegreen'
        vline_color = 'black'
    

    fig, ax = plt.subplots(1, 2, figsize=(14, 6))

    my_data_lw = 5
    carta_data_lw = 5

    if type_of_plot == "Line":
        # Major axis
        ax[0].plot(my_major_offset, my_major_data, 
                   color=slices_color_my_major, ls='-', lw=my_data_lw, 
                   label="Mine")

        ax[0].plot(carta_major_offset, carta_major_data, 
                   color=slices_color_carta_major, ls='--', lw=carta_data_lw, 
                   label="CARTA")

        # Minor axis
        ax[1].plot(my_minor_offset, my_minor_data, 
                   color=slices_color_my_minor, ls='-', lw=my_data_lw, 
                   label="Mine")

        ax[1].plot(carta_minor_offset, carta_minor_data, 
                   color=slices_color_carta_minor, ls='--', lw=carta_data_lw, 
                   label="CARTA")

    elif type_of_plot == "Scatter":
        # Major axis
        ax[0].scatter(my_major_offset, my_major_data, 
                      color=slices_color_my_major, marker='o', 
                      label="Mine")

        ax[0].scatter(carta_major_offset, carta_major_data, 
                      color=slices_color_carta_major, marker='x', 
                      label="CARTA")

        # Minor axis
        ax[1].scatter(my_minor_offset, my_minor_data, 
                      color=slices_color_my_minor, marker='o', 
                      label="Mine")

        ax[1].scatter(carta_minor_offset, carta_minor_data, 
                      color=slices_color_carta_minor, marker='x', 
                      label="CARTA")

    # Add vertical lines at x = 0
    if vline:
        ax[0].axvline(0, color = vline_color, linestyle='-', linewidth=3)
        ax[1].axvline(0, color = vline_color, linestyle='-', linewidth=3)
        
        
        
    # Titles and labels
    ax[0].set_title("Major Axis Data", fontsize=title_fs)
    ax[1].set_title("Minor Axis Data", fontsize=title_fs)

    for i in range(2):
        ax[i].set_xlabel("Offset (arcsec)", fontsize=axis_label_fs)
        ax[i].set_ylabel(y_label, fontsize=axis_label_fs)
        ax[i].legend(fontsize=legend_text_fs)
        ax[i].tick_params(axis="both", which="major", direction="in", length=7, labelsize=axis_num_fs)
        
    

    plt.tight_layout()

    return ax
# ----------------------------------------------------------------------------------------------------------------------



















# ----------------------------------------------------------------------------------------------------------------------
def analyze_gaussian_averages_vs_chi2(gaussian_values, phi_values, BMAJ_pix_values, BMIN_pix_values):
    """
    Analyzes chi-squared values from Gaussian fitting results and plots their averages
    as functions of phi, BMAJ, and BMIN.

    Parameters:
    -----------
    gaussian_values : list of tuples
        Each tuple should contain (phi, BMAJ, BMIN, chi^2).
    phi_values : array-like
        Unique values of phi used in the fit grid.
    BMAJ_pix_values : array-like
        Unique values of BMAJ used in the fit grid.
    BMIN_pix_values : array-like
        Unique values of BMIN used in the fit grid.
    """

    phi_grid, BMAJ_grid, BMIN_grid, chi_grid = make_gaussian_grids(gaussian_values, phi_values, BMAJ_pix_values, BMIN_pix_values)

    # Create dictionaries to group chi^2 values by phi, BMAJ, and BMIN
    phi_dict  = defaultdict(list)
    bmaj_dict = defaultdict(list)
    bmin_dict = defaultdict(list)

    # Populate the dictionaries
    for phi, BMAJ, BMIN, chi in gaussian_values:
        phi_dict[phi].append(chi)
        bmaj_dict[BMAJ].append(chi)
        bmin_dict[BMIN].append(chi)

    # Compute average chi^2 values for each parameter value
    phi_avg  = {phi: np.mean(chi_values) for phi, chi_values in phi_dict.items()}
    bmaj_avg = {BMAJ: np.mean(chi_values) for BMAJ, chi_values in bmaj_dict.items()}
    bmin_avg = {BMIN: np.mean(chi_values) for BMIN, chi_values in bmin_dict.items()}
    # Create figure
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))

    # Plot phi vs chi^2
    ax[0].plot(list(phi_avg.keys()), list(phi_avg.values()), marker='o', linestyle='-', color = 'blue', lw = 4, ms = 15)
    ax[0].set_xlabel(r'$\phi$', fontsize = axis_label_fs)
    ax[0].set_title(r'$\phi$ vs $\chi^2$', fontsize = title_fs)

    # Plot BMAJ vs chi^2
    ax[1].plot(list(bmaj_avg.keys()), list(bmaj_avg.values()), marker='o', linestyle='-', color = 'red', lw = 4, ms = 15)
    ax[1].set_xlabel('BMAJ', fontsize = axis_label_fs)
    ax[1].set_title(r'BMAJ vs $\chi^2$', fontsize = title_fs)

    # Plot BMIN vs chi^2
    ax[2].plot(list(bmin_avg.keys()), list(bmin_avg.values()), marker='o', linestyle='-', color = 'forestgreen', lw = 4, ms = 15)
    ax[2].set_xlabel('BMIN', fontsize = axis_label_fs)
    ax[2].set_title(r'BMIN vs $\chi^2$', fontsize = title_fs)

    for i in range (3):
        # they all have the same y-axis label as chi^2
        ax[i].set_ylabel(r'$\chi^2$', fontsize = axis_label_fs)

        # Set minor and major ticks
        ax[i].minorticks_on()
        ax[i].tick_params(axis="x", which="major", direction="in", length=7, labelsize=axis_num_fs)
        ax[i].tick_params(axis="y", which="major", direction="in", length=7, labelsize=axis_num_fs)

    plt.tight_layout()
# ----------------------------------------------------------------------------------------------------------------------









# ----------------------------------------------------------------------------------------------------------------------
def plot_2d_contours_for_gaussian(gaussian_values, phi_values, BMAJ_values_pix, BMIN_values_pix, constant_parameter_name,
                                 cmap):
    
    
    phi_grid, BMAJ_grid, BMIN_grid, chi_grid = make_gaussian_grids(gaussian_values, phi_values, BMAJ_values_pix, BMIN_values_pix)
        
        
    if constant_parameter_name == 'phi':
        constant_parameter_name = '$\phi$'
        const_param_values = phi_values
        x_axis_label = 'BMAJ'
        y_axis_label = 'BMIN'
        x_axis_grid = BMAJ_grid
        y_axis_grid = BMIN_grid

    elif constant_parameter_name == 'BMAJ':
        const_param_values = BMAJ_values_pix
        x_axis_label = 'BMIN'
        y_axis_label = '$\phi$'
        x_axis_grid = BMIN_grid
        y_axis_grid = phi_grid

    elif constant_parameter_name == 'BMIN':
        const_param_values = BMIN_values_pix
        x_axis_label = 'BMAJ'
        y_axis_label = '$\phi$'
        x_axis_grid = BMAJ_grid
        y_axis_grid = phi_grid

    else:
        print("The constant_parameter options are 'phi', 'BMAJ', and 'BMIN'")
        return
    
    
    
    n_cols = 2
    n_vals = len(const_param_values)
    n_rows = math.ceil(n_vals / n_cols)
    
    # Plotting the result
    fig, ax = plt.subplots(n_rows, n_cols, figsize=(20, 18))


    # Loop over the first `num_plots` values of phi
    for i, fixed_index in enumerate(range(len(const_param_values))): 
    

        if constant_parameter_name == '$\phi$':
            x_axis_vals  = x_axis_grid[fixed_index, :, :]
            y_axis_vals  = y_axis_grid[fixed_index, :, :]
            chi          = chi_grid[fixed_index, :, :]

        elif constant_parameter_name == 'BMAJ':
            x_axis_vals  = x_axis_grid[:, fixed_index, :]
            y_axis_vals  = y_axis_grid[:, fixed_index, :]
            chi          = chi_grid[:, fixed_index, :]

        elif constant_parameter_name == 'BMIN':
            x_axis_vals  = x_axis_grid[:, :, fixed_index]
            y_axis_vals  = y_axis_grid[:, :, fixed_index]
            chi          = chi_grid[:, :, fixed_index]



#         x_axis_vals  = x_axis_grid[fixed_index, :, :]
#         y_axis_vals  = y_axis_grid[fixed_index, :, :]
#         chi          = chi_grid[fixed_index, :, :]

        row = i // n_cols  # Row index (0 or 1)
        col = i % n_cols   # Column index (0 or 1)

        # Create the contour plot
        contour = ax[row, col].contourf(x_axis_vals, y_axis_vals, chi, 50, cmap = cmap)

        # Add colorbar
        cbar = fig.colorbar(contour, ax=ax[row, col], orientation='vertical')
        cbar.set_label(fr'$\chi^2$', fontsize=cbar_fs)
        cbar.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")
        cbar.ax.tick_params(which='minor', length=4, direction="in")

        # Set title with actual c-value
        ax[row, col].set_title(fr'{constant_parameter_name} = {const_param_values[fixed_index]}', fontsize=title_fs)  

        # Labels
        ax[row, col].set_xlabel(x_axis_label, fontsize = axis_label_fs)
        ax[row, col].set_ylabel(y_axis_label, fontsize = axis_label_fs)

        # Ticks
        ax[row, col].minorticks_on()
        ax[row, col].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
        ax[row, col].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)

    plt.tight_layout()
    plt.show()
    

# ----------------------------------------------------------------------------------------------------------------------






# ----------------------------------------------------------------------------------------------------------------------
background_colours = ['#f8f9fa', '#e9ecef', '#dee2e6', '#ced4da']

def plot_3d_contours(df,
                     z_col='chi_squared',
                     n_plots=4,
                     cmap='viridis',
                     figsize=(14, 10),
                     nrows=2,
                     ncols=2):
    
    # Mapping of z_col to labels and titles
    z_col_map = {
        'chi_squared':           ('$\\chi^2$', '3D Plots of $\\chi^2$'),
        'chi_squared_norm':      ('Norm. $\\chi^2$', '3D Plots of Normalized $\\chi^2$'),
        'inv_chi_squared':       ('$1/\\chi^2$', '3D Plots of $1/\\chi^2$'),
        'inv_chi_squared_norm':  ('Norm. $1/\\chi^2$', '3D Plots of Normalized $1/\\chi^2$')
    }

    if z_col not in z_col_map:
        raise ValueError(f"Invalid z_col '{z_col}'. Must be one of: {list(z_col_map.keys())}")

    z_col_label, title = z_col_map[z_col]

    phi_values_unique = sorted(df["phi"].unique())
    fig = plt.figure(figsize=figsize)

    for i, phi_val in enumerate(phi_values_unique[:n_plots]):
        df_phi = df[df["phi"] == phi_val]
        if df_phi.empty:
            continue

        pivot = df_phi.pivot_table(index='BMIN', columns='BMAJ', values=z_col)
        if pivot.empty:
            continue

        BMAJ_vals = pivot.columns.values
        BMIN_vals = pivot.index.values
        z_vals = pivot.values

        BMAJ_grid, BMIN_grid = np.meshgrid(BMAJ_vals, BMIN_vals)

        ax = fig.add_subplot(nrows, ncols, i+1, projection='3d')
        ax.set_facecolor(background_colours[i % len(background_colours)])

        # Surface plot
        surf = ax.plot_surface(BMAJ_grid, BMIN_grid, z_vals, cmap=cmap, edgecolor='none', alpha=0.9)

        # Contour on bottom
        ax.contour(BMAJ_grid, BMIN_grid, z_vals, zdir='z',
                   offset=z_vals.min(), cmap=cmap, linestyles='solid')

        sf = 0.75  # scaling factor for font sizes

        # Labels
        ax.set_title(rf'$\phi$ = {phi_val}', fontsize=axis_label_fs)
        ax.set_xlabel('BMAJ', fontsize=axis_label_fs * sf)
        ax.set_ylabel('BMIN', fontsize=axis_label_fs * sf)
        ax.set_zlabel(z_col_label, fontsize=axis_label_fs * sf)

        ax.tick_params(axis="x", which="major", direction="in", length=7, labelsize=axis_num_fs * sf)
        ax.tick_params(axis="y", which="major", direction="in", length=7, labelsize=axis_num_fs * sf)
        ax.tick_params(axis="z", which="major", direction="in", length=7, labelsize=axis_num_fs * sf)

        # Colorbar
        mappable = plt.cm.ScalarMappable(cmap=cmap)
        mappable.set_array(z_vals)
        cbar = fig.colorbar(mappable, ax=ax, shrink=0.8, aspect=10, pad=0.1)
        cbar.set_label(z_col_label, fontsize=cbar_fs)
        cbar.ax.tick_params(labelsize=axis_num_fs * sf)

    fig.suptitle(title, fontsize=title_fs)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.90, bottom=0.05, wspace=0.3, hspace=0.2)
    plt.show()

# ----------------------------------------------------------------------------------------------------------------------
