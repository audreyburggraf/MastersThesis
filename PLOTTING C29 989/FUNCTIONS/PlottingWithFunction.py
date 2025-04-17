# Import necessary packages
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns

# Import Functions
from FITS_Image_Functions import *
from PolarizationFunctions import *

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


# # Define the general 'read_in' function
# def read_in(group_dict):
#     """
#     Takes a dictionary of key-value pairs and unpacks them into variables.
#     """
#     # Unpack the dictionary keys and values dynamically into local variables
#     for key, value in group_dict.items():
#         globals()[key] = value


# # Make function
# def plot_fits_data(plotting_data, what_i_am_plotting, units, StokesI_header, StokesI_wcs, StokesI_vmin_custom, StokesI_vmax_custom, data_files, min_str, max_str, distance_pc, reference_length_AU, reference_fraction, max_length_pix, cmap, font_sizes):
    
#     # Now you can call read_in for either font_sizes or constants
#     read_in(font_sizes)  # This will create variables like title_fs, axis_label_fs, etc.
#     # read_in(constants)  # This will create variables like distance_pc, reference_length_AU, etc.

    
    
#     # Maximum and minimum of the plot
#     # ------------------------------------------------------------------------------------
#     RA_min_pix, Dec_min_pix = string_to_pixel(min_str, StokesI_wcs)
#     RA_max_pix, Dec_max_pix = string_to_pixel(max_str, StokesI_wcs)

#     xmin = RA_max_pix
#     xmax = RA_min_pix

#     ymin = Dec_min_pix
#     ymax = Dec_max_pix
#     # ------------------------------------------------------------------------------------


#     # Create a figure with the WCS projection
#     # ------------------------------------------------------------------------------------
#     fig, ax = plt.subplots(figsize=(14, 12), subplot_kw={'projection': StokesI_wcs})
#     # ------------------------------------------------------------------------------------

#     # Add data
#     # ------------------------------------------------------------------------------------
#     im = ax.imshow(plotting_data, 
#                    cmap = cmap) 
#     # ------------------------------------------------------------------------------------


#     # Colorbar
#     # ------------------------------------------------------------------------------------
#     cbar = plt.colorbar(im, ax=ax)

#     cbar.set_label(f'{what_i_am_plotting} ({units})', 
#                    fontsize = font_sizes['cbar_fs'])

#     cbar.ax.tick_params(labelsize= font_sizes['axis_num_fs'] , which='major', length=7, direction="in")
#     cbar.ax.tick_params(which='minor', length=4, direction="in")

#     # Update colorbar ticks and labels
#     if what_i_am_plotting == 'Stokes I':
        
#         normalized_cbar_ticks = np.array([0, 0.2, 0.4, 0.6, 0.8, 1])

       
#         StokesI_unstretched_cbar_ticks = unstretch(normalized_cbar_ticks, StokesI_vmin_custom, StokesI_vmax_custom)
        
       
        
#         cbar.set_ticks(normalized_cbar_ticks)  # Set tick locations in stretched space
        
#         cbar.set_ticklabels([f"{val:.2f}" for val in StokesI_unstretched_cbar_ticks])
# #     # ------------------------------------------------------------------------------------


#     # Add title and axis labels
#     # ------------------------------------------------------------------------------------
#     # ax.set_title('Title', fontsize=title_fs, fontweight='bold')
#     ax.set_xlabel('Right Ascension', fontsize = font_sizes['axis_label_fs'])
#     ax.set_ylabel('Declination', fontsize = font_sizes['axis_label_fs'])
#     # ------------------------------------------------------------------------------------


#     # Set x and y limits using ax.set_xlim() and ax.set_ylim()
#     # ------------------------------------------------------------------------------------------
#     ax.set_xlim(xmin, xmax)
#     ax.set_ylim(ymin, ymax)
#     # ------------------------------------------------------------------------------------------


#     # Add line and text for reference_length AU 
#     # ------------------------------------------------------------------------------------
#     reference_length_pix = length_in_pixels(reference_length_AU, distance_pc, StokesI_header)
#     # ------------------------------------------------------------------------------------------
#     # Plot the line in axes coordinates
#     line_x_pos = xmax - 0.05 * (xmax - xmin) # 5% in from the right
#     line_y_pos = ymax - 0.1  * (ymax - ymin) # 10% down from the top

#     ax.plot([(line_x_pos - reference_length_pix), (line_x_pos)], 
#             [line_y_pos, line_y_pos],
#             color='black', 
#             linewidth=3)

#     # Add the text label centered below the line
#     ax.text((line_x_pos - reference_length_pix/2), 
#             (line_y_pos - 2), 
#             f'{reference_length_AU} AU', 
#             fontsize = font_sizes['text_fs'], 
#             ha='center', 
#             va='top') 
#     # ------------------------------------------------------------------------------------


#     # Beam
#     # ----------------------------------------------------------------------------------------
#     # Get the beam information
#     beam_info = get_beam_info(StokesI_header)

#     BMAJ_deg, BMIN_deg, BMAJ_pix, BMIN_pix, BPA_astronomy_deg, BPA_deg_cartesian = beam_info
#     # ----------------------------------------------------------------------------------------
#     # Get the position of the beam
#     beam_x_pos = xmin - 0.1 * (xmin - xmax) # 10% in from the left
#     beam_y_pos = ymin - 0.1 * (ymin - ymax) # 10% up from the bottom
#     # ----------------------------------------------------------------------------------------
#     # Add the beam to the plot
#     beam = Ellipse(
#         (beam_x_pos, beam_y_pos),                      
#         width = BMAJ_pix, 
#         height = BMIN_pix,   
#         angle = BPA_deg_cartesian,
#         edgecolor='black',                 # Edge color
#         facecolor='none',                  # Fill color
#         alpha=1,
#         lw = 2)                           # Transparency
    
#     ax.add_patch(beam)
#     # ------------------------------------------------------------------------------------



#     # Vectors
#     # ------------------------------------------------------------------------------------
#     # Set parameters
#     step = 7  # Plot vectors every step-th pixel
#     nx, ny = data_files['StokesI_data_2d_mJy'].shape
#     # ------------------------------------------------------------------------------------
#     # Loop over values in x and y
#     for x in range(0, nx, step):
#         for y in range(0, ny, step):
#             # Check if the conditions are met
#             if (data_files['StokesI_data_2d_mJy'][y, x] / data_files['StokesIerr_data_2d_mJy'][y, x] > 3 and 
#                 data_files['calculated_polarized_intensity'][y, x] / data_files['PolarizedIntensity_err_data_2d_mJy'][y, x] > 3 and 
#                 data_files['PolarizationAngle_err_data_2d_deg'][y, x] < 10): 


#                 # Get the polarization angle at this pixel
#                 angle_rad = data_files['calculated_polarization_angle_rad'][y, x] + np.pi / 2

#                 # Get the polarization fraction at this pixel
#                 polarization_fraction = data_files['calc_polarized_frac'][y, x]  # Add your data source

#                 # Scale vector length by polarization fraction
#                 vector_length_pix = max_length_pix * polarization_fraction

#                 # Compute the vector components
#                 dx = vector_length_pix * np.cos(angle_rad)
#                 dy = vector_length_pix * np.sin(angle_rad)

#                 # Plot the vector
#                 ax.plot([x - dx / 2, x + dx / 2], [y - dy / 2, y + dy / 2], color='black')
#     # ------------------------------------------------------------------------------------



#     # Draw the reference vector line
#     # ------------------------------------------------------------------------------------
#     vector_x_pos = xmax - 0.05 * (xmax - xmin)  # 5% in from the right
#     vector_y_pos = ymin - 0.1 * (ymin - ymax)  # 10% up from the bottom

#     ax.plot([(vector_x_pos - max_length_pix * reference_fraction), vector_x_pos], 
#             [vector_y_pos, vector_y_pos],
#             color='black', 
#             linewidth=3)

#     # Add the text label centered below the line
#     ax.text(vector_x_pos - max_length_pix * reference_fraction / 2,  # Midpoint of the line
#             vector_y_pos - 2,  # Adjusted position below the line
#             f'{reference_fraction * 100:.0f}%',  # Format the fraction as a percentage
#             fontsize= font_sizes['text_fs'], 
#             ha='center', 
#             va='top')
#     # ------------------------------------------------------------------------------------


#     # Adjust ticks and numbers for bottom and left axes
#     # ------------------------------------------------------------------------------------
#     ax.minorticks_on()

#     ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize= font_sizes['axis_num_fs'])
    
#     ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=font_sizes['axis_num_fs'])
#     # ------------------------------------------------------------------------------------


#     # Tight layout and show the plot
#     # ------------------------------------------------------------------------------------------
#     # plt.tight_layout()
#     plt.show()
#     # ------------------------------------------------------------------------------------------

#     return 




# ---------------------------------------------------------------------------------------------

# def create_base_plot(StokesI_wcs, StokesI_stretched, soft_colormap_v2, 
#                      normalized_cbar_ticks, StokesI_unstretched_cbar_ticks, 
#                      xmin, xmax, ymin, ymax, reference_length_pix, reference_length_AU,
#                      text_fs, axis_label_fs, axis_num_fs, cbar_fs,
#                      BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
#                      max_length_pix, reference_fraction):
    
#     # Create a figure with the WCS projection
#     fig, ax = plt.subplots(figsize=(14, 12), subplot_kw={'projection': StokesI_wcs})

#     # Add data
#     im = ax.imshow(StokesI_stretched, cmap=soft_colormap_v2)

#     # Colorbar
#     cbar = plt.colorbar(im, ax=ax)
#     cbar.set_label(r'Stokes I (mJy/beam)', fontsize=cbar_fs)
#     cbar.ax.tick_params(labelsize=axis_num_fs, which='major', length=7, direction="in")
#     cbar.ax.tick_params(which='minor', length=4, direction="in")
#     cbar.set_ticks(normalized_cbar_ticks)
#     cbar.set_ticklabels([f"{val:.2f}" for val in StokesI_unstretched_cbar_ticks])

#     # Add axis labels
#     ax.set_xlabel('Right Ascension', fontsize=axis_label_fs)
#     ax.set_ylabel('Declination', fontsize=axis_label_fs)

#     # Set x and y limits
#     ax.set_xlim(xmin, xmax)
#     ax.set_ylim(ymin, ymax)

#     # Add line and text for 100 AU 
#     line_x_pos = xmax - 0.05 * (xmax - xmin) 
#     line_y_pos = ymax - 0.1  * (ymax - ymin) 

#     ax.plot([(line_x_pos - reference_length_pix), (line_x_pos)], 
#             [line_y_pos, line_y_pos],
#             color='black', linewidth=3)

#     ax.text((line_x_pos - reference_length_pix/2), (line_y_pos - 2), 
#             f'{reference_length_AU} AU', fontsize=text_fs, ha='center', va='top') 

#     # Create and add the circular beam
#     beam_x_pos = xmin - 0.1 * (xmin - xmax)
#     beam_y_pos = ymin - 0.1 * (ymin - ymax)

#     beam = Ellipse((beam_x_pos, beam_y_pos), width=BMAJ_pix, height=BMIN_pix,  
#                    angle=BPA_deg_cartesian, edgecolor='black', facecolor='none', alpha=1, lw=2)

#     ax.add_patch(beam)

# #     # Draw the reference vector line
# #     vector_x_pos = xmax - 0.05 * (xmax - xmin)
# #     vector_y_pos = ymin - 0.1 * (ymin - ymax)

# #     ax.plot([(vector_x_pos - max_length_pix * reference_fraction), vector_x_pos], 
# #             [vector_y_pos, vector_y_pos], color='black', linewidth=3)

# #     ax.text(vector_x_pos - max_length_pix * reference_fraction / 2,  
# #             vector_y_pos - 2,  
# #             f'{reference_fraction * 100:.0f}%',  
# #             fontsize=text_fs, ha='center', va='top')

#     # Adjust ticks
#     ax.minorticks_on()
#     ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
#     ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)

#     return fig, ax







# ---------------------------------------------------------------------------------------------------------------
def create_stokes_i_base_plot(StokesI_wcs, StokesI_stretched, soft_colormap_v2, 
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
    im = ax.imshow(StokesI_stretched, cmap=soft_colormap_v2)

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


def create_base_plot(StokesI_wcs, plotting_data, cbar_label, soft_colormap_v2, 
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
    im = ax.imshow(plotting_data, cmap=soft_colormap_v2)

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


def create_blank_grid(i, j, fs_scale, ax, data_plotting, StokesI_wcs, StokesI_stretched, soft_colormap_v2, 
                      normalized_cbar_ticks, StokesI_unstretched_cbar_ticks, 
                      xmin, xmax, ymin, ymax, reference_length_pix, reference_length_AU,
                      text_fs, axis_label_fs, axis_num_fs, cbar_fs,
                      BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
                      max_length_pix, reference_fraction):
    
    
    # Plot the data on the given ax
    im = ax.imshow(data_plotting, cmap=soft_colormap_v2)

    # Colorbar for every rightmost column (in a 6x2 grid) to avoid redundancy
    if j == 1:
        cbar = plt.colorbar(im, ax=ax, pad=0.02, fraction=0.05)
        cbar.set_label(r'Stokes I (mJy/beam)', fontsize=cbar_fs* fs_scale)
        cbar.ax.tick_params(labelsize=axis_num_fs* fs_scale, which='major', length=7, direction="in")
        cbar.ax.tick_params(which='minor', length=4, direction="in")
        cbar.set_ticks(normalized_cbar_ticks)
        cbar.set_ticklabels([f"{val:.2f}" for val in StokesI_unstretched_cbar_ticks])

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






# Plot CARTA slices along major and minor axis
# ----------------------------------------------------------------------------------------------------------------------
def plot_slices_along_axes(
    my_major_data, my_minor_data, 
    carta_major_data, carta_minor_data, 
    my_major_offset, my_minor_offset, 
    carta_major_offset, carta_minor_offset,
    y_label, type_of_plot="Line", 
    cb_friendly=False,
    vline = True
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
    """
    
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
                   label="My Major Axis Data")

        ax[0].plot(carta_major_offset, carta_major_data, 
                   color=slices_color_carta_major, ls='--', lw=carta_data_lw, 
                   label="CARTA Major Axis Data")

        # Minor axis
        ax[1].plot(my_minor_offset, my_minor_data, 
                   color=slices_color_my_minor, ls='-', lw=my_data_lw, 
                   label="My Minor Axis Data")

        ax[1].plot(carta_minor_offset, carta_minor_data, 
                   color=slices_color_carta_minor, ls='--', lw=carta_data_lw, 
                   label="CARTA Minor Axis Data")

    elif type_of_plot == "Scatter":
        # Major axis
        ax[0].scatter(my_major_offset, my_major_data, 
                      color=slices_color_my_major, marker='o', 
                      label="My Major Axis Data")

        ax[0].scatter(carta_major_offset, carta_major_data, 
                      color=slices_color_carta_major, marker='x', 
                      label="CARTA Major Axis Data")

        # Minor axis
        ax[1].scatter(my_minor_offset, my_minor_data, 
                      color=slices_color_my_minor, marker='o', 
                      label="My Minor Axis Data")

        ax[1].scatter(carta_minor_offset, carta_minor_data, 
                      color=slices_color_carta_minor, marker='x', 
                      label="CARTA Minor Axis Data")

    # Add vertical lines at x = 0
    if vline:
        ax[0].axvline(0, color = vline_color, linestyle='-', linewidth=1)
        ax[1].axvline(0, color = vline_color, linestyle='-', linewidth=1)
        
        
        
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

