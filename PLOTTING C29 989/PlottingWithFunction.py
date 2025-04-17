# Import necessary packages
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# Import Functions
from FITS_Image_Functions import *
from PolarizationFunctions import *


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

def create_base_plot(StokesI_wcs, StokesI_stretched, soft_colormap_v2, 
                     normalized_cbar_ticks, StokesI_unstretched_cbar_ticks, 
                     xmin, xmax, ymin, ymax, reference_length_pix, reference_length_AU,
                     text_fs, axis_label_fs, axis_num_fs, cbar_fs,
                     BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
                     max_length_pix, reference_fraction):
    
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
    ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
    ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize=axis_num_fs)

    return fig, ax


# ---------------------------------------------------------------------------------------------


def create_blank_grid(i, j, fs_scale, ax, StokesI_wcs, StokesI_stretched, soft_colormap_v2, 
                      normalized_cbar_ticks, StokesI_unstretched_cbar_ticks, 
                      xmin, xmax, ymin, ymax, reference_length_pix, reference_length_AU,
                      text_fs, axis_label_fs, axis_num_fs, cbar_fs,
                      BMAJ_pix, BMIN_pix, BPA_deg_cartesian, 
                      max_length_pix, reference_fraction):
    
    
    # Plot the data on the given ax
    im = ax.imshow(StokesI_stretched, cmap=soft_colormap_v2)

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
def plot_grids(data_list, soft_colormap_v2, StokesI_wcs, axis_label_fs, axis_num_fs, cbar_fs, xmin, xmax, ymin, ymax):
    grid_titles = ['100% U', '50% U 50% A', '100% A']
    cb_titles = [
        "Polarization Angle", "Polarization Angle", "Polarization Angle",
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










def plot_grids_4x3(data_list, soft_colormap_v2, StokesI_wcs, axis_label_fs, axis_num_fs, cbar_fs, xmin, xmax, ymin, ymax):
    grid_titles = ['Actual Data', '100% U', '50% U 50% A', '100% A']  # Updated to 4 columns
    cb_titles = [
        "Polarization Angle", "Polarization Angle", "Polarization Angle", "Polarization Angle",
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

