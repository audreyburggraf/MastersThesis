# Import the necessary packages
# ------------------------------------------------------------
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter
from matplotlib.lines import Line2D
from PlottingWithFunction import * 
# ------------------------------------------------------------
#
#
#
#
## Import the constants
# ------------------------------------------------------------
import sys

sys.path.append("/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/") 

import constants
# ------------------------------------------------------------
#
#
#
#
#
# ------------------------------------------------------------
# Custom colours from the Daniel Lin 2025 paper
# ------------------------------------------------------------
glitterinColours = {
    # Green
    1: "#65c2a5",
    '870 micron': "#65c2a5",
    'teal': "#65c2a5",
    # Orange
    29: "#fc8d63",
    '3.1 mm': "#fc8d63",
    'orange' : "#fc8d63",
    # Blue
    90: "#8da0cb",
    '7 mm': "#8da0cb",
    'blue': "#8da0cb",

}
# ------------------------------------------------------------
#
#
#
#
#
# ----------------------------------------------------------
# Function to replicate Figure 3
# ----------------------------------------------------------
def Lin2025Figure3(x_vol, x_enc, omega, 
                text_fs = constants.text_fs , 
                axis_label_fs = constants.axis_label_fs, 
                axis_num_fs = constants.axis_num_fs, 
                cbar_pad = 0.1,
                cbar_shrink = 1,
                fig_size_x = 5,
                fig_size_y = 5,
                x_num = True,
                y_num = True,
                full_axis_labels = True):
   
    
    # Create a figure with the WCS projection
    fig, ax = plt.subplots(figsize=(fig_size_x, fig_size_y))
   
    
    # Label the axes
    ax.set_xlabel(rf'$x_{{vol}}$', fontsize=axis_label_fs)
    ax.set_ylabel(rf'$\omega$',    fontsize=axis_label_fs)
        
   
    # Plot the data 
    ax.plot(x_vol, omega, color = 'mediumaquamarine', lw = 5, alpha = 0.75)

    # Make axes log-space
    ax.set_yscale('log')
    ax.set_xscale('log')


    # Adjust ticks
    add_min_major_ticks(ax)
#     ax.minorticks_on()
#     ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
#     ax.tick_params(axis="y", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)



    return fig, ax
# ------------------------------------------------------------
#
#
#
#
#
# ----------------------------------------------------------
# Function to replicate Figure 4 and 5
# ----------------------------------------------------------
def Lin2025Figure45(x_vol, 
                    result, 
                    text_fs = constants.text_fs , 
                    axis_label_fs = constants.axis_label_fs, 
                    axis_num_fs = constants.axis_num_fs, 
                    cbar_pad = 0.1,
                    cbar_shrink = 1,
                    fig_size_x = 17,
                    fig_size_y = 10,
                    x_num = True,
                    y_num = True,
                    full_axis_labels = True,
                    x_axis_labels = 'bottom', 
                    sharex = True):
   
    plotting_labels = ['Q_ext_vol',
                       'Q_abs_vol',
                       'albedo',
                       'Q_ext_pja',
                       'Q_abs_pja',
                       'epsilon']
    
    # Create a figure with the WCS projection
    fig, axes = plt.subplots(2, 3, figsize=(fig_size_x, fig_size_y), sharex = sharex)
   
    
    panel_labels = ["a)", "b)", "c)", "d)", "e)", "f)"]
    y_axis_labels = ["$Q_{ext}^{vol}$", 
                     "$Q_{abs}^{vol}$", 
                     "$\omega$", 
                     "$Q_{ext}^{pja}$", 
                     "$Q_{abs}^{pja}$", 
                     "$\epsilon$"]
    
    for i, ax in enumerate(axes.flatten()):
        
        ax.text(0.05, 0.85, f"{panel_labels[i]}",
                       transform=ax.transAxes,
                       fontsize=text_fs)
        
        if x_axis_labels == 'bottom':
            # Add x-labels to bottom row:
            if i in (3, 4, 5):
                ax.set_xlabel(rf'$x_{{vol}}$', fontsize=axis_label_fs)
        
        elif x_axis_labels == 'both':
              ax.set_xlabel(rf'$x_{{vol}}$', fontsize=axis_label_fs)
            
        # Label the y-axes
        ax.set_ylabel(rf'{y_axis_labels[i]}',    fontsize=axis_label_fs)
        
        # PLot data
        label = plotting_labels[i]

        if label in result:
            ax.plot(x_vol, result[label], color=glitterinColours[1], lw=3)

        # Make axes log-space
        ax.set_yscale('log')
        ax.set_xscale('log')


        # Adjust ticks
        add_min_major_ticks(ax)
        
        # Add the horizontal lines that he has on his plots
        # ------------------------------------------------------------
        ax_top = ax.twiny()
        ax_top.set_xscale('log')

        ax_top.plot(result['x_enc'], result[label], alpha=0)  # invisible, just to sync limits
        ax_top.set_xlim(result['x_enc'].min(), result['x_enc'].max())
        
        
        if i in (0, 1, 2):  # top row only
            ax_top.set_xlabel(r'$x_{\mathrm{enc}}$', fontsize=axis_label_fs)
            add_min_major_ticks(ax_top)
        else:
            ax_top.set_xticklabels([])
    # ------------------------------------------------------------

    # Add the horizontal lines that he has on his plots
    # ------------------------------------------------------------
    # Top right
    axes[0, 0].axhline(1, color = 'darkgrey',  ls = '-', lw = 3)
    axes[0, 0].axhline(2, color = 'black', ls = ':', lw = 3)
    # ------------------------------------------------------------
    # Bottom left
    axes[1, 0].axhline(2, color = 'darkgrey',  ls = '-', lw = 3)
    # ------------------------------------------------------------
    # Bottom right
    axes[1, 2].axhline(1e0, color = 'black', ls = ':', lw = 3)
    # ------------------------------------------------------------
    
    
    # Adjust the space between the plots
    fig.subplots_adjust(hspace = 0.3, wspace = 0.35)

    return fig, axes
# ----------------------------------------------------------
#
#
#
#
#
# ----------------------------------------------------------
# Function to replicate Figure 8
# ----------------------------------------------------------
def Lin2025Figure8(text_fs = constants.text_fs , 
                axis_label_fs = constants.axis_label_fs, 
                axis_num_fs = constants.axis_num_fs, 
                cbar_pad = 0.1,
                cbar_shrink = 1,
                fig_size_x = 15,
                fig_size_y = 10,
                x_num = True,
                y_num = True,
                full_axis_labels = True):
   
    
    # Create a figure with the WCS projection
    fig, axes = plt.subplots(2, 3, figsize=(fig_size_x, fig_size_y), sharex = True)
   
    
    panel_labels = ["a)", "b)", "c)", "d)", "e)", "f)"]
    y_axis_labels = ["$2 \pi Z_11 / C_{sca}$", "$-Z_{12}/Z_{11}$", "$Z_{22}/Z_{11}$", "$Z_{33}/Z_{11}$", "$Z_{34}/Z_{11}$", "$Z_{44}/Z_{11}$"]
    
    for i, ax in enumerate(axes.flatten()):
        
        ax.text(0.05, 0.85, f"{panel_labels[i]}",
                       transform=ax.transAxes,
                       fontsize=text_fs)
        
        # Add x-labels to bottom row:
        if i in (3, 4, 5):
            ax.set_xlabel(rf'$\theta$ [$^\circ$]', fontsize=axis_label_fs)
        
        # Label the y-axes
        ax.set_ylabel(rf'{y_axis_labels[i]}',    fontsize=axis_label_fs)


#     # Make axes log-space
#     ax.set_yscale('log')
#     ax.set_xscale('log')


#     # Adjust ticks
#        add_min_major_ticks(ax)
#     ax.minorticks_on()
#     ax.tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)
#     ax.tick_params(axis="y", which="major", direction="in", bottom=True, top=True, length=7, labelsize=axis_num_fs)



    return fig, ax
# ----------------------------------------------------------
#
#
#
#
#
# ----------------------------------------------------------
# Function to replicate Figure 10 and 11
# ----------------------------------------------------------
def Lin2025Figure1011(theta, results, 
                      idx = -1, 
                      text_sf = 1, 
                    text_fs = None, 
                    axis_label_fs = None, 
                    axis_num_fs = None, 
                    fig_size_x = 17,
                    fig_size_y = 10,
                    x_num = True,
                    y_num = True,
                    full_axis_labels = True,
                    x_axis_labels = 'bottom', 
                    sharex = True,
                    fig_style = 10,
                     hspace = 0.3, 
                      wspace = 0.5,
                     hlines = True):
    
    # Assign defaults inside the function
    if text_fs is None:
        text_fs = constants.text_fs * text_sf
    if axis_label_fs is None:
        axis_label_fs = constants.axis_label_fs * text_sf
    if axis_num_fs is None:
        axis_num_fs = constants.axis_num_fs * text_sf
        
    
    # Create a figure 
    fig, axes = plt.subplots(2, 3, figsize=(fig_size_x, fig_size_y), sharex = sharex)
    
    
    ylims_by_style = {
        10: [
            [(0.4, 1e3), (-0.2, 0.2), (0.3, 1.05)],
            [(-1.1, 1.1), (-0.7, 0.15), (-1.1, 1.05)]
        ],
        11: [
            [None, (-0.14, 0.22), (0.41, 1.05)],
            [(-1.1, 1.1), (-0.05, 0.35), (-1.1, 1.1)]
        ]
    }
    
    ylims = ylims_by_style.get(fig_style)

    if ylims is not None:
        flat_ylims = [item for row in ylims for item in row]

        for ax, ylim in zip(axes.flatten(), flat_ylims):
            if ylim is not None:
                ax.set_ylim(*ylim)
    
    axhlines = {
        10: [None, 0.2, 0.4, -0.75, -0.1, -0.25],
        11: [None, 0.1, 0.5, -0.5, 0.15, -0.25]
    }
    
       
    idx_90 = np.argmin(np.abs(theta - 90))
    print(rf'Sanity check: theta[idx_90] = {theta[idx_90]}')

   
    # Arrays for the loop
    # ------------------------------------------------------
    panel_labels = ["a)", "b)", "c)", "d)", "e)", "f)"]
    # ------------------------------------------------------
    y_axis_labels = ["normalized $Z_{11}$", 
                     "$- Z_{12}/Z_{11}$", 
                     "$Z_{22}/Z_{11}$", 
                     "$Z_{33}/Z_{11}$", 
                     "$Z_{34}/Z_{11}$", 
                     "$Z_{44}/Z_{11}$"]
    


    # ------------------------------------------------------
    plot_keys = [None, 'N12_eff', 'N22_eff', 'N33_eff', 'N34_eff', 'N44_eff'] 
    # plot_keys = [None, 'N12', 'N22', 'N33', 'N34', 'N44'] 
    # plot_keys = [None, 'N12_fig10', 'N22_fig10', 'N33_fig10', 'N34_fig10', 'N44_fig10'] 
    
    # ------------------------------------------------------
    
    for i, ax in enumerate(axes.flatten()):
        
        lines = axhlines.get(fig_style)

        if hlines:
            if lines[i] is not None:
                ax.axhline(lines[i], color='lightsteelblue', lw=3, ls='--')
        
        ax.text(0.05, 0.85, f"{panel_labels[i]}",
                       transform=ax.transAxes,
                       fontsize=text_fs)
        
        if x_axis_labels == 'bottom':
            # Add x-labels to bottom row:
            if i in (3, 4, 5):
                ax.set_xlabel(rf'$\theta$ [deg]', fontsize=axis_label_fs)
        
        elif x_axis_labels == 'both':
              ax.set_xlabel(rf'$\theta$ [deg]', fontsize=axis_label_fs)
            
        # Label the y-axes
        ax.set_ylabel(rf'{y_axis_labels[i]}',    fontsize=axis_label_fs)


        # Adjust ticks
        add_min_major_ticks(ax, sf = text_sf)
        
        # Add lines
        ax.axvline(90, color = 'silver', lw = 2, ls = '-')
        
        if i in (1, 3, 4, 5):
            ax.axhline(0, color = 'silver', lw = 2, ls = '-')
            
            
        ax.set_xlim(0, 180)
        
            
    
        # Plot the data 
        # ------------------------------------------------------------
        key = plot_keys[i]
        
        colors = (glitterinColours[1], 'red', 'orange', 'green', 'indigo', 'pink')
        
        if key is not None:
            for j in range(len(results)):
                ax.plot(theta, results[j][key][idx], color = colors[j], lw=3, label = rf'{j + 1}')
                
                axes[0, 0].plot(theta, results[j]['Z11'][:, idx]/results[j]['Z11'][idx_90, idx], color = colors[j], lw=3)
                
        # ------------------------------------------------------
            
        
        
    
        
    # Make axes log-space
    # ------------------------------------------------------------
    axes[0, 0].set_yscale('log')
    # ------------------------------------------------------------
    
    


    
    # Adjust the space between the plots
    fig.subplots_adjust(hspace = hspace, wspace = wspace)

    return fig, axes
# ----------------------------------------------------------
#
#
#
#
#
# ----------------------------------------------------------
# Function to replicate Figure 12
# ----------------------------------------------------------
def Lin2025Figure12(results_all,
                    text_fs = constants.text_fs , 
                    axis_label_fs = constants.axis_label_fs, 
                    axis_num_fs = constants.axis_num_fs, 
                    fig_size_x = 15,
                    fig_size_y = 10,
                    x_num = True,
                    y_num = True,
                    full_axis_labels = True):
   
    
    # Create a figure 
    fig, axes = plt.subplots(2, 2, figsize=(fig_size_x, fig_size_y), sharex = True)
   
    
    panel_labels = ["a)", "b)", "c)", "d)"]
    
    
    y_axis_labels = ["$\kappa_{abs, 1 mm}$ [cm$^2$ g$^{-1}]$", 
                     "$\\beta (870 \mu$m - 3.11 mm)",
                     "$\kappa_{sca, 1 mm}$ [cm$^2$ g$^{-1}]$", 
                     "$P \omega$"]
    
    for i, ax in enumerate(axes.flatten()):
        
        ax.text(0.05, 0.85, f"{panel_labels[i]}",
                       transform=ax.transAxes,
                       fontsize=text_fs)
        
        # Add x-labels to bottom row:
        if i in (2, 3):
            ax.set_xlabel(rf'Maximum $r_{{vol}}$ [cm]', fontsize=axis_label_fs)
        
        # Label the y-axes
        ax.set_ylabel(rf'{y_axis_labels[i]}',    fontsize=axis_label_fs)


        # Make axes log-space
        if i in (0, 2):
            ax.set_yscale('log')
        ax.set_xscale('log')

        # Adjust ticks
        add_min_major_ticks(ax)
        
        # Set the x-lims of the plots
        ax.set_xlim(1e-3, 1e0)
        
    # Plot the data 
    # --------------------------------------------- 
    axes[0, 0].plot(results_all['1 mm']['r_vol'], results_all['1 mm']['k_abs_eff'])
#     axes[0, 1].plot(results_all['1 mm']['r_vol'], results_all['1 mm']['Q_ext_vol'])
    axes[1, 0].plot(results_all['1 mm']['r_vol'], results_all['1 mm']['k_sca_eff'])
    
    labels = {
        '870 micron': '870 $\mu$m', 
        '3.1 mm': '3.1 mm',
        '7 mm': '7 mm', 
    }

    for lambdas in ('870 micron', '3.1 mm', '7 mm'):
        axes[1, 1].plot(results_all[lambdas]['r_vol'], 
                        results_all[lambdas]['P_omega_90'], 
                        color = glitterinColours[lambdas],
                        label = labels[lambdas],
                        lw = 3)
    # --------------------------------------------- 
    
    # Customize the bottom right plot    
    # --------------------------------------------- 
    axes[1, 1].axhline(0, color = 'silver')
    axes[1, 1].legend(fontsize = 20)
    # ---------------------------------------------    
    
    # Set the x-lims of the plots
    # --------------------------------------------- 
    axes[1, 0].set_ylim(1e-1, 90)
    # --------------------------------------------- 
    
    # Adjust the space between the plots
    fig.subplots_adjust(hspace = 0.1, wspace = 0.25)

    return fig, axes
# ----------------------------------------------------------
#
#
#
#
#
# ----------------------------------------------------------
# Function to look at P and omega twin axes
# ----------------------------------------------------------
def P_vs_omega_twin_axis(results_all, wavelength_labels, wavelengths_cm):
    
    # Make the figure 
    fig, axes = plt.subplots(2, 2, figsize=(20, 15)) 

    
    # Loop over all axes
    for i, ax in enumerate(axes.flatten()):
        
        # Set how you index the wavelength 
        idx = wavelength_labels[i]
        results = results_all[idx]
        
        # Make the twin axis
        ax2 = ax.twinx()
        
        # Label the x- and y-axis labels
        ax.set_ylabel('P', fontsize=axis_label_fs)
        ax2.set_ylabel(rf'$\omega$', fontsize=axis_label_fs)
        ax.set_xlabel(rf'Maximum $r_{{vol}}$ [cm]', fontsize=axis_label_fs)
        
        
        # Make the x-axes log scaling
        ax.set_xscale('log')
        
        # Add the wavelength text
        ax.text(0.05, 0.90,
                f"{wavelengths_cm[idx]} cm",
                transform=ax.transAxes,
                fontsize=text_fs)
        
        # Adjust ticks
        add_min_major_ticks(ax)
        add_min_major_ticks(ax2)
        
        # Plot the data 
        ax_P, = ax.plot(results['r_vol'], results['P'],  color='red',  label='P', lw = 3)
        ax2_omega, = ax2.plot(results['r_vol'], results['omega'], color='blue', label=r'$\omega$', lw = 3)
    

        # Add a legend
        ax.legend(handles = [ax_P, ax2_omega], fontsize = legend_text_fs, loc = 'upper right')


    # Add spacing between subplots
    plt.subplots_adjust(wspace=0.4)


    
    
    return fig, axes