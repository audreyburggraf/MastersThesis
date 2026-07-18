# Import the necessary packages
# ------------------------------------------------------------
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter
from matplotlib.lines import Line2D
from PlottingWithFunction import * 
from UnitConversion import * 
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
def Lin2025Figure1011(theta, mat_dist, 
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
            [(0.6, 40), (-0.14, 0.12), (0.41, 1.05)],
            [(-1.1, 1.1), (-0.05, 0.32), (-1.1, 1.1)]
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
#     plot_keys = [None, 'N12_eff', 'N22_eff', 'N33_eff', 'N34_eff', 'N44_eff'] 
    plot_keys = [None, 'N12', 'N22', 'N33', 'N34', 'N44'] 
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
        
        if key is not None:
            ax.plot(theta, mat_dist[key], color = glitterinColours[1], lw=3)

            axes[0, 0].plot(theta, mat_dist['Z11']/mat_dist['Z11'][idx_90], glitterinColours[1], lw=3)

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
def Lin2025Figure12(results, text_fs = constants.text_fs , 
                    axis_label_fs = constants.axis_label_fs, 
                    axis_num_fs = constants.axis_num_fs, 
                    fig_size_x = 15,
                    fig_size_y = 10,
                    x_num = True,
                    y_num = True,
                    full_axis_labels = True):
    
   
    # Setting up the figure and what we want to plot
    # ---------------------------------------------
    # Create a figure 
    fig, axes = plt.subplots(2, 2, figsize=(fig_size_x, fig_size_y), sharex = True)
    # ---------------------------------------------
    panel_labels = ["a)", "b)", "c)", "d)"]
    # ---------------------------------------------
    y_axis_labels = ["$\kappa_{abs, 1 mm}$ [cm$^2$ g$^{-1}]$", 
                 "$\\beta (870 \mu$m - 3.11 mm)",
                 "$\kappa_{sca, 1 mm}$ [cm$^2$ g$^{-1}]$", 
                 "$P \omega$"]
    # ---------------------------------------------
    ylims = [
        [(0.5, 3.75), None],
        [(1e-1, 90), (-0.1, 1)]
    ]

    # flatten
    flat_ylims = [item for row in ylims for item in row]

    # apply
    for ax, ylim in zip(axes.flatten(), flat_ylims):
        if ylim is not None:
            ax.set_ylim(*ylim)
    # ---------------------------------------------
    
    # Loop over each plot
    # ---------------------------------------------
    
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
        
    
    labels = {
        '870 micron': '870 $\mu$m', 
        '3.1 mm': '3.1 mm',
        '7 mm': '7 mm', 
    }

    for lambdas in ('870 micron', '3.1 mm', '7 mm'):
        
        res = results[lambdas]
        
        axes[1, 1].plot(res['rvol_max_cm'], 
                        res['P_omega'], 
                        color = glitterinColours[lambdas],
                        label = labels[lambdas],
                        lw = 3)
#     # --------------------------------------------- 
    
    
    
    # Plot the data for the left plots
    # --------------------------------------------- 
    res1 = results['1 mm']
    axes[0, 0].plot(res1['rvol_max_cm'], res1['k_abs'], color = 'black')
    axes[1, 0].plot(res1['rvol_max_cm'], res1['k_sca'], color = 'black')
    # --------------------------------------------- 
    
    
    # Customize the bottom right plot    
    # --------------------------------------------- 
    axes[1, 1].axhline(0, color = 'silver')
    axes[1, 1].legend(fontsize = 10)
    # ---------------------------------------------    
    
    # Set the x-lims of the plots
    # --------------------------------------------- 
#     axes[1, 0].set_ylim(1e-1, 90)
    # --------------------------------------------- 
    
    # Adjust the space between the plots
    fig.subplots_adjust(hspace = 0.1, wspace = 0.25)

    return fig, axes
# ----------------------------------------------------------
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
        ax.set_xlabel('Maximum grain size [$\mu$m]', fontsize=axis_label_fs)
        
        
        # Make the x-axes log scaling
        ax.set_xscale('log')
        
        # Add the wavelength text
        ax.text(0.05, 0.90,
                f"{constants.lambda_mm[idx]} mm",
                transform=ax.transAxes,
                fontsize=text_fs)
        
        # Remove scientific notation from x-axis
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
        ax.ticklabel_format(style='plain', axis='x')
        
        # Adjust ticks
        add_min_major_ticks(ax)
        add_min_major_ticks(ax2)
        
        # Plot the data 
        rvol_max_micron = cm_to_micron(results['rvol_max_cm'])
        ax_P, = ax.plot(rvol_max_micron, results['P'],  color='red',  label='P', lw = 3)
        ax2_omega, = ax2.plot(rvol_max_micron, results['omega'], color='blue', label=r'$\omega$', lw = 3)
    

    
        # Set y-lim
        ax.set_ylim(-0.2, 1.4)
        ax2.set_ylim(-0.2, 1.4)
            
        # Add a legend
        ax.legend(handles = [ax_P, ax2_omega], fontsize = legend_text_fs, loc = 'upper right')


    # Add spacing between subplots
    plt.subplots_adjust(wspace=0.4)


    
    
    return fig, axes
# ----------------------------------------------------------
#
#
#
#
#



def Pomega_Plot(results, BandsIndex,
                text_fs = constants.text_fs , 
                axis_label_fs = constants.axis_label_fs, 
                axis_num_fs = constants.axis_num_fs, 
                legend_text_fs = constants.legend_text_fs,
                fig_size_x = 8,
                fig_size_y = 5):
    
   
    # Setting up the figure and what we want to plot
    # ---------------------------------------------
    # Create a figure 
    fig, ax = plt.subplots(figsize=(fig_size_x, fig_size_y))
    # ---------------------------------------------
    
    
    # Label the x- and y-axis labels
    ax.set_ylabel('P $\omega$', fontsize=axis_label_fs)
    ax.set_xlabel('Maximum grain size [$\mu$m]', fontsize=axis_label_fs)

    
    # Make x-axis log
    ax.set_xscale('log')
    
     # Adjust ticks
    add_min_major_ticks(ax)
    

    # Loop over each wavelength
    # ---------------------------------------------
    for i, b in enumerate(BandsIndex):
        print(b)
        
        res = results[b]
        rvol_max_micron = cm_to_micron(res['rvol_max_cm'])
        
        ax.plot(rvol_max_micron, 
                res['P_omega'], 
                color = constants.alma_band_colors[b],
                label = constants.alma_band_plot_labels[b],
                lw = 3)
    # ---------------------------------------------
    
    # Add a legend
    ax.legend(fontsize = legend_text_fs)
        
    
    

    return fig, ax
# ----------------------------------------------------------
# ----------------------------------------------------------
def plot_scale_factor_resultsGlitterin(bands,
                                       bands_included_in_fit, 
                                       results, 
                                       sf_results, 
                                       df_POLF, 
                                       fig_x_size = 8,
                                       fig_y_size = 5,
                                       plot_sf = 1,
                                       custom_lw = None,
                                       chi_sq_precision = 4,
                                       band_label_fs = 15,
                                       marker_legend_size = 15,
                                      for_poster = False):
    
    sf = sf_results['best_sf']
    rvol_max_best_micron = sf_results['best_rvol_max_micron']
    
    

    # ----------------------------------------------------------------
    band_ls = [
    '-' if b in bands_included_in_fit else '--'
    for b in bands
    ]
    
    if custom_lw is None:
        # default behavior
        bands_lw = [
            3 if b in bands_included_in_fit else 2
            for b in bands
        ]
    else:
        # override where specified
        bands_lw = [
            custom_lw.get(b, 3 if b in bands_included_in_fit else 2)
            for b in bands
        ]
    # ----------------------------------------------------------------
    # Set the colors and marker size for each band 
    # ----------------------------------------------------------------
    band_colors = [alma_band_colors[b] for b in bands]
    ms = [constants.alma_band_ms[b] for b in bands]
    # ----------------------------------------------------------------
    
    
    # Get POLF values
    # ----------------------------------------------------------------
    POLF_markers = dict(zip(df_POLF["Band"], df_POLF["POLF_maxStokesI"]))
    # ----------------------------------------------------------------
    
    
    # Make the figure
    # Later I will make this so you can customize how many figures you need
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size))
    
    # Axes
    # ----------------------------------------------------------------
    # Make x-axis log
    ax.set_xscale('log')
    
    # Label the x- and y-axis labels
    ax.set_ylabel('P $\omega$', fontsize=axis_label_fs)
    ax.set_xlabel('Maximum grain size [$\mu$m]', fontsize=axis_label_fs * plot_sf)
    
    # Adjust ticks
    add_min_major_ticks(ax)
    # ----------------------------------------------------------------
    

    # Set the x-axis values
     # ----------------------------------------------------------------
    rvol_max_micron = results[bands[0]]['rvol_max_micron']
     # ----------------------------------------------------------------
    


    marker_handles = []
    
    # Loop over each band 
    for j, b in enumerate(bands):

        idx = bands.index(b)

        # Set the face and edge color of the POLF markers
        # ---------------------------------------------------
        fc = constants.alma_band_colors[b]
        ec = constants.alma_band_colors[b]
        marker_lw = 0 

        # Customize Band 5
        if b in ["Band 5"]:
            fc = 'none'
            marker_lw = 5
        # ---------------------------------------------------
        # Plot the POLF values at the best amax
        ax.scatter(rvol_max_best_micron, 
                   POLF_markers[b],
                   marker=ms[j], 
                   s=100,
                   facecolors=fc,
                   edgecolors=ec,
                   linewidths=marker_lw)
        # ---------------------------------------------------


        marker_handles.append(
        Line2D([0], [0],
               marker=ms[j],
               color='none',
               markerfacecolor=fc,
               markeredgecolor=ec,
               markeredgewidth=marker_lw,
               markersize=15,
               linestyle='None'))

        # Plot rvol_max grid and P omega multiplied by the scale factor
        # ---------------------------------------------------
        if b not in ["Band 5 robust -1", "Band 5 robust -2"]:
            print(f'[{constants.alma_band_plot_labels[b]}]')
            ax.plot(rvol_max_micron,
                    results[b]['P_omega'] * sf,
                       color=band_colors[j],
                       #label= f'{constants.alma_band_plot_labels[b]}',
                       lw=bands_lw[j],
                       ls=band_ls[j]
                      )
        # ---------------------------------------------------
    
        # Calculate chi_sq
        # ---------------------------------------------
        chi_sq_running_sum = 0

        best_idx = sf_results['best_idx_arr']

        for b in bands:

            # optional: skip Band 5 if you want
            if b == "Band 5":
                continue

            POLF_obs = POLF_markers[b]
            POLF_model = results[b]['P_omega'][best_idx] * sf

            chi_sq_running_sum += (POLF_obs - POLF_model)**2

        chi_sq = chi_sq_running_sum
        # ---------------------------------------------
        
        
        
                
    x_pos = 0.4 
    y_pos = 0.05
    y_gap = 0.1
        
        
    fs = 15
    ax.text(x_pos, y_pos + y_gap,  f'$\chi^2$ = {chi_sq:.{chi_sq_precision}f}', transform=ax.transAxes, fontsize = fs)


    if for_poster == False:
        ax.text(x_pos, y_pos + 2*y_gap, f'sf = {sf:.2f}', transform=ax.transAxes, fontsize = fs)

    # Add the a_max to the text
    ax.text(x_pos, y_pos, rf'$a_{{\mathrm{{max}}}}$ = {rvol_max_best_micron:.0f} $\mu$m', 
               transform=ax.transAxes, 
               fontsize = fs)
        
        
        
    
    
    # Add coloured text for each Band 
    # ----------------------------------------------------
    start = 0.9 # y
    c = 0
    marker_legend_handles = []
    bands_labels = ['Band 4', 'Band 5', 'Band 5 robust -1', 'Band 6', 'Band 7']
    # Loop over each band 
    for j, b_label in enumerate(bands_labels):
        
        if "robust" in b_label:
            label = b_label.replace(" robust ", "\nrob. ")
            c = c + 0.075
        else:
            label = b_label

#         ax.text(0.01, start - c,
#                    label,
#                    transform=ax.transAxes,
#                    fontsize=15,
#                    color=band_colors[j])
        
        handle = Line2D(
        [0], [0],
        marker = ms[j],  
        color='none',
        markerfacecolor=band_colors[j],
        markeredgecolor=band_colors[j],
        markersize=marker_legend_size,
        linestyle='None',
        label=label)
            
        c = c + 0.1
        marker_legend_handles.append(handle)
        # ----------------------------------------------------

    
    
    # Add a legend for the POLF markers
    # ----------------------------------------------------
    ax.legend(handles=marker_handles, 
                 fontsize = 15, 
                 frameon=False, 
#                  title="Obs. POLF", 
                 title_fontsize = 15,
                 loc = 'upper right',
                 handletextpad=-3.8,
                 handlelength=0.2)   
    # ----------------------------------------------------
    
    
    
    # Add the names of the bands in the correct colour
    # ----------------------------------------------------
    legend = ax.legend(
        handles=marker_legend_handles,
        fontsize = band_label_fs,
        frameon=False,
        loc='upper right',
        handletextpad=0.5, # Cotrols distance between marker and text 
        handlelength=1, # Cotrols distance between marker and text 
        labelspacing=0.6 # Controls vertical gap between the markers
    )
    
    for text, color in zip(legend.get_texts(), band_colors):
        text.set_color(color)
    # ----------------------------------------------------





    return fig, ax
    

# --------------------------------------------------------------------------------------