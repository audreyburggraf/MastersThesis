import numpy as np
import dsharp_opac as do
import matplotlib.pyplot as plt


# Import Functions
from PlottingWithFunction import *

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




def do_things(file):
    
    d = file
    
    a      = d['a']
    lam    = d['lam']
    k_abs  = d['k_abs']
    k_sca  = d['k_sca']
    gsca   = d['g']
    rho_s  = d['rho_s']
    m   = 4 * np.pi / 3 * rho_s * a**3
    S1     = d['S1']
    S2     = d['S2']
    theta =  d['theta']
    g     =  d['g']
    
    
    lam_calc = np.array([lam])

    S1_calc = S1[:, None, :]      # (258,1,5)
    S2_calc = S2[:, None, :]      # (258,1,5)
    k_sca_calc = k_sca[:, None]   # (258,1)

    # index where theta = 90 degree
    i90 = theta.searchsorted(90)
    zscat = do.calculate_mueller_matrix(
                                        lam_calc,
                                        m,
                                        S1_calc,
                                        S2_calc,
                                        theta=theta,
                                        k_sca=k_sca_calc,
                                    )["zscat"]
    
#     print('in do_things')
#     print(zscat.shape)
#     print(zscat[100, 0, i90, 0])   # Z11
#     print(zscat[100, 0, i90, 1])   # Z12


    # make a loop over the maximum particle size
    k_a = np.zeros_like(k_abs)
    k_s = np.zeros_like(k_sca)
    k_s_eff = np.zeros_like(k_sca)
    g_s = np.zeros_like(gsca)

    Z11eff = np.zeros([len(a)])
    Z12eff = np.zeros_like(Z11eff)
    omega = np.zeros_like(Z11eff)
    omega_eff = np.zeros_like(Z11eff)

    for ia, _a in enumerate(a):

        # make a normalized MRN distribution

        distri = a**0.5
        distri[a > _a] = 0
        distri /= distri.sum()

        # average the opacities to get the average albedo

        k_a[ia] = np.sum(k_abs * distri)
        k_s[ia] = np.sum(k_sca * distri)

        #k_s_eff = k_s*(1-g_s)
        # The way caculate g here is wrong!!!        
        #g_s[ia] = np.sum(gsca * distri[:, None], 0)


        k_s_eff[ia] = np.sum(k_sca * (1 - g) * distri)

                # average the scattering matrix elements

        Z11eff[ia] = np.sum(zscat[:, :, i90, 0] * distri[:, None], 0)
        Z12eff[ia] = np.sum(zscat[:, :, i90, 1] * distri[:, None], 0)


    # calculate polarization fraction at 90 degree

    P = np.abs(Z12eff / Z11eff)

    #chi = k_a + k_s_eff
    #chi_iso = k_a + k_s
    omega_eff = (k_s_eff)/(k_s_eff + k_a)
    omega = (k_s_eff)/(k_s + k_a)
        
    return {
        'a': a,
        'lam': lam,
        'k_abs': k_abs,
        'k_sca': k_sca,
        'g': g,
        'rho_s': rho_s,
        'S1': S1,
        'S2': S2,
        'theta': theta,
        'P': P,
        'omega': omega,
        'omega_eff': omega_eff,
        'k_s_eff': k_s_eff,
        'Z11eff': Z11eff,
        'Z12eff': Z12eff,
        'P_omega': P * omega,
        'P_omega_eff': P * omega_eff,
    }
        
# ---------------------------------------------------------------------   
    
    
    
    
    
    
    
    
# ---------------------------------------------------------------------   
# ---------------------------------------------------------------------   
def Pw_base_plot(data, f_values, wavelengths, 
                 colours,
                 y_axis_plot_P = True,
                 y_axis_plot_omega = True,
                 y_axis_plot_Pomega = True,
                 x_axis = 'a', # The other option is 'af'
                 ymin = -0.1,  ymax = 1.7, 
                 xmin = -0.1,  xmax = 1.7, 
                 plot_sf = 1.5, # This controls how much the text labels are scaled by
                 fig_x = 30, fig_y = 8,
                 text_x_pos = 0.05,
                 lw = 3):
    
    
    # Make figures
    fig, axes = plt.subplots(1, len(f_values), figsize = (fig_x, fig_y))
    
    # Loop over the f values you want to plot
    #########################################
    for i, f in enumerate(f_values): 
        
        ax = axes[i]
        
        # y-axis
        # -----------------------------------------------------------------
        # Label the y-axis on far left 
        if i == 0:
            ax.set_ylabel('P $\omega$', fontsize=axis_label_fs * plot_sf)
        else:
            ax.set_yticklabels([])
            ax.set_ylabel('')
            ax.tick_params(left=False)
        # -----------------------------------------------------------------
        # Set y min and max
        ax.set_ylim(ymin, ymax)
        # -----------------------------------------------------------------
        
        
        # Loop over wavelengths and plot data
        # -----------------------------------------------------------------
        for j, wav in enumerate(wavelengths):
            a = data[wav][str(f)]['a']
            
            if x_axis == 'log af':
                x_data = np.log10(a*f)
                ax.set_xlabel(r'$\log_{10}(a_{max}f/cm)$', fontsize=axis_label_fs*plot_sf)
                ax.set_xlim(-3, 1)
                ax.set_xticks([-3, -2, -1, 0, 1])
                ax.set_xticklabels([-3, -2, -1, 0, 1], fontsize=16)
            
            elif x_axis == 'a':
                x_data = a
                ax.set_xlabel(r'$a_{\mathrm{max}}$ $[\mu m]$', fontsize=axis_label_fs*plot_sf)
                ax.set_xscale("log")
                
                
            if y_axis_plot_P:
                ax.plot(x_data, data[wav][str(f)]['P'], ls = 'dashdot', lw = lw, color = colours[j])
            
            if y_axis_plot_omega:
                ax.plot(x_data, data[wav][str(f)]['omega_eff'], ls = ':', lw = lw,  color = colours[j])
            
            if y_axis_plot_Pomega: 
                ax.plot(x_data, data[wav][str(f)]['P_omega_eff'], ls = '-', lw = lw,  color = colours[j])
                
            
        # -----------------------------------------------------------------
            
            
        # Add ticks
        add_min_major_ticks(ax, sf = plot_sf)
        
        # Add f label in top left
        ax.text(text_x_pos, 0.9, f'$f$ = {f:.2f}', transform=ax.transAxes, fontsize = 30)
            
    #########################################
    
    return fig, axes
# ---------------------------------------------------------------------   
# --------------------------------------------------------------------- 




# ---------------------------------------------------------------------   
# ---------------------------------------------------------------------   
# --------------------------------------------------------------------------------------
# def Pw_sf_plot_with_files(bands,                  # These are the bands we are looking at
#                           bands_labels, 
#                           bands_included_in_fit, # These are the values that were fit
#                           f_values,            # These are the porosity values we are testing
#                           a_max_dist_micron_by_f,
#                           P_times_omega_by_f,
#                           sf_by_f,             # This is the values the data is being scaled by 
#                           best_idx_by_f, 
#                           best_a_max_by_f,     # These are the best a_max values for each f
#                           POLF_markers,        # These are the values for each band we are plotting the marker at
#                           chi_sq_by_f, 
#                           plot_sf = 1.5,       # This controls how much the text labels are scaled by
#                           ymin = -0.1, 
#                           ymax = 1.7, 
#                           custom_lw = None,
#                           custom_text_x = None,
#                           chi_sq_precision = 3,
#                           fig_x = 30, fig_y = 8,
#                           for_poster = False,
#                           f_fs = 45, 
#                           f_x = 0.05, f_y = 0.9,
#                           spine_width = 1,
#                           included_lw = 5,
#                           marker_border = False,
#                           BandLegend_bbox_to_anchor=(1.0, 1.0),
#                           BandLegend_handletextpad=0.5,
#                           BandLegend_handlelength=1,
#                           BandLegend_labelspacing=0.6
#                           ):
#     band_ls = [
#     '-' if b in bands_included_in_fit else '--'
#     for b in bands
#     ]
    
#     # Set default line widths
#     # ------------------------------------------------

#     excluded_lw = included_lw * (2/3)

#     if custom_lw is None:
#         bands_lw = [
#             included_lw if b in bands_included_in_fit else excluded_lw
#             for b in bands
#         ]
#     else:
#         bands_lw = [
#             custom_lw.get(
#                 b,
#                 included_lw if b in bands_included_in_fit else excluded_lw
#             )
#             for b in bands
#         ]
#     # ------------------------------------------------
        
#     # Change sizes of things if for_poster or not
#     # ---------------------------
#     if for_poster:
#         num_fs = constants.poster_axis_num_fs
#         xy_axis_fs = constants.poster_axis_label_fs 
#         stats_fs = constants.poster_text2_fs
#         legend_marker_size = constants.poster_colour_bands_legend_ms
#         band_legend_fs = constants.poster_colour_bands_fs
#         marker_size = constants.poster_colour_bands_ms
#     # ---------------------------
#     else: 
#         num_fs = axis_num_fs*plot_sf
#         xy_axis_fs = axis_label_fs * plot_sf
#         stats_fs = 30
#         legend_marker_size = 20
#         band_legend_fs = 30
#         marker_size = 400
#     # ---------------------------
    
#     # Set the colors and marker size for each band 
#     band_colors = []
#     ms = []

#     for b in bands:
#         if for_poster and b == "Band 5 robust -1":
#             band_colors.append(alma_band_colors["Band 5"])
#             ms.append(constants.alma_band_ms["Band 5"])
#         else:
#             band_colors.append(alma_band_colors[b])
#             ms.append(constants.alma_band_ms[b])
    
    
    
#     # Make the figure
#     # Later I will make this so you can customize how many figures you need
#     fig, ax = plt.subplots(1, 4, figsize=(fig_x, fig_y))
    
    
    
    
#     # Loop over the f values you want to plot
#     for i, f in enumerate(f_values): 
        
#         for spine in ax[i].spines.values():
#             spine.set_linewidth(spine_width)
        
       

#         if i == 0:
#             ax[i].set_ylabel('P $\omega$', fontsize = xy_axis_fs)
#         else:
#             ax[i].set_yticklabels([])
#             ax[i].set_ylabel('')
#             ax[i].tick_params(left=False)

#         # Make the x-axis log scale 
#         ax[i].set_xscale("log")


            
#         if for_poster:
#             ax[i].set_xlabel('Maximum grain size [$\mu$m]', fontsize = xy_axis_fs)
        
#         else:
#             ax[i].set_xlabel(r'$a_{\mathrm{max}}  [\mu\mathrm{m}]$', fontsize = xy_axis_fs)
            
    
#         # Set ticks
#         ax[i].minorticks_on()
#         ax[i].tick_params(axis="x", which="minor", direction="in", left=True, right=True, length=4, width=spine_width)
#         ax[i].tick_params(axis="y", which="minor", direction="in", left=True, right=True, length=4, width=spine_width)
        
        
#         # Control the font size of the numbers on the axes
#         # ---------------------------------------------------------------------------------------
            
#         ax[i].tick_params(axis="x", which="major", direction="in", bottom=True, top=True, length=7, labelsize = num_fs, width=spine_width)
#         ax[i].tick_params(axis="y", which="major", direction="in", left=True, right=True, length=7, labelsize = num_fs, width=spine_width)
#         # ---------------------------------------------------------------------------------------
        
        
#         # ---------------------------------------------------------------------------------------
#         # ---------------------------------------------------------------------------------------

#         # f labels on plot
#         # ---------------------------------------------------------------------------------------
#         if for_poster == False:
#             ax[i].set_title(f'$f$ = {f:.2f}', fontsize=50)
            
#         else:
#             ax[i].text(f_x, f_y,  f'$f$ = {f:.2f}', transform=ax[i].transAxes, fontsize = constants.poster_text1_fs)
            
        
# #         if for_poster == True:
# #             ax[i].text(x_pos, 0.15,  f'$\chi^2$ = {chi_sq:.{chi_sq_precision}f}', transform=ax[i].transAxes, fontsize = 30)
#         # ---------------------------------------------------------------------------------------

#         ax[i].set_ylim(ymin, ymax)
#         # ------------------------------------------
        
        
#         marker_handles = []
#         # plot each band for this f
#         for j, b in enumerate(bands):
            
#             idx = bands.index(b)
            
            
# #             print(rf'At band {b}, the POLF_obs value is: {POLF_obs[f][idx]}')
# #             print(' ')

#             fc = band_colors[j]
    
#             if marker_border:
#                 ec = 'black'
#                 marker_lw = 1.5
#             else:
#                 ec = band_colors[j]
#                 marker_lw = 0
#             #ec = band_colors[j]
# #             marker_lw = 0 

#             if for_poster and b == "Band 5 robust -1":
#                 fc = alma_band_colors["Band 5"]
#                 ec = alma_band_colors["Band 5"]

#             if b == "Band 5" and not for_poster:
#                 fc = 'none'
#                 marker_lw = 5

#             # Plot the POLF markers
#             plot_marker = True

#             if for_poster and b == "Band 5":
#                 plot_marker = False

#             if plot_marker:
#                 ax[i].scatter(
#                     best_a_max_by_f[f],
#                     POLF_markers[b],
#                     marker=ms[j], 
#                     s=marker_size,
#                     facecolors=fc,
#                     edgecolors=ec,
#                     linewidths=marker_lw
#                 )
#             # ---------------------------------------------------------------------

        
#             marker_handles.append(
#             Line2D([0], [0],
#                    marker=ms[j],
#                    color='none',
#                    markerfacecolor=fc,
#                    markeredgecolor=ec,
#                    markeredgewidth=marker_lw,
#                    markersize=15,
#                    linestyle='None'))

#             if b not in ["Band 5 robust -1", "Band 5 robust -2"]:

#                 ax[i].plot(a_max_dist_micron_by_f[f],
#                                P_times_omega_by_f[f][:, idx] * sf_by_f[f],
#                                color=band_colors[j],
#                                label=b,
#                                lw=bands_lw[j],
#                                ls=band_ls[j])
        
                
#         x_pos = 0.05 
        
#         if custom_text_x is not None:
#             x_pos = x_pos + custom_text_x[i]
        
        
        
#         ax[i].text(x_pos, 0.15,  f'$\chi^2$ = {chi_sq_by_f[f]:.{chi_sq_precision}f}', transform=ax[i].transAxes, fontsize = stats_fs)

#         # Only add sf if for_poster = False
#         if for_poster == False:
#             ax[i].text(x_pos, 0.25, f'sf = {sf_by_f[f]:.2f}', transform=ax[i].transAxes, fontsize = stats_fs)


#         # Add the a_max to the text
#         ax[i].text(x_pos, 0.05, rf'$a_{{\mathrm{{max}}}}$ = {best_a_max_by_f[f]:.0f} $\mu$m', 
#                    transform=ax[i].transAxes, 
#                    fontsize = stats_fs)
        
        
        
#     #plt.tight_layout()
    
#     plt.subplots_adjust(
#         left=0.05,
#         right=0.98,
#         bottom=0.15,
#         top=0.95,
#         wspace=0.02
#     )
    
    
#     # Add coloured text for each Band in the far right plot
#     # ----------------------------------------------------
#     # ----------------------------------------------------
#     #start = start # This is the starting vertical position
#     c = 0
#     marker_legend_handles = []
    
#     # Choose legend entries
#     if for_poster:
#         legend_bands = ["Band 4", "Band 5 robust -1", "Band 6", "Band 7"]
#     else:
#         legend_bands = bands_labels

#     # Loop over legend entries
#     for b_label in legend_bands:

#         if for_poster and b_label == "Band 5 robust -1":
#             label = "Band 5"
#             color = alma_band_colors["Band 5"]
#             marker = constants.alma_band_ms["Band 5"]

#         else:
#             label = b_label
#             color = alma_band_colors[b_label]
#             marker = constants.alma_band_ms[b_label]

#         handle = Line2D(
#             [0], [0],
#             marker=marker,
#             color='none',
#             markerfacecolor=color,
#             markeredgecolor=color,
#             markersize=legend_marker_size,
#             linestyle='None',
#             label=label
#         )

#         marker_legend_handles.append(handle)
        
        
        
#     legend_elements = [
#         Line2D([0], [0], color='black', linestyle='-', lw = 5, label='Included'),
#         Line2D([0], [0], color='black', linestyle='--', lw = 5, label='Excluded')
#     ]
#     # ----------------------------------------------------
#     # ----------------------------------------------------
        
#     legend = ax[3].legend(
#             handles=marker_legend_handles,
#             fontsize= band_legend_fs,
#             frameon=False,
#             loc='upper right',
#             bbox_to_anchor=BandLegend_bbox_to_anchor,
#             handletextpad=BandLegend_handletextpad, # Cotrols distance between marker and text 
#             handlelength=BandLegend_handlelength, # Cotrols distance between marker and text 
#             labelspacing=BandLegend_labelspacing # Controls vertical gap between the markers
#         )
    
#     if for_poster:
#         legend_colors = [
#             alma_band_colors["Band 4"],
#             alma_band_colors["Band 5"],
#             alma_band_colors["Band 6"],
#             alma_band_colors["Band 7"]
#         ]
#     else:
#         legend_colors = band_colors

#     for text, color in zip(legend.get_texts(), legend_colors):
#         text.set_color(color)
# #     # ----------------------------------------------------
# #     # ----------------------------------------------------



# # 

#     return fig, ax


# ---------------------------------------------------------------------   
# ---------------------------------------------------------------------   