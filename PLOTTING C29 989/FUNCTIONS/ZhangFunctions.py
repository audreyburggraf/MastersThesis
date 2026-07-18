import numpy as np
import dsharp_opac as do
import matplotlib.pyplot as plt


##############################################################################################################
def get_dsharp_mix_v2(fm_ice=0.2, porosity=0.0, rule='Bruggeman'):
    """
    This method calculates the mixed mie coefficients for the DSHARP project.

    The densities and volume fractions are based on D'Alessio et al. 2001, 2006.
    However newer optical constants are used and the water fraction comes from
    Rosetta measurements.

    |                 | water ice | astrosilicates | troilite    | organics  |
    |:---------------:|:---------:|:--------------:|:-----------:|:---------:|
    | volume fraction | 0.3642    | 0.167          | 0.02578     | 0.443     |
    | solid densities | 0.92 g/cc | 3.3 g/cc       | 4.83 g/cc   | 1.5 g/cc  |

    Keywords:
    ---------

    fm_ice : float
        mass fraction of water ice.

    porosity : float
        porosity = vacuum volume fraction as float between [0, 1].
        Vaccum will be mixed in a second step with the rest using the MG-Rule.

    rule : str
        'Bruggeman' or 'Maxwell-Garnett'. Ricci et al. 2010 used 'Bruggeman'.

    Output:
    -------

    diel_constants : object of class diel_const
        the mixed dielectric constants

    rho_s : float
        the material density of the particles in g/cm**3
    """
    # define arrays for optical constants, bulk densities, and volume fractions for each species

    constants = [
        do.diel_warrenbrandt08(),
        do.diel_draine2003('astrosilicates'),
        do.diel_henning('troilite'),
        do.diel_henning('organics', refractory=True),
    ]
    # material densities

    densities = np.array([
        0.92,
        3.30,
        4.83,
        1.50])

    # fm_rest is the normalized mass fractions of the rest (adding up to 1)

    fm_rest = np.array([0.41127, 0.09292, 0.49581])
    f_mass = np.hstack((
        fm_ice,
        (1 - fm_ice) * fm_rest))

    # calculate the mean density, needed to get opacity in units of cm^2/g

    rho_s = 1.0 / (f_mass / densities).sum()

    f_vol = rho_s / densities * f_mass
    f_vol = f_vol / f_vol.sum()

    length = max([len(c.material_str) for c in constants])
    length = max([length, 16])

    print('| material'.ljust(length + 2) +
          '| volume fractions | mass fractions |')
    print('|' + (length + 1) * '-' + '|' + 18 * '-' + '|' + 16 * '-' + '|')
    for c, fv, fm in zip(constants, f_vol, f_mass):
        print('| ' + c.material_str.ljust(length) +
              '| {:.4}'.format(fv).ljust(19) + '| {:.4}'.format(fm).ljust(17) + '|')

    # mix the optical constants using the Bruggeman rule

    diel_const = do.diel_mixed(constants, f_vol, rule=rule)

    if porosity > 0:
        diel_const = do.diel_mixed([do.diel_vacuum(), diel_const], [
                                porosity, (1 - porosity)], rule='Bruggeman')
        rho_s *= 1 - porosity

    return diel_const, rho_s
##############################################################################################################




####################################################################################
def single_particle(f, a_grid_cm, lam_cm):
    
    # Make empty arrays
    Z11_90 = np.zeros(len(a_grid_cm))
    Z12_90 = np.zeros(len(a_grid_cm))
    
    # Optical constants
    oc, rho_g_cm3 = do.get_dsharp_mix(porosity = 1 - f, rule='Bruggeman')
        
    # Loop over each a value in the grid
    for j, a_max_cm in enumerate(a_grid_cm):

        # Make it an array to give to DSHARP
        a_cm = np.array([a_max_cm])

        # Compute for ALL grain sizes
        res = do.get_opacities(a_cm, lam_cm, 
                               rho_g_cm3, oc, 
                               extrapolate_large_grains=False)

        # Calculate mass
        mass_grams = (4/3) * np.pi * a_cm**3 * rho_g_cm3


        zscat = do.calculate_mueller_matrix(lam_cm, mass_grams, 
                                            res['S1'], res['S2'], 
                                            theta=res['theta'], 
                                            k_sca=res['k_sca'])['zscat']
        theta = res['theta']

        # Index closest to 90 degrees
        i90 = np.argmin(np.abs(theta - 90))

        # Average the Mueller matrix elements
        Z11_90[j] = np.sum(zscat[:, 0, i90, 0])
        Z12_90[j] = np.sum(zscat[:, 0, i90, 1])
    
        
    data = {'Z11': Z11_90, 
            'Z12': Z12_90}
    
    return data
####################################################################################
def single_particle_new(f, a_grid_cm, lam_cm):
    
    # Make empty arrays
    Z11_90 = np.zeros(len(a_grid_cm))
    Z12_90 = np.zeros(len(a_grid_cm))
    
    # Optical constants
    oc, rho_g_cm3 = get_dsharp_mix_v2(porosity = 1 - f)
        
    # Loop over each a value in the grid
    for j, a_max_cm in enumerate(a_grid_cm):

        # Make it an array to give to DSHARP
        a_cm = np.array([a_max_cm])

        # Compute for ALL grain sizes
        res = do.get_opacities(a_cm, lam_cm, 
                               rho_g_cm3, oc, 
                               extrapolate_large_grains=False)

        # Calculate mass
        mass_grams = (4/3) * np.pi * a_cm**3 * rho_g_cm3


        zscat = do.calculate_mueller_matrix(lam_cm, mass_grams, 
                                            res['S1'], res['S2'], 
                                            theta=res['theta'], 
                                            k_sca=res['k_sca'])['zscat']
        theta = res['theta']

        # Index closest to 90 degrees
        i90 = np.argmin(np.abs(theta - 90))

        # Average the Mueller matrix elements
        Z11_90[j] = np.sum(zscat[:, 0, i90, 0])
        Z12_90[j] = np.sum(zscat[:, 0, i90, 1])
    
        
    data = {'Z11': Z11_90, 
            'Z12': Z12_90}
    
    return data
####################################################################################







####################################################################################
def distribution_vs_amax(f, a_grid_cm, lam_cm):
    
    # Make empty arrays
    kabs_avg  = np.zeros(len(a_grid_cm))
    ksca_eff  = np.zeros(len(a_grid_cm))
    omega_eff = np.zeros(len(a_grid_cm))
    Z11_90    = np.zeros(len(a_grid_cm))
    Z12_90    = np.zeros(len(a_grid_cm))
    P_90      = np.zeros(len(a_grid_cm))
    
        
    # Optical constants
    oc, rho_g_cm3 = do.get_dsharp_mix(porosity=1-f)

    # Compute for ALL grain sizes
    res = do.get_opacities(a_grid_cm, lam_cm, 
                           rho_g_cm3, oc, 
                           extrapolate_large_grains=False)

    # Calculate mass
    mass_grams = (4/3) * np.pi * a_grid_cm**3 * rho_g_cm3

    # 
    zscat = do.calculate_mueller_matrix(lam_cm, 
                                        mass_grams, 
                                        res['S1'], 
                                        res['S2'], 
                                        theta=res['theta'], 
                                        k_sca=res['k_sca'])['zscat']
    k_abs = res["k_abs"]
    k_sca = res["k_sca"]
    g = res["g"]
    theta = res['theta']
    

    # Index closest to 90 degrees
    i90 = np.argmin(np.abs(theta - 90))

    # Loop over maximum grain size
    for j, amax in enumerate(a_grid_cm):

        # MRN mass weighting (q = 3.5)
        distri = a_grid_cm**0.5
        distri[a_grid_cm > amax] = 0
        distri /= distri.sum()

        # Average opacity
        kabs_avg[j] = np.sum(k_abs * distri[:, None], axis=0).item()

        ksca_avg = np.sum(k_sca * distri[:, None], axis=0).item()

        # Average asymmetry parameter
        g_avg = (np.sum(g * k_sca * distri[:, None], axis=0).item()/ ksca_avg)

        # Effective scattering opacity
        ksca_eff[j] = (1 - g_avg) * ksca_avg

        # Effective albedo
        omega_eff[j] = (
            ksca_eff[j] /
            (kabs_avg[j] + ksca_eff[j])
        )

        # Average Mueller matrix elements
        Z11_90[j] = np.sum(zscat[:, 0, i90, 0] * distri)
        Z12_90[j] = np.sum(zscat[:, 0, i90, 1] * distri)

        # Polarization fraction
        P_90[j] = np.abs(Z12_90[j] / Z11_90[j])
        
    data = {"kabs": kabs_avg,
            "ksca_eff": ksca_eff,
            "omega_eff": omega_eff,
            "Z11": Z11_90,
            "Z12": Z12_90,
            "P": P_90,
            "P_omega_eff": P_90 * omega_eff,}

    return data
####################################################################################
def distribution_vs_amax_v2(f, a_grid_cm, lam_cm):
    
    # Make empty arrays
    kabs_avg  = np.zeros(len(a_grid_cm))
    ksca_eff  = np.zeros(len(a_grid_cm))
    omega_eff = np.zeros(len(a_grid_cm))
    Z11_90    = np.zeros(len(a_grid_cm))
    Z12_90    = np.zeros(len(a_grid_cm))
    P_90      = np.zeros(len(a_grid_cm))
    
        
    # Optical constants
    oc, rho_g_cm3 = do.get_dsharp_mix(porosity=1-f)

    # Compute for ALL grain sizes
    res = do.get_opacities(a_grid_cm, lam_cm, 
                           rho_g_cm3, oc, 
                           extrapolate_large_grains=False)

    # Calculate mass
    mass_grams = (4/3) * np.pi * a_grid_cm**3 * rho_g_cm3

    # 
    zscat = do.calculate_mueller_matrix(lam_cm, 
                                        mass_grams, 
                                        res['S1'], 
                                        res['S2'], 
                                        theta=res['theta'], 
                                        k_sca=res['k_sca'])['zscat']
    k_abs = res["k_abs"]
    k_sca = res["k_sca"]
    g = res["g"]
    theta = res['theta']
    

    # Index closest to 90 degrees
    i90 = np.argmin(np.abs(theta - 90))

    # Loop over maximum grain size
    for j, amax in enumerate(a_grid_cm):

        # MRN mass weighting (q = 3.5)
        distri = a_grid_cm**0.5
        distri[a_grid_cm > amax] = 0
        distri /= distri.sum()

        # Average opacity
        kabs_avg[j] = np.sum(k_abs * distri[:, None], axis=0).item()

        ksca_avg = np.sum(k_sca * distri[:, None], axis=0).item()

        # Average asymmetry parameter
        g_avg = (np.sum(g * k_sca * distri[:, None], axis=0).item()/ ksca_avg)

        # Effective scattering opacity
        ksca_eff[j] = (1 - g_avg) * ksca_avg

        # Effective albedo
        omega_eff[j] = (
            ksca_eff[j] /
            (kabs_avg[j] + ksca_eff[j])
        )

        # Average Mueller matrix elements
        Z11_90[j] = np.sum(zscat[:, 0, i90, 0] * distri)
        Z12_90[j] = np.sum(zscat[:, 0, i90, 1] * distri)

        # Polarization fraction
        P_90[j] = np.abs(Z12_90[j] / Z11_90[j])
        
    data = {"kabs": kabs_avg,
            "ksca_eff": ksca_eff,
            "omega_eff": omega_eff,
            "Z11": Z11_90,
            "Z12": Z12_90,
            "P": P_90,
            "P_omega_eff": P_90 * omega_eff,}

    return data
####################################################################################





####################################################################################
def distribution_vs_lambda(f, a_grid_cm, lam_cm):
    

    # Make empty arrays
    kabs_avg = np.zeros(len(lam_cm))
    ksca_eff = np.zeros(len(lam_cm))
    omega_eff = np.zeros(len(lam_cm))
    Z11_90 = np.zeros(len(lam_cm))
    Z12_90 = np.zeros(len(lam_cm))
    P_90 = np.zeros(len(lam_cm))
    
        
    # Optical constants
    oc, rho_g_cm3 = do.get_dsharp_mix(porosity=1-f)

    # Compute for ALL grain sizes
    res = do.get_opacities(a_grid_cm, lam_cm, rho_g_cm3, oc,
                           extrapol=True, extrapolate_large_grains=True)

    # Calculate mass
    mass_grams = (4/3) * np.pi * a_grid_cm**3 * rho_g_cm3

    # 
    zscat = do.calculate_mueller_matrix(lam_cm, 
                                        mass_grams, 
                                        res['S1'], 
                                        res['S2'], 
                                        theta=res['theta'], 
                                        k_sca=res['k_sca'])['zscat']
    k_abs = res["k_abs"]
    k_sca = res["k_sca"]
    g = res["g"]
    theta = res['theta']
    

    # Index closest to 90 degrees
    i90 = np.argmin(np.abs(theta - 90))

    # MRN mass weighting (q = 3.5)
    distri = a_grid_cm**0.5
    distri /= distri.sum()

    # Average opacity
    kabs_avg = np.sum(k_abs * distri[:, None], axis=0)
    ksca_avg = np.sum(k_sca * distri[:, None], axis=0)

    # Average asymmetry parameter
    g_avg = np.sum(g * k_sca * distri[:, None], axis=0) / ksca_avg

    # Effective scattering opacity
    ksca_eff = (1 - g_avg) * ksca_avg

    # Effective albedo
    omega_eff = ksca_eff / (kabs_avg + ksca_eff)

    # Average Mueller matrix elements
    Z11_90 = np.sum(zscat[:,:,i90,0]*distri[:,None], axis=0)
    Z12_90 = np.sum(zscat[:,:,i90,1]*distri[:,None], axis=0)

    # Polarization fraction
    P_90 = np.abs(Z12_90 / Z11_90)

    data = {'kabs': kabs_avg,
            'ksca_eff': ksca_eff,
            'omega_eff': omega_eff,
            'Z11': Z11_90,
            'Z12': Z12_90,
            'P': P_90,
            'P_omega_eff': P_90 * omega_eff}

    return data
####################################################################################
def distribution_vs_lambda_new(f, a_grid_cm, lam_cm):
    

    # Make empty arrays
    kabs_avg = np.zeros(len(lam_cm))
    ksca_eff = np.zeros(len(lam_cm))
    omega_eff = np.zeros(len(lam_cm))
    Z11_90 = np.zeros(len(lam_cm))
    Z12_90 = np.zeros(len(lam_cm))
    P_90 = np.zeros(len(lam_cm))
    
        
    # Optical constants
    oc, rho_g_cm3 = get_dsharp_mix_v2(porosity=1-f, rule='Bruggeman')

    # Compute for ALL grain sizes
    res = do.get_opacities(a_grid_cm, lam_cm, rho_g_cm3, oc,
                           extrapol=True, extrapolate_large_grains=True)

    # Calculate mass
    mass_grams = (4/3) * np.pi * a_grid_cm**3 * rho_g_cm3

    # 
    zscat = do.calculate_mueller_matrix(lam_cm, 
                                        mass_grams, 
                                        res['S1'], 
                                        res['S2'], 
                                        theta=res['theta'], 
                                        k_sca=res['k_sca'])['zscat']
    k_abs = res["k_abs"]
    k_sca = res["k_sca"]
    g = res["g"]
    theta = res['theta']
    

    # Index closest to 90 degrees
    i90 = np.argmin(np.abs(theta - 90))

    # MRN mass weighting (q = 3.5)
    distri = a_grid_cm**0.5
    distri /= distri.sum()

    # Average opacity
    kabs_avg = np.sum(k_abs * distri[:, None], axis=0)
    ksca_avg = np.sum(k_sca * distri[:, None], axis=0)

    # Average asymmetry parameter
    g_avg = np.sum(g * k_sca * distri[:, None], axis=0) / ksca_avg

    # Effective scattering opacity
    ksca_eff = (1 - g_avg) * ksca_avg

    # Effective albedo
    omega_eff = ksca_eff / (kabs_avg + ksca_eff)

    # Average Mueller matrix elements
    Z11_90 = np.sum(zscat[:,:,i90,0]*distri[:,None], axis=0)
    Z12_90 = np.sum(zscat[:,:,i90,1]*distri[:,None], axis=0)

    # Polarization fraction
    P_90 = np.abs(Z12_90 / Z11_90)

    data = {'kabs': kabs_avg,
            'ksca_eff': ksca_eff,
            'omega_eff': omega_eff,
            'Z11': Z11_90,
            'Z12': Z12_90,
            'P': P_90,
            'P_omega_eff': P_90 * omega_eff}

    return data
####################################################################################






####################################################################################


####################################################################################
def plot_fig3(data, f_values, bands_cm):

    # Make the figure size 
    fig, axes = plt.subplots(2, 4, figsize=(15, 5), sharex=True)


    # Label y-axes
    axes[0, 0].set_ylabel('$\log_{10}(\kappa_{abs})$', fontsize = 20)
    axes[1, 0].set_ylabel('$P \omega_{eff}, P, \omega_{eff}$', fontsize = 20)

    # axes[0, 0].legend(fontsize = 15, frameon = False)


    # title_list = ['160 $\mu$m', '1 mm', '1 cm', '10 cm']


    colours = ['red', 'blue', 'black', 'magenta']

    # Loop over the columns and rows 
    for col in range(4):

        f = f_values[col]
        

        for row in range(2):

            ax = axes[row, col]
            ax.set_xlim(-3, 1)


            ax.tick_params(which='major', direction='in', 
                           top=True, right=True, length=5, width=1.2, labelsize=16)


            # Top row
            if row == 0:
                ax.set_title(rf"f = {f}", fontsize = 20)

                # Loop over bands
                for i, lam in enumerate(bands_cm):
                    ax.plot(np.log10(data[lam][f]['a_grid_cm'] * f), 
                            np.log10(data[lam][f]['data']['kabs']), color = colours[i])


            # Bottom row
            if row == 1:
                tol = 0.1
                ax.set_ylim(0 - tol, 1 + tol)
                ax.set_xlabel("$\log_{10}(a_{max}f/cm)$", fontsize = 20)

                # Loop over bands
                for i, lam in enumerate(bands_cm):
                    log_af_grid = np.log10(data[lam][f]['a_grid_cm'] * f)

                    ax.plot(log_af_grid, data[lam][f]['data']['omega_eff'], color = colours[i], ls = ':')
                    ax.plot(log_af_grid, data[lam][f]['data']['P'], color = colours[i], ls = 'dashdot')
                    ax.plot(log_af_grid, data[lam][f]['data']['P_omega_eff'], color = colours[i], ls = '-')

            if col in (1, 2, 3):
                ax.tick_params(labelleft=False)
                
                
         



    axes[0, 0].legend(fontsize = 15, frameon = False)

    axes[-1, 0].set_xticks([-3, -2, -1, 0, 1])
    axes[-1, 0].set_xticklabels([-3, -2, -1, 0, 1], fontsize=16)

    for ax in axes[0, :]:
        ax.set_yticks([-2, -1, 0])
        ax.set_yticklabels([-2, -1, 0], fontsize=16)

    plt.tight_layout()
    #plt.show()
    
    return fig, ax
###################################################################################
def plot_fig3_with_files(data, f_values_idx, wavelengths):

    # Make the figure size 
    fig, axes = plt.subplots(3, 4, figsize=(15, 8), sharex=True)


    # Label y-axes
    axes[0, 0].set_ylabel('$\log_{10}(\kappa_{abs})$', fontsize = 20)
    axes[1, 0].set_ylabel('$Z_{11}, Z_{12}$', fontsize = 20)
    axes[2, 0].set_ylabel('$P \omega_{eff}, P, \omega_{eff}$', fontsize = 20)

    # axes[0, 0].legend(fontsize = 15, frameon = False)


    # title_list = ['160 $\mu$m', '1 mm', '1 cm', '10 cm']


    colours = ['red', 'blue', 'black', 'magenta']

    # Loop over the columns and rows 
    for col, f_idx in enumerate(f_values_idx):
        
        f = float(f_idx)
        print(f'f = {f} and f_idx = {f_idx}')

        
        

        for row in range(3):

            ax = axes[row, col]
            ax.set_xlim(-3, 1)


            ax.tick_params(which='major', direction='in', 
                           top=True, right=True, length=5, width=1.2, labelsize=16)


            # Top row
            if row == 0:
                ax.set_ylim(-2.5, 0.5)
                ax.set_title(rf"f = {f}", fontsize = 20)

                # Loop over bands
                for i, lam in enumerate(wavelengths):
                    ax.plot(np.log10(data[lam][f_idx]['a'] * f), 
                            np.log10(data[lam][f_idx]['k_abs']), color = colours[i])
                    
            # Middle row
            if row == 1:
#                 tol = 0.1
#                 ax.set_ylim(0 - tol, 1 + tol)
#                 ax.set_xlabel("$\log_{10}(a_{max}f/cm)$", fontsize = 20)

                # Loop over bands
                for i, lam in enumerate(wavelengths):
                    log_af_grid = np.log10(data[lam][f_idx]['a'] * f)

                    ax.plot(log_af_grid, data[lam][f_idx]['Z11eff'], color = colours[i], ls = ':')
                    ax.plot(log_af_grid, data[lam][f_idx]['Z12eff'], color = colours[i], ls = '-')


            # Bottom row
            if row == 2:
                tol = 0.1
                ax.set_ylim(0 - tol, 1 + tol)
                ax.set_xlabel("$\log_{10}(a_{max}f/cm)$", fontsize = 20)

                # Loop over bands
                for i, lam in enumerate(wavelengths):
                    log_af_grid = np.log10(data[lam][f_idx]['a'] * f)

                    ax.plot(log_af_grid, data[lam][f_idx]['omega_eff'], color = colours[i], ls = ':')
                    ax.plot(log_af_grid, data[lam][f_idx]['P'], color = colours[i], ls = 'dashdot')
                    ax.plot(log_af_grid, data[lam][f_idx]['P_omega_eff'], color = colours[i], ls = '-')
                    
                    
            

            if col in (1, 2, 3):
                ax.tick_params(labelleft=False)
                
                
         



    axes[0, 0].legend(fontsize = 15, frameon = False)

    axes[-1, 0].set_xticks([-3, -2, -1, 0, 1])
    axes[-1, 0].set_xticklabels([-3, -2, -1, 0, 1], fontsize=16)

    for ax in axes[0, :]:
        ax.set_yticks([-2, -1, 0])
        ax.set_yticklabels([-2, -1, 0], fontsize=16)

    plt.tight_layout()
    #plt.show()
    
    return fig, ax
###################################################################################
def plot_fig_15(single_data, distribution_data, f_values, a_grid_cm):
    # Make the plot
    fig, ax = plt.subplots(2, 2, figsize=(15, 9), sharex = True)

    lw = 2
    for i, (a, f) in enumerate(zip(ax.ravel(), f_values)):

        data1 = single_data[f]
        data2 = distribution_data[f]

        a.plot(a_grid_cm, data1['Z11'], ls = '-', lw = lw, color = 'blue', 
               label = '$Z_{11}$' if i == 0 else None)
        a.plot(a_grid_cm, data1['Z12'], ls = '-', lw = lw, color = 'green', 
               label = '$Z_{12}$' if i == 0 else None)

        a.plot(a_grid_cm, data2['Z11'], ls = '--', lw = lw, color = 'blue')
        a.plot(a_grid_cm, data2['Z12'], ls = '--', lw = lw, color = 'green')

        a.set_xscale('log')


        a.set_xlim(1e-5,10)
    #     a.set_ylim(1e-3,2e3)

        a.text(0.05, 0.95, f'f = {f_values[i]}', transform=a.transAxes, 
               fontsize=20, va='top', ha='left')

        a.tick_params(
            which='major',
            direction='in',
            top=True,
            right=True,
            length=5,
            width=1.2,
            labelsize=16)




        # Bottom row 
        for a in ax[-1, :]:
            a.set_xlabel(r'$a\ [{\rm cm}]$', fontsize=18)

        # Left column
        for a in ax[:, 0]:
            a.set_ylabel(r'$Z_{11}$ or $Z_{12}$ at $90^\circ$', fontsize=18)



    # Top left 
    # ---------------------------------------------------------------------------------------------------------
    ax = ax[0, 0]
    ax.text(0.05, 0.85, f'$\lambda$ = 1 mm', transform=ax.transAxes, 
            fontsize=20, va='top', ha='left')
    leg1 = ax.legend(loc = 'lower left', frameon = False, fontsize = 15)  

    # Dummy handles for line styles
    solid, = ax.plot([], [], 'k-', lw=2, label='Single size')
    dashed, = ax.plot([], [], 'k--', lw=2, label=r'$q=3.5$')

    leg2 = ax.legend(handles=[solid, dashed],
                     loc='lower left',
                     bbox_to_anchor=(0.25, 0),   # adjust this
                     frameon=False, fontsize = 15)

    # Add the first legend back
    ax.add_artist(leg1)

    sigma = 1.5
    ax.set_ylim(-5 - sigma, 5 + sigma)
    # ---------------------------------------------------------------------------------------------------------

    plt.subplots_adjust(wspace=0.2, hspace=0.2)


    plt.show()
###################################################################################