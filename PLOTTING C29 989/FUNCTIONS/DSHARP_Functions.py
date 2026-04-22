import numpy as np 
import miepython
from scipy.interpolate import interp1d
import dsharp_opac as do


import pandas as pd



def compute_P_vs_angle(scatt_angles_deg, lambda_bands_cm, a_max_list, a_grid, zscat, theta, res):
    Z11eff = np.zeros((len(a_max_list), len(lambda_bands_cm), len(scatt_angles_deg)))
    Z12eff = np.zeros_like(Z11eff)

    distri_template = a_grid**0.5
    distri_template /= distri_template.sum()

    for ia, amax in enumerate(a_max_list):
        distri = distri_template.copy()
        distri[a_grid > amax] = 0
        distri /= distri.sum()

        for iw in range(len(lambda_bands_cm)):
            z11_avg = np.sum(zscat[:, iw, :, 0] * distri[:, None], axis=0)
            z12_avg = np.sum(zscat[:, iw, :, 1] * distri[:, None], axis=0)

            f_z11 = interp1d(theta, z11_avg, kind='linear', fill_value="extrapolate")
            f_z12 = interp1d(theta, z12_avg, kind='linear', fill_value="extrapolate")

            Z11eff[ia, iw, :] = f_z11(scatt_angles_deg)
            Z12eff[ia, iw, :] = f_z12(scatt_angles_deg)

    return -Z12eff / Z11eff




def compute_P_and_omega_vs_amax(lambda_bands_cm, a_max_list, zscat, res, theta):
    i90 = theta.searchsorted(90)
    distri_template = a_max_list**0.5
    distri_template /= distri_template.sum()

    Z11eff = np.zeros((len(a_max_list), len(lambda_bands_cm)))
    Z12eff = np.zeros_like(Z11eff)
    omega = np.zeros_like(Z11eff)

    for ia, amax in enumerate(a_max_list):
        distri = distri_template.copy()
        distri[a_max_list > amax] = 0
        distri /= distri.sum()

        Z11eff[ia] = np.sum(zscat[:, :, i90, 0] * distri[:, None], axis=0)
        Z12eff[ia] = np.sum(zscat[:, :, i90, 1] * distri[:, None], axis=0)

        k_a = np.sum(res['k_abs'] * distri[:, None], axis=0)
        k_s = np.sum(res['k_sca'] * distri[:, None], axis=0)
        omega[ia] = k_s / (k_a + k_s)

    P = np.abs(Z12eff / Z11eff)
    
    P_times_omega = P * omega
    
    return P, omega, P_times_omega





# Run DSHARP
# --------------------------------------------------------------------------------------
def run_DSHARP(f, a_max_test_cm, a_max_dist_cm, lambda_bands_cm, lambda_dist_cm):
    
    oc, rho_g_cm3 = do.get_dsharp_mix(porosity = 1 - f)
    
    res = do.get_opacities(a_max_test_cm, lambda_dist_cm, rho_g_cm3, oc)
    
    k_abs  = res['k_abs']
    k_sca  = res['k_sca']
    S1     = res['S1']
    S2     = res['S2']
    theta  = res['theta']

    mass_grams = calculate_grain_mass(a_max_dist_cm, rho_g_cm3, "cm")
    
    res_scatter = do.get_opacities(a_max_dist_cm, lambda_bands_cm, rho_g_cm3, oc, n_angle=100)

    zscat = do.calculate_mueller_matrix(lambda_bands_cm, mass_grams, 
                                    res_scatter['S1'], res_scatter['S2'], 
                                    theta=res_scatter['theta'], 
                                    k_sca=res_scatter['k_sca'])['zscat']
    theta = res_scatter['theta']
    
    P, omega, P_times_omega = compute_P_and_omega_vs_amax(lambda_bands_cm, a_max_dist_cm, zscat, res_scatter, theta)

    
    return P, omega, P_times_omega
# --------------------------------------------------------------------------------------


