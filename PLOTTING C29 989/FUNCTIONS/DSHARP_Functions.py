import numpy as np 
import miepython
from scipy.interpolate import interp1d
import dsharp_opac as do


import pandas as pd

from DustModelFunctions import * 



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
            


    return  -Z12eff / Z11eff




def compute_P_and_omega_vs_amax(lambda_bands_cm, a_max_list, zscat, res, theta):
    i90 = theta.searchsorted(90)
    
    print(rf'in compute_P_and_omega_vs_amax: theta[i90] = {theta[i90]}')
    
    
    distri_template = a_max_list**0.5
    distri_template /= distri_template.sum()

    Z11eff = np.zeros((len(a_max_list), len(lambda_bands_cm)))
    Z12eff = np.zeros_like(Z11eff)
    omega = np.zeros_like(Z11eff)
    omega_eff = np.zeros_like(Z11eff)

    for ia, amax in enumerate(a_max_list):
        distri = distri_template.copy()
        distri[a_max_list > amax] = 0
        distri /= distri.sum()

        Z11eff[ia] = np.sum(zscat[:, :, i90, 0] * distri[:, None], axis=0)
        Z12eff[ia] = np.sum(zscat[:, :, i90, 1] * distri[:, None], axis=0)

        g = np.sum(res['g'] * distri[:, None], axis=0)
        
        k_a = np.sum(res['k_abs'] * distri[:, None], axis=0)
        k_s = np.sum(res['k_sca'] * distri[:, None], axis=0)
        #k_s_eff = (1-g) * k_s
        k_s_eff = np.sum((1 - res['g']) * res['k_sca'] * distri[:, None], axis=0)
        
        omega[ia] = k_s / (k_a + k_s)
        omega_eff[ia] = k_s_eff / (k_a + k_s_eff)
        

    P = np.abs(Z12eff / Z11eff)
    
#     print("\nFigure 3 P values:")
#     print(P[:, 1])   # choose one wavelength index
     
    P_times_omega = P * omega
    P_times_omega_eff = P * omega_eff
    
    
    return P, omega, omega_eff, P_times_omega, P_times_omega_eff


    
    
    

# Run DSHARP
# --------------------------------------------------------------------------------------
def run_DSHARP(f, a_max_test_cm, a_max_dist_cm, lambda_bands_cm, lambda_dist_cm, debug_label=None):
    


    oc, rho_g_cm3 = do.get_dsharp_mix(porosity = 1 - f)
    
    
    print("rho =", rho_g_cm3)
    print("Band wavelengths (cm):", lambda_bands_cm)
    
    res = do.get_opacities(a_max_test_cm, lambda_dist_cm, rho_g_cm3, oc)
    
    k_abs  = res['k_abs']
    k_sca  = res['k_sca']
    S1     = res['S1']
    S2     = res['S2']
    theta  = res['theta']

    mass_grams = calculate_grain_mass(a_max_dist_cm, rho_g_cm3, "cm")
    
    res_scatter = do.get_opacities(a_max_dist_cm, lambda_bands_cm, rho_g_cm3, oc, n_angle=100)
    
    print("S1 shape:", res_scatter["S1"].shape)
    print("S2 shape:", res_scatter["S2"].shape)
    print("theta shape:", res_scatter["theta"].shape)

    print("Any NaNs in S1?", np.isnan(res_scatter["S1"]).any())
    print("Any NaNs in S2?", np.isnan(res_scatter["S2"]).any())
    print("Max |S1|:", np.nanmax(np.abs(res_scatter["S1"])))
    print("Max |S2|:", np.nanmax(np.abs(res_scatter["S2"])))


    zscat = do.calculate_mueller_matrix(lambda_bands_cm, mass_grams, 
                                    res_scatter['S1'], res_scatter['S2'], 
                                    theta=res_scatter['theta'], 
                                    k_sca=res_scatter['k_sca'])['zscat']
    
    print("zscat shape:", zscat.shape)
    print("Any NaNs?", np.isnan(zscat).any())
    print("Max Z11:", np.nanmax(zscat[...,0]))
    print("Min Z11:", np.nanmin(zscat[...,0]))

    theta = res_scatter['theta']
    
    
    # ----------------------------------------------------
# DEBUG: inspect single-grain Mueller matrix
# ----------------------------------------------------

    iw = 1      # Band 5

    # Print summary for several grain sizes
    for target in [1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 1]:

        igrain = np.argmin(np.abs(a_max_dist_cm - target))

        Z11 = zscat[igrain, iw, :, 0]
        Z12 = zscat[igrain, iw, :, 1]
        P = -Z12 / Z11

        print(
            f"target={target:.3e} cm   "
            f"actual={a_max_dist_cm[igrain]:.3e} cm   "
            f"P90={P[np.argmin(np.abs(theta-90))]:.3f}   "
            f"Pmax={np.nanmax(P):.3f}   "
            f"Pmin={np.nanmin(P):.3f}"
        )


    # ---------------- Plot ONE grain ----------------

    target = 3e-1                      # choose which grain to inspect
    igrain = np.argmin(np.abs(a_max_dist_cm - target))

    Z11 = zscat[igrain, iw, :, 0]
    Z12 = zscat[igrain, iw, :, 1]
    P_single = -Z12 / Z11

    print("\nSingle grain")
    print("----------------")
    print(f"Requested amax = {target:.3e} cm")
    print(f"Nearest grid point = {a_max_dist_cm[igrain]:.3e} cm")

    i90 = np.argmin(np.abs(theta - 90))

    print(f"P(90°) = {P_single[i90]:.4f}")
    print(f"Pmax   = {np.nanmax(P_single):.4f}")
    print(f"Pmin   = {np.nanmin(P_single):.4f}")
    print(f"theta(Pmax) = {theta[np.nanargmax(P_single)]:.1f} deg")

    import matplotlib.pyplot as plt

    plt.figure(figsize=(7,5))

    plt.plot(theta, Z11 / np.max(Z11), label="Z11 / max(Z11)")
    plt.plot(theta, -Z12 / np.max(np.abs(Z12)), label="-Z12 / max(|Z12|)")
    plt.plot(theta, P_single, linewidth=3, label="P = -Z12/Z11")

    plt.axvline(90, color="k", linestyle=":")

    plt.xlabel("Scattering angle (deg)")
    plt.ylabel("Normalized value")
    plt.title(f"Single grain: amax = {a_max_dist_cm[igrain]:.3e} cm")
    plt.legend()
    plt.grid()
    plt.show()
    
    
    
    P, omega, omega_eff, P_times_omega, P_times_omega_eff = compute_P_and_omega_vs_amax(lambda_bands_cm, a_max_dist_cm, zscat, res_scatter, theta)

#     if debug_label is not None:
#     debug_print_DSHARP_inputs_outputs(
#         debug_label,
#         f,
#         a_max_test_cm,
#         a_max_dist_cm,
#         res_scatter,
#         zscat,
#         P,
#         omega,
#         P_times_omega
#     )
    
    return P, omega, omega_eff, P_times_omega, P_times_omega_eff
# --------------------------------------------------------------------------------------


