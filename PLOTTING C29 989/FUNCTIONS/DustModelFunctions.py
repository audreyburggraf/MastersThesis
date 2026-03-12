import numpy as np 
import miepython
from scipy.interpolate import interp1d
import dsharp_opac as do









# Function to convert micron to cm
# --------------------------------------------------------------------------------------
def micron_to_cm(value_micron):
    value_micron = np.array(value_micron)  # convert input to numpy array
    return value_micron * 1e-4  # 1 micron = 1e-4 cm
# --------------------------------------------------------------------------------------

# Function to convert mm to cm
# --------------------------------------------------------------------------------------
def mm_to_cm(value_mm):
    value_mm = np.array(value_mm)  # convert input to numpy array
    return value_mm * 1e-1  # 1 mm = 0.1 cm = 1e-1 cm
# --------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------
def calculate_grain_mass(a, density_g_cm3, a_units):
    if a_units == "cm":
        a_cm = a
    elif a_units == "micron":
        a_cm = micron_to_cm(a)
    else:
        raise ValueError("Only 'cm' and 'micron' units are supported for grain size")
    
    # Calculate the volume of the sphere in cm³
    volume_cm3 = (4/3) * np.pi * a_cm**3
    
    # Calculate the mass in grams (mass = density × volume)
    mass_grams = density_g_cm3 * volume_cm3
    
    return mass_grams
# --------------------------------------------------------------------------------------


# Size distribution of n(a) prop to a^-3.5 from Kataoka et al, 2015
# --------------------------------------------------------------------------------------
# def grain_size_distribution(xmin, xmax, N):
#     p = 3.5 
    
#     # Normalization constant for f(x) = C * x^-p over [xmin, xmax]
#     C = (1 - p) / (xmax**(1 - p) - xmin**(1 - p))
    
#     # Define x values
#     x = np.linspace(xmin, xmax, N)

#     # Calculate PDF
#     pdf = C * x**(-p)
    
#     # Normalize weights to sum to 1 (optional but useful)
#     weights = pdf / np.sum(pdf)
    
#     return x, pdf, weights
# --------------------------------------------------------------------------------------
# def sample_grain_sizes(a_min, a_max, N, p=3.5):
#     # Power-law exponent related term
#     inv_exp = 1 - p
    
#     # Uniform random numbers from 0 to 1
#     u = np.random.uniform(0, 1, N)
    
#     # Inverse CDF transform to get grain sizes distributed as a^-p
#     sizes = (u * (a_max**inv_exp - a_min**inv_exp) + a_min**inv_exp) ** (1 / inv_exp)
    
#     return sizes
# --------------------------------------------------------------------------------------

# Function to the size parameter of x = 2*pi*a/lambda 
# # --------------------------------------------------------------------------------------
# def size_parameter(a_micron, lambda_micron):
    
#     # Note: This can read in any units, just as long as the units of a and lambda are the same
    
#     x_micron = (2 * np.pi * a_micron**2) / lambda_micron 
    
#     return x_micron
# # --------------------------------------------------------------------------------------






# This function is taken from Kataoka et al. 2015, page 2
# --------------------------------------------------------------------------------------
def calculate_omega(kappa_abs, kappa_sca):
    
    omega = kappa_sca / (kappa_abs + kappa_sca)
    
    return omega
# --------------------------------------------------------------------------------------






# # Calcluate degree of polarization
# # --------------------------------------------------------------------------------------
# # From Kataoka et al, 2015 they find Z11 and Z12 using a scattering matrix
# # --------------------------------------------------------------------------------------
# def calculate_P_scatt_angle(lambda_micron, scatt_angles_deg_array, a_max_micron, a_min_micron=0.001, N=100):
#     # Grain size array (microns)
#     grain_sizes_micron = np.logspace(np.log10(a_min_micron), np.log10(a_max_micron), N)

#     # Compute PDF from grain size distribution and normalize
#     grainsize_pdf = sample_grain_sizes(a_min_micron, a_max_micron, N, p=3.5)

#     # Empty arrays for Z11 and Z12
#     Z11_total = np.zeros_like(scatt_angles_deg_array, dtype=float)
#     Z12_total = np.zeros_like(scatt_angles_deg_array, dtype=float)

#     # Refractive index from Lin et al. (2023)
#     m = find_complex_m_average(lambda_micron)

#     for j, a in enumerate(grain_sizes_micron):
#         x = size_parameter(a, lambda_micron)

#         for i, theta_deg in enumerate(scatt_angles_deg_array):
#             theta_rad = np.deg2rad(theta_deg)
#             Z = miepython.phase_matrix(m, x, theta_rad)

#             Z11_total[i] += Z[0][0] 
#             Z12_total[i] += Z[0][1] 

#     # Normalize (optional, for visualization)
#     Z11_total /= np.sum(Z11_total)
#     Z12_total /= np.sum(Z11_total)

#     pol_deg = -Z12_total / Z11_total

#     return pol_deg
# # --------------------------------------------------------------------------------------
# def calculate_P_grainsize(lambda_micron, scatt_angle_deg, a_max_micron_array, a_min_micron=0.001, N=100):
#     m = find_complex_m_average(lambda_micron)

#     # Convert angle to radians
#     theta_rad = np.deg2rad(scatt_angle_deg)

#     # Output array
#     pol_deg_array = np.zeros_like(a_max_micron_array, dtype=float)

#     for i, a_max in enumerate(a_max_micron_array):
#         grainsize_pdf = sample_grain_sizes(a_min_micron, a_max_micron, N, p=3.5)


#         Z11_total = 0.0
#         Z12_total = 0.0

#         for a, weight in zip(grain_sizes, weights):
#             x = size_parameter(a, lambda_micron)
            
#             Z = miepython.phase_matrix(m, x, theta_rad)
            
#             Z11_total += Z[0][0] 
#             Z12_total += Z[0][1] 
#         pol_deg_array[i] = -Z12_total / Z11_total if Z11_total != 0 else 0.0

#     return pol_deg_array
# # --------------------------------------------------------------------------------------






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




# sum of squares function 
# --------------------------------------------------------------------------------------
# The function we want to minimize is the sum of squares:
# sum of squares = $\sum_{i=1}^n(X_i - X)^2$, where

# $X_i$ is the ith item

# $X$ is the mean 

# so, for this scenario $X_i$ is the model $P \omega$ times some scaling factor, and $X$ is the POLF value. 

# We want to take the sum of squares for all bands, and find the scaling factor that minimizes it.
# --------------------------------------------------------------------------------------
def dust_sum_of_sq(sf, model_Pw, obs_POLF):
    total = 0.0
    
    for i in range(len(model_Pw)):
        diff = np.sum(sf * model_Pw[i] - obs_POLF[i]) ** 2
        total += diff
    
    return total
# --------------------------------------------------------------------------------------







def find_sf_v2(a_max_f_dist_micron, P_times_omega, df_POLF, print_results=False):

    # We are looking at the POLF from the mean/average method, none others 
    POLF_obs = df_POLF['POLF_Gaussian']

    # Make arrays to store the median and standard deviation values 
    sf_medians = []
    sf_stds = []

    # Loop over each a_max grid point
    for i in range(len(a_max_f_dist_micron)):

        Pw_model = P_times_omega[i, :]  # 4 bands
        
        # Avoid divide by zero
        valid = Pw_model > 0
        
        sf_i = POLF_obs[valid] / Pw_model[valid]

        sf_medians.append(np.median(sf_i))
        sf_stds.append(np.std(sf_i))

    sf_medians = np.array(sf_medians)
    sf_stds = np.array(sf_stds)

    # Choose best a_max where scatter is minimized
    best_idx = np.argmin(sf_stds)

    best_a_max_f = a_max_f_dist_micron[best_idx]
    best_sf = sf_medians[best_idx]

    # Compute predicted POLF at best solution
#     best_POLF = P_times_omega[best_idx, :] * best_sf

    if print_results:
        print("Best a_max:", best_a_max_f)
        print("Best SF:", best_sf)

    # return values
    return best_a_max_f, POLF_obs, best_sf

