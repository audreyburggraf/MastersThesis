# ------------------------------------------------------------
# Import the necessary packages
# ------------------------------------------------------------
import glitterin
import numpy as np
import matplotlib.pyplot as plt
from UnitConversion import *
from glitterin_LabDistributionFunc import *
from tqdm import tqdm

from DustModelJustAMAX import *
# ------------------------------------------------------------
#
#
#
#
# ------------------------------------------------------------
# Get the producer
# ------------------------------------------------------------
# Make the path to the neural network models 
nndir = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/glitterin/nnmodels/"

# This producer is an object that can compute scattering properties using the neural network models in 
# the nndir folder 

# Create the scattering wrapper
producer = glitterin.user.ScatteringProducer(nndir=nndir)
from glitterin.user import extrapolate_small_a
# ------------------------------------------------------------
#
#
#
#
# ------------------------------------------------------------
# Function to get r from x using Equation 9 from Lin 2025
# ------------------------------------------------------------
def r_from_x(x, w):
    
    # w is wavelength 
    
    return (x * w) / (2 * np.pi)
# ------------------------------------------------------------
# Function to get x from r using Equation 9 from Lin 2025
# ------------------------------------------------------------
def x_from_r(r, w):
    
    # w is wavelength 
    
    return (2 * np.pi * r) / w
# ------------------------------------------------------------
#
#
#
#
# ------------------------------------------------------------
# Function to get r_vol from r_enc using Equation 6from Lin 2025
# ------------------------------------------------------------
def rvol_from_renc(r_enc, f = 0.236):
    
    
    
    return (f * r_enc**3)**(1/3)
# ------------------------------------------------------------
#
#
#
#
# ------------------------------------------------------------
# def CalculateNormalizedZ_fig10(result, theta):
    
#     # 1. Extraction (Using dynamic keys to save space)
#     Z11 = result['Z11']
    
#     idx_90 = np.argmin(np.abs(theta - 90))

#     Z11_norm = result['Z11']/result['Z11'][idx_90]
#     Z12 = result['Z12']
#     Z22 = result['Z22']
#     Z33 = result['Z33']
#     Z34 = result['Z34']
#     Z44 = result['Z44']
    
#     # 2. Calculations
# #     N11 = (2 * np.pi * Z11) # / result['Csca']
# #     print('')
#     N12 = -Z12 / Z11_norm
#     N22 = Z22 / Z11_norm
#     N33 = Z33 / Z11_norm
#     N34 = Z34 / Z11_norm
#     N44 = Z44 / Z11_norm  # Fixed the typo here
    
#     # 3. Saving Results
    
# #     result['N11_calc'] = N11 
#     result['N12_fig10'] = N12 
#     result['N22_fig10'] = N22 
#     result['N33_fig10'] = N33 
#     result['N34_fig10'] = N34 
#     result['N44_fig10'] = N44
# ------------------------------------------------------------
# ------------------------------------------------------------
# Get the normalized quantities from Lin 2025 Equations 5
# Note: Used Google gemini to help with the suffix part
# ------------------------------------------------------------
def CalculateNormalizedZ(result):
    
    # 1. Extraction (Using dynamic keys to save space)
    Z11 = result['Z11_eff']
    Z12 = result['Z12_eff']
    Z22 = result['Z22_eff']
    Z33 = result['Z33_eff']
    Z34 = result['Z34_eff']
    Z44 = result['Z44_eff']
    
    Z11_safe = np.where(Z11 == 0, 1e-12, Z11)
    
    # 2. Calculations
#     N11 = (2 * np.pi * Z11) # / result['Csca']
#     print('')
    N12 = -Z12 / Z11_safe
    N22 = Z22 / Z11_safe
    N33 = Z33 / Z11_safe
    N34 = Z34 / Z11_safe
    N44 = Z44 / Z11_safe 
    
    # 3. Saving Results
    
#     result['N11_calc'] = N11 
    result['N12_eff'] = N12 
    result['N22_eff'] = N22 
    result['N33_eff'] = N33 
    result['N34_eff'] = N34 
    result['N44_eff'] = N44
# ------------------------------------------------------------

#
#
#
#
# ------------------------------------------------------------
# Get n and k from wavelength
# ------------------------------------------------------------
def get_Re_Im(w, w_units):
    # Note: This is only using astronomical silicate 
    
    # w is the wavelength
    # w is the units 
    
    
    if w_units == 'micron':
        w_micron = w 
    
    elif w_units == 'mm':
        w_micron = mm_to_micron(w)
        
    elif w_units == 'cm':
        w_micron = cm_to_micron(w)
    
    else:
        return print("The only units allowed are 'micron', 'mm', 'cm'")
    
    
    # -----------------------------
    # Load file (skip header lines)
    # -----------------------------
    data = np.loadtxt('callindex.out_silD03.txt', skiprows=6)

    lam = data[:, 0]           # wavelength [micron]
    Re_n = data[:, 3] + 1      # Re(n) = (Re(n)-1) + 1
    Im_n = data[:, 4]          # Im(n)
    
    
    # -----------------------------
    # Sort wavelength (IMPORTANT)
    # -----------------------------
    idx = np.argsort(lam)
    lam = lam[idx]
    Re_n = Re_n[idx]
    Im_n = Im_n[idx]

    # -----------------------------
    # Interpolate to your wavelength
    # -----------------------------
    n = np.interp(w_micron, lam, Re_n)
    k = np.interp(w_micron, lam, Im_n)

    return n, k
# ------------------------------------------------------------
#
#
#
#
# ------------------------------------------------------------
# Get the result
# ------------------------------------------------------------
def run_glitterin(x_vol, # Size parameter x = 2 pi a / lambda 
                  re_m,  # Refractive index m = n + ik 
                  im_m, # re_m = n and im_m = k
                  theta,  # Scatting angles 
                  w,  # Wavelength [cm]
                  xtype = 'vol',
                  outtype = 'dict'):



    #print("This function currently only works for xtype = 'xvol'")
    #print(' ')
    
    producer.setup(['Cext', 'Cabs', 'Csca', 'albedo', 'Z11', 'Z12', 'Z22', 'Z33', 'Z34', 'Z44', 
                                              'N11', 'N12', 'N22', 'N33', 'N34', 'N44'])
    
    result = producer(x_vol, re_m, im_m, theta, w, xtype = xtype, outtype = outtype)
    
    # This is code from glitterin intro files:
    # ---------------------------------------------------------------------------------------
    # this is the volume equivalent grain size 
    # This converts the size parameter (x) to the physical grain size (a)
    avol = w / 2 / np.pi * x_vol

    # These are the Q efficiency factors. 
    # It's easier in the user side to explicitly determine what kind of radius we are using for the efficiency. 
    # For radiation transfer, it doesn't matter. Only the actual Cext, etc, matters
    # Q = C = pi a^2
    result['Q_ext_vol'] = result['Cext'] / (np.pi * avol**2)
    result['Q_abs_vol'] = result['Cabs'] / (np.pi * avol**2)

    # This tells us what range of grain sizes is safe to use this model
    # The minimum xenc of the training data is 0.1
    # We can figure result4 what the volume equivalent x that is with
#     min_x_vol = glitterin.user.xvol_from_xenc(0.1)

    # the maximum xvol of the training data depends on the refractive index
#     max_x_vol = glitterin.user.max_xvol_from_nk(re_m[0], im_m[0])

#     print(rf'The range of grain sizes that is safe to use is {min_x_vol:.2f} to {max_x_vol:.2f}')

#     a_min = w / (2 * np.pi) * min_x_vol
#     a_max = w / (2 * np.pi) * max_x_vol

#     print(f"Safe grain size range: {a_min} to {a_max} cm")
    # ---------------------------------------------------------------------------------------

    # Change 'albedo' to 'omega'
    result['omega'] = result['albedo']
    
    # Scattering cross-section from page 3 of the paper
#     result['Csca'] = result['Cext'] - result['Cabs'] 
    
    # epsilon from equation 4 of the paper (page 3)
    result['epsilon'] = 1 - result['omega'] 
    
    # Get the r values from x and wavelength, using equation 9 from the paper
    # Note: w is in cm
    result['r_vol'] = (x_vol * w) / (2 * np.pi)
    
    # CHECK THIS
    result['Q_sca_vol'] = result['Q_ext_vol'] - result['Q_abs_vol']
    
    # kappa values using equation 20 from the paper
    rho_s = 1.675 # [grams / cm^3]
    result['k_ext'] = result['Cext'] / (4/3 * np.pi * rho_s * result['r_vol']**3)
    result['k_abs'] = result['Cabs'] / (4/3 * np.pi * rho_s * result['r_vol']**3)
    result['k_sca'] = result['Csca'] / (4/3 * np.pi * rho_s * result['r_vol']**3)
    
    # Looking at enc
    # ---------------------------------------------------------------------------------------
    # Get r_enc from r_vol using Equation 6 from Lin 2025
    f = 0.236 # Constant from Lin 2025, based on Zubko 2024
    result['r_enc'] = (result['r_vol']**3 / f)**(1/3)
    
    # Get x_enc from r_enc using Equation 9 in Lin 2025
    result['x_enc'] = (result['r_enc'] * 2 * np.pi) / w
    # ---------------------------------------------------------------------------------------
    
    # Looking at pja 
    # ---------------------------------------------------------------------------------------
    # Get r_pja from r_enc using Equation 9 from Lin 2025
    g = 0.611 # Constant from Lin 2025, based on Zubko 2024
    result['r_pja'] = (g * result['r_enc']**2)**(1/2)
    
    # Using Equation 8 from Lin 2025 (https://arxiv.org/pdf/2511.09668)
    # C_i = Q_i^j pi r_j^2,
    # where i can be ext, abs, sca and j can be vol, pja, enc
    result['Q_ext_pja'] = result['Cext'] / (np.pi * result['r_pja']**2)
    result['Q_abs_pja'] = result['Cabs'] / (np.pi * result['r_pja']**2)
    result['Q_sca_pja'] = result['Q_ext_pja'] - result['Q_abs_pja']
    # ---------------------------------------------------------------------------------------
    
    
    add_P_to_result(result, theta)
    
    
    #CalculateNormalizedZ_fig10(result, theta)
    

#     print(result.keys())
    

    return result
# ------------------------------------------------------------
#
#
#
#
#
# ------------------------------------------------------------
# Get effective Z values and P
# ------------------------------------------------------------
# def add_P_to_result(result, w, theta):
    
#     # Start the distriution
#     r = result['r_vol']
    
#     weights = r**(-2.5) 

#     # Normalize
#     weights /= weights.sum()
    
#     print(np.sum(weights))
#     print(np.max(weights))
    
# #     if show_weight_plot:
# #         # Plot the weights
# #         plt.figure()
# #         plt.loglog(r, weights / weights.max())
# #         plt.xlabel('r [cm]')
# #         plt.ylabel('normalized weights')
# #         plt.title('Shape of distribution')
# #         plt.show()
    

#     # Get the effective matrix elements 
#     # Summing the grain grain size * n(r)
#     result['Z11_eff'] = np.sum(result['Z11'] * weights, axis=1)
#     result['Z12_eff'] = np.sum(result['Z12'] * weights, axis=1)
#     result['Z22_eff'] = np.sum(result['Z22'] * weights, axis=1)
#     result['Z33_eff'] = np.sum(result['Z33'] * weights, axis=1)
#     result['Z34_eff'] = np.sum(result['Z34'] * weights, axis=1)
#     result['Z44_eff'] = np.sum(result['Z44'] * weights, axis=1)
    
   
#     CalculateNormalizedZ(result)
    
    
#     # Calculate P 
#     result['P'] = - result['Z12_eff'] / result['Z11_eff']
    
#     # Get the effective kappa values
#     result['k_abs_eff'] = np.sum(result['k_abs'] * weights)
#     result['k_sca_eff'] = np.sum(result['k_sca'] * weights)
    
#     # Calculate effective omega
#     result['omega_eff'] = result['k_sca_eff'] / (result['k_abs_eff'] + result['k_sca_eff'])
    
#     # Get the index of the 90 degree value in theta
#     idx_90 = np.argmin(np.abs(theta - 90))
#     print(rf'Sanity check: theta[idx_90] = {theta[idx_90]}')
    
#     # Get P_90
#     result['P_90'] = result['P'][idx_90]
    
#     # calculate P times omega
#     result['P_omega_90'] = result['P_90'] * result['omega_eff']

def add_P_to_result(result, theta):

    r = result['r_vol']
    idx_90 = np.argmin(np.abs(theta - 90))

    base_weights = r**(-3)

    N_amax = len(r)
    N_theta = len(theta)

    # Allocate all effective Z's
    result['Z11_eff'] = np.zeros((N_amax, N_theta))
    result['Z12_eff'] = np.zeros((N_amax, N_theta))
    result['Z22_eff'] = np.zeros((N_amax, N_theta))
    result['Z33_eff'] = np.zeros((N_amax, N_theta))
    result['Z34_eff'] = np.zeros((N_amax, N_theta))
    result['Z44_eff'] = np.zeros((N_amax, N_theta))
    
    

    # Scalars vs a_max
    result['k_abs_eff'] = np.zeros(N_amax)
    result['k_sca_eff'] = np.zeros(N_amax)
    result['omega_eff'] = np.zeros(N_amax)
    result['P_90'] = np.zeros(N_amax)

    for i, rmax in enumerate(r):

        # --- weights ---
        weights = base_weights.copy()
        weights[r > rmax] = 0
        weights /= weights.sum()

        # --- effective Z ---
        result['Z11_eff'][i] = np.sum(result['Z11'] * weights, axis=1)
        result['Z12_eff'][i] = np.sum(result['Z12'] * weights, axis=1)
        result['Z22_eff'][i] = np.sum(result['Z22'] * weights, axis=1)
        result['Z33_eff'][i] = np.sum(result['Z33'] * weights, axis=1)
        result['Z34_eff'][i] = np.sum(result['Z34'] * weights, axis=1)
        result['Z44_eff'][i] = np.sum(result['Z44'] * weights, axis=1)
        
       

        # --- polarization curve ---
        P_theta = - result['Z12_eff'][i] / result['Z11_eff'][i]
        result['P_90'][i] = P_theta[idx_90]

        # --- opacities ---
        k_abs_eff = np.sum(result['k_abs'] * weights)
        k_sca_eff = np.sum(result['k_sca'] * weights)

        result['k_abs_eff'][i] = k_abs_eff
        result['k_sca_eff'][i] = k_sca_eff
        result['omega_eff'][i] = k_sca_eff / (k_abs_eff + k_sca_eff)
        
    CalculateNormalizedZ(result)

    # --- final combined quantity ---
    result['P_omega_90'] = result['P_90'] * result['omega_eff']

    return result
    
# ------------------------------------------------------------
# Get P*omega
# ------------------------------------------------------------
# ChatGPT helped convert my DSHARP function to get Pw to get it from glitterin
# def add_P_to_results_all(results_all, wavelength_labels, theta):

#     i90 = np.argmin(np.abs(theta - 90)) # Note: I checked this and this gives "theta[i90] = 90.0"

#     for label in wavelength_labels:

#         res = results_all[label]

#         a = res['r_vol']   # your grain sizes (this is your x-axis!)
       

#         # distribution template
#         distri_template = a**(0.5)
#         distri_template /= distri_template.sum()
        

#         Z11_all = res['Z11']
#         Z12_all = res['Z12']

#         # normalize by scattering cross section (or k_sca)
#         Z11_norm = Z11_all / res['k_sca']
#         Z12_norm = Z12_all / res['k_sca']

#         Z11 = Z11_norm[i90, :]
#         Z12 = Z12_norm[i90, :]
        
#         P = np.zeros(len(a))
#         omega = np.zeros(len(a))
#         k_abs_eff = np.zeros(len(a))
#         k_sca_eff = np.zeros(len(a))

#         for ia, amax in enumerate(a):

#             distri = distri_template.copy()
#             distri[a > amax] = 0
#             distri /= distri.sum()


#             k_a = np.sum(res['k_abs'] * distri)
#             k_s = np.sum(res['k_sca'] * distri)

#             omega[ia] = k_s / (k_a + k_s)
            
#             P_per_a = -Z12 / Z11   # polarization per grain size

#             P[ia] = np.sum(P_per_a * distri)
            
#             k_abs_eff[ia] = np.sum(res['k_abs'] * weights) / np.sum(weights)
#             k_sca_eff[ia] = np.sum(res['k_sca'] * weights) / np.sum(weights)

#         # store back into your dictionary
#         res['P'] = P
#         res['omega_dist'] = omega
#         res['P_omega'] = P * omega
#         res['k_abs_dist'] = k_abs_eff
#         res['k_sca_dist'] = k_sca_eff
# ------------------------------------------------------------
def add_P_to_results_all(results_all, wavelength_labels, theta):

    i90 = np.argmin(np.abs(theta - 90))

    for label in wavelength_labels:

        res = results_all[label]
        a = res['r_vol']

        # distribution (mass-weighted MRN)
        distri_template = a**(-0.5)
        distri_template /= distri_template.sum()

        # Mueller elements (already per grain size)
        Z11 = res['Z11'][i90, :]
        Z12 = res['Z12'][i90, :]

        # polarization per grain
        P_per_a = -Z12 / Z11

        P = np.zeros(len(a))
        omega = np.zeros(len(a))
        k_abs_eff = np.zeros(len(a))
        k_sca_eff = np.zeros(len(a))

        for ia, amax in enumerate(a):

            distri = distri_template.copy()
            distri[a > amax] = 0
            distri /= distri.sum()

            k_a = np.sum(res['k_abs'] * distri)
            k_s = np.sum(res['k_sca'] * distri)

            omega[ia] = k_s / (k_a + k_s)

            # polarization averaging
            P[ia] = np.sum(P_per_a * distri)

            # opacity averaging
            k_abs_eff[ia] = k_a
            k_sca_eff[ia] = k_s

        res['P'] = P
        res['omega_dist'] = omega
        res['P_omega'] = P * omega
        res['k_abs_dist'] = k_abs_eff
        res['k_sca_dist'] = k_sca_eff
# ------------------------------------------------------------










# Function to run glitterin for the lab distribution
# ------------------------------------------------------------
def RunGlitterinLabDist(particle):
    # Path to the lab data text files
    particle_file_loc = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/glitterin/DustModels/"
    
    
    # Set things based on the particle type
    # -------------------------------
    if particle == 'feldspar':
        fname = 'sizes_feldspar.txt'
        w_um = 0.4416 # [micron]
        n = 1.5
        k = 1e-4

    # -------------------------------
    elif particle == 'hematite':
        fname = 'size_hematite.txt'
        w_um = 0.6328 # [micron]
        n = 3
        k = 0.03

    # -------------------------------
    else:
        return print("The function currently only works for the particle 'feldspar' or 'hematite'")
    # -------------------------------
    

    # Get the particle size
    sizes = read_gals_size_table(particle_file_loc + fname, skiprows=2, out_type='container')
    lab_x = sizes.x(w_um)
    
    # prepare inputs for the neural network
    x_inp = lab_x
    n_inp = n + np.zeros_like(x_inp)
    k_inp = k + np.zeros_like(x_inp)
    w_inp = w_um*1e4 + np.zeros_like(x_inp) # [cm]
    theta = np.linspace(0, 180, 181)
    

    # Now we can produce the scattering quantities for each grain size
    producer.setup(['Cext', 'Cabs', 'Csca', 'albedo', 'Z11', 'Z12', 'Z22', 'Z33', 'Z34', 'Z44', 
                                              'N11', 'N12', 'N22', 'N33', 'N34', 'N44'])
    
    mat = producer(x_inp, n_inp, k_inp, theta, w_inp, xtype='enc', outtype='MuellerMatrix')
    
    


        
    # Make an additional place to store the values averaged over the distribution
    mat_dist = {}

    mat_dist['Cext'] = average_over_size(mat.Cext, sizes.N_log_r, sizes.d_log_r)

    mat_dist['Z11'] = average_over_size(mat.Z11.T, sizes.N_log_r, sizes.d_log_r)
    mat_dist['Z12'] = average_over_size(mat.Z12.T, sizes.N_log_r, sizes.d_log_r)
    mat_dist['Z22'] = average_over_size(mat.Z22.T, sizes.N_log_r, sizes.d_log_r)
    mat_dist['Z33'] = average_over_size(mat.Z33.T, sizes.N_log_r, sizes.d_log_r)
    mat_dist['Z34'] = average_over_size(mat.Z34.T, sizes.N_log_r, sizes.d_log_r)
    mat_dist['Z44'] = average_over_size(mat.Z44.T, sizes.N_log_r, sizes.d_log_r)
    

    # 3. Saving Results
    
#     result['N11_calc'] = N11 
    mat_dist['N12'] = - mat_dist['Z12'] / mat_dist['Z11']
    mat_dist['N22'] = mat_dist['Z22'] / mat_dist['Z11']
    mat_dist['N33'] = mat_dist['Z33'] / mat_dist['Z11']
    mat_dist['N34'] = mat_dist['Z34'] / mat_dist['Z11']
    mat_dist['N44'] = mat_dist['Z44'] / mat_dist['Z11']

    return mat, mat_dist, theta
# ------------------------------------------------------------
#
#
#
#
#
# -----------------------------------------------------------
# Function to run glitterin using the DSHARP size distrbution
# ------------------------------------------------------------
# def make_dsharp_sizes(amax_um, amin_um=0.1, Nsize=200):
#     sizes = SizeContainer()

#     # uniform grid in log10(r)
#     sizes.log_r = np.linspace(np.log10(amin_um), np.log10(amax_um), Nsize)

#     # r in micron
#     r = sizes.r

#     # DSHARP-style n(r) ∝ r^-3.5
#     # convert to N_log_r ∝ r * n(r) ∝ r^-2.5
#     sizes.N_log_r = r**(-2.5)

#     # not needed for averaging, but can fill with dummy arrays
#     sizes.S_log_r = np.zeros_like(r)
#     sizes.V_log_r = np.zeros_like(r)

#     return sizes
# ------------------------------------------------------------
def RunGlitterinPowerLawDist(amax_cm, w, w_units, n, k, amin_cm = 1e-5, Nsize = 200):
    
    # Note: I want everything in units of cm
    
    if w_units == 'micron':
    
        w_cm = w / 10000
        
    elif w_units == 'cm':
        w_cm = w 

    else:
        return print('The function currently only works for w_units of micron or cm')
        
    sizes = SizeContainer()
    sizes.log_r = np.linspace(np.log10(amin_cm), np.log10(amax_cm), Nsize)
    sizes.N_log_r = sizes.r**(-2.5)   # because n(r) ∝ r^-3.5
    sizes.S_log_r = np.zeros_like(sizes.r)
    sizes.V_log_r = np.zeros_like(sizes.r)

    x_inp = 2 * np.pi * sizes.r / w_cm
    n_inp = n + np.zeros_like(x_inp)
    k_inp = k + np.zeros_like(x_inp)
    w_inp = w_cm + np.zeros_like(x_inp)
    theta = np.linspace(0, 180, 181)
    
    producer.setup(['Cext', 'Cabs', 'Csca', 'albedo', 'Z11', 'Z12', 'Z22', 'Z33', 'Z34', 'Z44', 
                                          'N11', 'N12', 'N22', 'N33', 'N34', 'N44'])

    mat = producer(x_inp, n_inp, k_inp, theta, w_inp,
                   xtype='vol',
                   outtype='MuellerMatrix')
    
    # Calculate some things
    # ------------------------------------------------------------------
    mat.r_vol = (x_inp * w_inp) / (2 * np.pi)
    # ------------------------------------------------------------------
    # kappa values using equation 20 from the paper
    rho_s = 1.675 # [grams / cm^3]
    mass = (4/3) * np.pi * rho_s * mat.r_vol**3
    # ------------------------------------------------------------------
    mat.k_ext = mat.Cext / mass
    mat.k_abs = mat.Cabs / mass
    mat.k_sca = mat.Csca / mass
    # ------------------------------------------------------------------
    
    
    # Look at distrbution

    mat_dist = {}
    for key in ['Z11', 'Z12', 'Z22', 'Z33', 'Z34', 'Z44']:
        mat_dist[key] = average_over_size(getattr(mat, key).T,
                                          sizes.N_log_r,
                                          sizes.d_log_r)
        

    Cext_avg = average_over_size(mat.Cext, sizes.N_log_r, sizes.d_log_r)
    Cabs_avg = average_over_size(mat.Cabs, sizes.N_log_r, sizes.d_log_r)
    Csca_avg = average_over_size(mat.Csca, sizes.N_log_r, sizes.d_log_r)
    mass_avg = average_over_size(mass, sizes.N_log_r, sizes.d_log_r)

    mat_dist['k_ext'] = Cext_avg / mass_avg
    mat_dist['k_abs'] = Cabs_avg / mass_avg
    mat_dist['k_sca'] = Csca_avg / mass_avg



    return theta, sizes, mat, mat_dist

# ------------------------------------------------------------ 
def RunGlitterinPowerLawDistAveraged(indicies, data_w_cm, data_n, data_k, max_rvol_grid_cm):

    # Make main storage dictionary
    results = {}


    # Loop over all of the wavelength index
    for i, idx in enumerate(indicies):

        # Get the index
        print(rf'Now working on: {idx}')

        # Get n, k, wavelength
        n = data_n[idx]
        k = data_k[idx]
        w_cm = data_w_cm[idx]
        print(rf'n = {n:.3f} and k = {k:.3f}')
        print(rf'lambda = {w_cm} cm')


        # Make empty arrays for this wavelength
        k_abs = np.zeros_like(max_rvol_grid_cm)
        k_sca = np.zeros_like(max_rvol_grid_cm)
        omega = np.zeros_like(max_rvol_grid_cm) 
        P = np.zeros_like(max_rvol_grid_cm) 
        P_omega = np.zeros_like(max_rvol_grid_cm) 

        # Loop over each rvol_max value for 
        for i, rvol_max_cm in enumerate(tqdm(max_rvol_grid_cm)):
            theta, sizes, mat, mat_dist = RunGlitterinPowerLawDist(rvol_max_cm, w_cm, 'cm', n, k)

            k_abs[i] = mat_dist['k_abs']
            k_sca[i] = mat_dist['k_sca']

            omega[i] = mat_dist['k_sca'] / (mat_dist['k_sca'] + mat_dist['k_abs'])

            idx_90 = np.argmin(np.abs(theta - 90))
            P[i] = - mat_dist['Z12'][idx_90] / mat_dist['Z11'][idx_90]

            P_omega[i] = P[i] * omega[i]

        # Save everything for this wavelength
        results[idx] = {
            'w_cm': w_cm,
            'n': n,
            'k': k,
            'rvol_max_cm': max_rvol_grid_cm,
            'rvol_max_micron': cm_to_micron(max_rvol_grid_cm),
    #         'rvol_max_um': max_rvol_grid_cm * 1e4,
            'k_abs': k_abs,
            'k_sca': k_sca,
    #         'omega': omega,
    #         'g': g,
            'omega': omega, 
            'P': P,
            'P_omega': P_omega,
        }
        
    return results
# ------------------------------------------------------------
#
#
#
#
# 
# I am updating this so we can say what POLF we want to use 
def find_sf_glitterin(results, AllBands, BandsInFit, df_POLF, POLF_index, print_results = True,
                      trouble=False,
                      trouble_index=0,
                      choose_trouble_index=False):
    
    
    # POLF Stuff
    # --------------------------------------------------------------------------
    # Name the columns of the POLF values
    POLF_columns = {
        'gaussian': 'POLF_Gaussian',
        'max Stokes I': 'POLF_maxStokesI',
        'POLI': 'POLF_maxPOLI',
        'mean': 'POLF_mean'
    }
    
        
    POLF_err_columns = {
    #'gaussian': 'POLF_err_Gaussian',
    'max Stokes I': 'POLF_err_maxStokesI',
    'POLI': 'POLF_err_maxPOLI',
    #'mean': 'POLF_err_mean'   # if you have this column
    }
    
    
    # POLF_index is the POLF we want to use to compare
    # Make sure we have that data 
    if POLF_index not in POLF_columns:
        raise ValueError(
            "POLF_index must be 'gaussian', 'poli', 'avg', 'mean', or 'max Stokes I'"
        )
    print(rf'The POLF value we are considering is: {POLF_index}')
    
    # Get all the POLF values together
    POLF_obs_all = df_POLF[POLF_columns[POLF_index]].values
    POLF_err_obs_all = df_POLF[POLF_err_columns[POLF_index]].values
    # --------------------------------------------------------------------------    
        
    
    print(f'\n The Bands we are looking at are: {AllBands}')
    print(rf"The Bands included in the fit are: {BandsInFit}")
 
    # Band names → indices (0, 1, 2, etc)
    band_to_index = {band: i for i, band in enumerate(AllBands)}
    fit_indices = [band_to_index[b] for b in BandsInFit]


    # Make an array for the rvol_max and make sure it is the same for all bands
    # --------------------------------------------------------------------------    
    rvol_max_cm = results[AllBands[0]]['rvol_max_cm']
#     rvol_max_micron = cm_to_micron(rvol_max_cm)
    
    ref = rvol_max_cm
    
    # Pick first band as reference
    ref = results[AllBands[0]]['rvol_max_cm']

    # Check all match
    all_same = all(
        np.allclose(results[band]['rvol_max_cm'], ref)
        for band in AllBands
    )

    print("\n All rvol_max_cm are the same:", all_same)
    # --------------------------------------------------------------------------    

    
    # Save the P omega data in the way that we need it to be saved
    # We want it saved like Pw = (N_rvol, N_bands), 
    # or one column per band, num r_vol_max long
    # --------------------------------------------------------------------------   
    N_rvol_max = len(rvol_max_cm)
    N_bands = len(AllBands)
    
    Pw = np.zeros((N_rvol_max, N_bands))
    
    for j, band in enumerate(AllBands):
        Pw[:, j] = results[band]['P_omega']
    # --------------------------------------------------------------------------   
    

    # Calculate the median and standard deviations
    # --------------------------------------------------------------------------   
    # Make arrays to store the median and standard deviation values 
    sf_medians = []
    sf_stds = []

    # Loop over each rvol_max grid point
    for i in range(len(rvol_max_cm )):
        Pw_model = Pw[i, :]  # all bands for current rvol_max

         # ONLY use selected bands for fit
        Pw_fit = Pw_model[fit_indices]
        POLF_fit = POLF_obs_all[fit_indices]

        # Avoid divide by zero
        valid = Pw_fit > 0

        sf_i = POLF_fit[valid] / Pw_fit[valid]

        sf_medians.append(np.median(sf_i))
        sf_stds.append(np.std(sf_i))

    sf_medians = np.array(sf_medians)
    sf_stds = np.array(sf_stds)
    
    best_idx = np.argmin(sf_stds)
    
    
    # --------------------------------------------------------------------------  
    
    
    # Trouble index
    # ------------------------------------------------------------

    if trouble:

        sorted_indices = np.argsort(sf_stds)

        print('\nTrouble case')
        print('Best few candidates:')

        for rank, idx in enumerate(sorted_indices[:4], start=1):
            print(
                f'{rank}: '
                f'a_max = {cm_to_micron(rvol_max_cm[idx]):.0f} micron, '
                f'sf = {sf_medians[idx]:.2f}, '
                f'sf_std = {sf_stds[idx]:.4f}'
            )

        if choose_trouble_index:

            if trouble_index >= len(sorted_indices):
                raise ValueError(
                    f"trouble_index={trouble_index} is too large. "
                    f"Only {len(sorted_indices)} valid candidates exist."
                )

            best_idx = sorted_indices[trouble_index]

            print(
                f'\nUsing trouble index {trouble_index}: '
                f'a_max = {cm_to_micron(rvol_max_cm[best_idx]):.0f} micron'
            )
            


    # Choose best a_max where scatter is minimized
    # ---------------------------------------------
    

    best_rvol_max_cm = rvol_max_cm[best_idx]
#     best_idx_arr = best_idx
    best_sf = sf_medians[best_idx]
    
    
    chi_sq = calculate_chi_squared_for_sf(
            Pw,
            best_idx,
            best_sf,
            POLF_obs_all,
            POLF_err_obs_all,
            AllBands,
        )
    
    POLF_obs = np.array([
        POLF_obs_all[band_to_index[b]] for b in AllBands
    ])
    # ---------------------------------------------
        
    if print_results:
        print("\nbest_idx:", best_idx)
        print("Best rvol_max_cm:", best_rvol_max_cm)
        print("Best SF:", best_sf)
        print("Best SF std:", sf_stds[best_idx])
    
    
    sf_results = {
        'rvol_max_cm': rvol_max_cm,
        'Pw': Pw,
        'best_rvol_max_cm': best_rvol_max_cm,
        'best_rvol_max_micron': cm_to_micron(best_rvol_max_cm),
        'POLF_obs': POLF_obs,
        'best_sf': best_sf,
        'best_idx_arr': best_idx,
        'sf_medians': sf_medians,
        'sf_stds': sf_stds,
        'chi_sq': chi_sq
    }

    return sf_results