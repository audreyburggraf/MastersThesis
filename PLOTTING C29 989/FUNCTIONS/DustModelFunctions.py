import numpy as np 
import miepython

def get_average_POLF_where_scattered(UniformRatios, POLF, tolerance=0.001):
    # Convert inputs to numpy arrays in case they aren't already
    UniformRatios = np.array(UniformRatios)
    POLF = np.array(POLF)
    
    # Create a mask where the ratio is within ±tolerance of 1
    mask = np.abs(UniformRatios - 1) <= tolerance
    
    # Use the mask to filter POLF_mJy values
    matching_POLF = POLF[mask]
    
    if matching_POLF.size == 0:
        return np.nan  # or raise an exception if preferred
    
    # Compute and return the average
    POLF_average = np.mean(matching_POLF)
    return POLF_average






# Function to convert micron to cm
# --------------------------------------------------------------------------------------
def micron_to_cm(value_micron):

    return value_micron / 1000
# --------------------------------------------------------------------------------------



# Function to find the mass of a spherical dust grain based on size and density
# --------------------------------------------------------------------------------------
def calculate_grain_mass(a_micron, density_g_cm3):
    
    # Convert a_micron to a_cm
    a_cm = micron_to_cm(a_micron)
    
    # Calculate the volume
    volume_cm3 = (4/3) * np.pi * a_cm**3
    
    # Find mass using p = m/V
    mass_grams = density_g_cm3 * volume_cm3 
    
    return mass_grams 
# --------------------------------------------------------------------------------------



# Size distribution of n(a) prop to a^-3.5 from Kataoka et al, 2015
# --------------------------------------------------------------------------------------
def grain_size_distribution(xmin, xmax, N):
    p = 3.5 
    
    # Normalization constant for f(x) = C * x^-p over [xmin, xmax]
    C = (1 - p) / (xmax**(1 - p) - xmin**(1 - p))
    
    # Define x values
    x = np.linspace(xmin, xmax, N)

    # Calculate PDF
    pdf = C * x**(-p)
    
    # Normalize weights to sum to 1 (optional but useful)
    weights = pdf / np.sum(pdf)
    
    return x, pdf, weights
# --------------------------------------------------------------------------------------
def sample_grain_sizes(a_min, a_max, N, p=3.5):
    # Power-law exponent related term
    inv_exp = 1 - p
    
    # Uniform random numbers from 0 to 1
    u = np.random.uniform(0, 1, N)
    
    # Inverse CDF transform to get grain sizes distributed as a^-p
    sizes = (u * (a_max**inv_exp - a_min**inv_exp) + a_min**inv_exp) ** (1 / inv_exp)
    
    return sizes
# --------------------------------------------------------------------------------------

# Function to the size parameter of x = 2*pi*a/lambda 
# --------------------------------------------------------------------------------------
def size_parameter(a_micron, lambda_micron):
    
    # Note: This can read in any units, just as long as the units of a and lambda are the same
    
    x_micron = (2 * np.pi * a_micron**2) / lambda_micron 
    
    return x_micron
# --------------------------------------------------------------------------------------





# --------------------------------------------------------------------------------------
import pandas as pd
m_data_folder_path = "/Users/audreyburggraf/Desktop/QUEEN'S/THESIS RESEARCH/PLOTTING C29 989/M FILES/"
# --------------------------------------------------------------------------------------
def get_complex_m(wavelengths_micron, material):
# Map material to filename
    if material == 'wi':
        file_name = "m_table_water_ice.txt"
    elif material == 'tr':
        file_name = "m_table_troilite.txt"
    elif material == 'as':
        file_name = "m_table_astronomical_silicate.txt"
    elif material == 'ro':
        file_name = "m_table_refractory_organics.txt"
    else:
        raise ValueError(f"Unknown material '{material}'. Choose from 'wi', 'tr', 'as', 'ro'.")

    # Handle the special case for astronomical silicate with extra header lines and extra columns
    if material == 'as':
        df = pd.read_csv(
            m_data_folder_path + file_name,
            delim_whitespace=True,
            header=None,
            skiprows=10,  # adjust if needed
            comment='#'
        )
        df.columns = ["wavelength (microns)", "Re_eps_minus_1", "Im_eps", "Re_m_minus_1", "Im_m"]
        # Compute actual refractive index from m-1 columns
        df["m_re"] = df["Re_m_minus_1"] + 1.0
        df["m_im"] = df["Im_m"]
        df_trim = df[["wavelength (microns)", "m_re", "m_im"]]
    else:
        # For other materials, simple three-column format assumed
        df = pd.read_csv(m_data_folder_path + file_name, delim_whitespace=True, header=None)
        df.columns = ["wavelength (microns)", "m_re", "m_im"]
        df_trim = df

    # Make sure wavelengths input is array-like for interpolation
    wavelengths = np.atleast_1d(wavelengths_micron)

    # Interpolate
    m_re = np.interp(wavelengths, df_trim["wavelength (microns)"], df_trim["m_re"])
    m_im = np.interp(wavelengths, df_trim["wavelength (microns)"], df_trim["m_im"])

    m_complex = m_re + 1j * m_im

    # Return scalar if input was scalar
    if np.isscalar(wavelengths_micron):
        return m_complex[0]

    return m_complex
# --------------------------------------------------------------------------------------
def find_complex_m_average(wavelengths_micron):
   
    
    m_wi = get_complex_m(wavelengths_micron, 'wi')
    m_tr = get_complex_m(wavelengths_micron, 'tr')
    m_as = get_complex_m(wavelengths_micron, 'as')
    m_ro = get_complex_m(wavelengths_micron, 'ro')
        
    # Mass fractions from DSHARP composition (sum should be 1.0)
    frac_wi = 0.20
    frac_as = 0.33
    frac_tr = 0.07
    frac_ro = 0.40

    # Weighted average of complex refractive indices
    m_complex_average = (frac_wi * m_wi +
                         frac_as * m_as +
                         frac_tr * m_tr +
                         frac_ro * m_ro)

    return m_complex_average
# --------------------------------------------------------------------------------------




# kappa_abs and kappa_sca
# --------------------------------------------------------------------------------------
# In Kataoka et al. 2015 page 2 they say:
# "We use Mie theory to calculate the absorption and scattering mass opacities"

# I will use the equations from Lin et al. 2023 on page 3:
# kappa_abs = Q_abs * pi * a^2 /m_g
# kappa_sca = Q_sca * pi * a^2 /m_g,
# where m_g is the mass of the grain
# --------------------------------------------------------------------------------------
# Information about miepython.efficiencies can be found at
# https://miepython.readthedocs.io/en/3.0.0/api/miepython.core.efficiencies.html
# m is the complex index of refraction of the sphere 
# x or d is the diameter of the sphere 
# lambda 0 is the wavelength in a vacuum (same units as lambda0)
# --------------------------------------------------------------------------------------
def calculate_kappa_wavelength(a_max_micron, lambda_micron_array, density_g_cm3, a_min_micron=0.001, N_dustgrains=100):
    
    # Grain size grid and PDF weights
    a_grid = np.logspace(np.log10(a_min_micron), np.log10(a_max_micron), N_dustgrains)
    _, pdf, _ = grain_size_distribution(a_min_micron, a_max_micron, N_dustgrains)
    weights = pdf / np.trapz(pdf, a_grid)  # normalize weights
    
    # Precompute mass and radius in cm for each grain size
    a_cm_array = micron_to_cm(a_grid)
    m_g_array = np.array([calculate_grain_mass(a, density_g_cm3) for a in a_grid])

    # Initialize outputs
    kappa_abs = np.zeros_like(lambda_micron_array, dtype=float)
    kappa_sca = np.zeros_like(lambda_micron_array, dtype=float)
    
    m_array = find_complex_m_average(lambda_micron_array)
    
    # Loop over wavelengths
    for i, (lambda_micron, m) in enumerate(zip(lambda_micron_array, m_array)):
        kappa_abs_sum = 0.0
        kappa_sca_sum = 0.0
        
        
        # Loop over grain sizes and weights
        for a_micron, a_cm, m_g, weight in zip(a_grid, a_cm_array, m_g_array, weights):
            x_micron = size_parameter(a_micron, lambda_micron)
            Q_ext, Q_sca, Q_back, g = miepython.efficiencies(m, x_micron, lambda_micron)
            Q_abs = Q_ext - Q_sca
            
            kappa_abs_sum += weight * (Q_abs * np.pi * a_cm**2) / m_g
            kappa_sca_sum += weight * (Q_sca * np.pi * a_cm**2) / m_g
        
        kappa_abs[i] = kappa_abs_sum
        kappa_sca[i] = kappa_sca_sum

    return kappa_abs, kappa_sca
# --------------------------------------------------------------------------------------
def calculate_kappa_grainsize(a_max_micron, lambda_micron, density_g_cm3, a_min_micron=0.001, N_dustgrains=100):
    
    m = find_complex_m_average(lambda_micron)
    
    
    kappa_abs = np.zeros_like(a_max_micron_array, dtype=float)
    kappa_sca = np.zeros_like(a_max_micron_array, dtype=float)

   
    for i, a_max_micron in enumerate(a_max_micron_array):
       
        a_grid = np.logspace(np.log10(a_min_micron), np.log10(a_max_micron), N_dustgrains)
        _, pdf, _ = grain_size_distribution(a_min_micron, a_max_micron, N_dustgrains)
        weights = pdf / np.trapz(pdf, a_grid)

        kappa_abs_sum = 0.0
        kappa_sca_sum = 0.0

        for a_micron, weight in zip(a_grid, weights):
            
            x_micron = size_parameter(a_micron, lambda_micron)
           
            Q_ext, Q_sca, Q_back, g = miepython.efficiencies(m, x_micron, lambda_micron)
            Q_abs = Q_ext - Q_sca

            m_g = calculate_grain_mass(a_micron, density_g_cm3)
            a_cm = micron_to_cm(a_micron)

            kappa_abs_sum += weight * (Q_abs * np.pi * a_cm**2) / m_g # units of cm^2 per gram of dist grains
            kappa_sca_sum += weight * (Q_sca * np.pi * a_cm**2) / m_g # units of cm^2 per gram of dist grains

        kappa_abs[i] = kappa_abs_sum
        kappa_sca[i] = kappa_sca_sum

    return kappa_abs, kappa_sca
# --------------------------------------------------------------------------------------





# This function is taken from Kataoka et al. 2015, page 2
# --------------------------------------------------------------------------------------
def calculate_omega(kappa_abs, kappa_sca):
    
    omega = kappa_sca / (kappa_abs + kappa_sca)
    
    return omega
# --------------------------------------------------------------------------------------






# Calcluate degree of polarization
# --------------------------------------------------------------------------------------
# From Kataoka et al, 2015 they find Z11 and Z12 using a scattering matrix
# --------------------------------------------------------------------------------------
def calculate_P_scatt_angle(lambda_micron, scatt_angles_deg_array, a_max_micron, a_min_micron=0.001, N=100):
    # Grain size array (microns)
    grain_sizes_micron = np.logspace(np.log10(a_min_micron), np.log10(a_max_micron), N)

    # Compute PDF from grain size distribution and normalize
    grainsize_pdf = sample_grain_sizes(a_min_micron, a_max_micron, N, p=3.5)

    # Empty arrays for Z11 and Z12
    Z11_total = np.zeros_like(scatt_angles_deg_array, dtype=float)
    Z12_total = np.zeros_like(scatt_angles_deg_array, dtype=float)

    # Refractive index from Lin et al. (2023)
    m = find_complex_m_average(lambda_micron)

    for j, a in enumerate(grain_sizes_micron):
        x = size_parameter(a, lambda_micron)

        for i, theta_deg in enumerate(scatt_angles_deg_array):
            theta_rad = np.deg2rad(theta_deg)
            Z = miepython.phase_matrix(m, x, theta_rad)

            Z11_total[i] += Z[0][0] 
            Z12_total[i] += Z[0][1] 

    # Normalize (optional, for visualization)
    Z11_total /= np.sum(Z11_total)
    Z12_total /= np.sum(Z11_total)

    pol_deg = -Z12_total / Z11_total

    return pol_deg
# --------------------------------------------------------------------------------------
def calculate_P_grainsize(lambda_micron, scatt_angle_deg, a_max_micron_array, a_min_micron=0.001, N=100):
    m = find_complex_m_average(lambda_micron)

    # Convert angle to radians
    theta_rad = np.deg2rad(scatt_angle_deg)

    # Output array
    pol_deg_array = np.zeros_like(a_max_micron_array, dtype=float)

    for i, a_max in enumerate(a_max_micron_array):
        grainsize_pdf = sample_grain_sizes(a_min_micron, a_max_micron, N, p=3.5)


        Z11_total = 0.0
        Z12_total = 0.0

        for a, weight in zip(grain_sizes, weights):
            x = size_parameter(a, lambda_micron)
            
            Z = miepython.phase_matrix(m, x, theta_rad)
            
            Z11_total += Z[0][0] 
            Z12_total += Z[0][1] 
        pol_deg_array[i] = -Z12_total / Z11_total if Z11_total != 0 else 0.0

    return pol_deg_array
# --------------------------------------------------------------------------------------