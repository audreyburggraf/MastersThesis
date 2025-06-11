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
# I got ChatGPT to write this
# --------------------------------------------------------------------------------------
def size_distribution_weights(grain_sizes):
    """
    Returns normalized weights for a power-law n(a) ∝ a^-3.5
    """
    weights = grain_sizes**-3.5
    return weights / np.sum(weights)  # Normalize so the weights sum to 1

# --------------------------------------------------------------------------------------



# Function to the size parameter of x = 2*pi*a/lambda 
# --------------------------------------------------------------------------------------
def size_parameter(a_micron, lambda_micron):
    
    # Convert a_micron to a_cm
    a_cm = micron_to_cm(a_micron)
    
    # Convert lambda_micron to cm
    lambda_cm = micron_to_cm(lambda_micron)
    
    x_cm = (2 * np.pi * a_cm) / lambda_cm
    
    return x_cm
# --------------------------------------------------------------------------------------


# kappa_abs and kappa_sca
# --------------------------------------------------------------------------------------
# In Kataoka et al. 2015 page 2 they say:
# "We use Mie theory to calculate the absorption and scattering mass opacities"

# I will use the equations from Lin et al. 2023 on page 3:
# kappa_abs = Q_abs * pi * a^2 /m_g
# kappa_sca = Q_sca * pi * a^2 /m_g,
# where m_g is the mass of the grain

# I used ChatGPT to help break the function up by if wavelength od grain size is changing,
# and to add a docstring
# --------------------------------------------------------------------------------------
def calculate_kappa(lambda_micron, a_max_micron, density_g_cm3, func_of='wavelength'):
    """
    Calculate absorption and scattering mass opacities (kappa_abs, kappa_sca)
    either as a function of wavelength (at fixed grain size distribution)
    or as a function of grain size (at fixed wavelength), including size distribution.

    Parameters:
    - lambda_micron: scalar or array of wavelengths [micron]
    - a_max_micron: maximum grain size [micron] for distribution
    - density_g_cm3: grain density [g/cm^3]
    - func_of: 'wavelength' or 'grainsize'

    Returns:
    - kappa_abs: absorption opacity [cm^2/g]
    - kappa_sca: scattering opacity [cm^2/g]
    """

    # Refractive index from Lin et al, 2023
    m = complex(1.65, 1e-5)

    a_min = 0.01  # micron
    N = 100
    grain_sizes = np.logspace(np.log10(a_min), np.log10(a_max_micron), N)
    weights = size_distribution_weights(grain_sizes)

    if func_of == 'wavelength':
        lambda_micron = np.atleast_1d(lambda_micron)
        kappa_abs = np.zeros_like(lambda_micron, dtype=float)
        kappa_sca = np.zeros_like(lambda_micron, dtype=float)

        for i, lam in enumerate(lambda_micron):
            kappa_abs_sum = 0.0
            kappa_sca_sum = 0.0

            for a_micron, weight in zip(grain_sizes, weights):
                x_cm = size_parameter(a_micron, lam)
                Q_ext, Q_sca, Q_back, g = miepython.efficiencies(m, x_cm)
                Q_abs = Q_ext - Q_sca

                m_g = calculate_grain_mass(a_micron, density_g_cm3)
                a_cm = micron_to_cm(a_micron)

                kappa_abs_sum += weight * (Q_abs * np.pi * a_cm**2) / m_g
                kappa_sca_sum += weight * (Q_sca * np.pi * a_cm**2) / m_g

            kappa_abs[i] = kappa_abs_sum
            kappa_sca[i] = kappa_sca_sum

    elif func_of == 'grainsize':
        lam = lambda_micron if np.isscalar(lambda_micron) else lambda_micron[0]
        kappa_abs = np.zeros_like(grain_sizes, dtype=float)
        kappa_sca = np.zeros_like(grain_sizes, dtype=float)

        for i, a_micron in enumerate(grain_sizes):
            x_cm = size_parameter(a_micron, lam)
            Q_ext, Q_sca, Q_back, g = miepython.efficiencies(m, x_cm)
            Q_abs = Q_ext - Q_sca

            m_g = calculate_grain_mass(a_micron, density_g_cm3)
            a_cm = micron_to_cm(a_micron)

            kappa_abs[i] = (Q_abs * np.pi * a_cm**2) / m_g
            kappa_sca[i] = (Q_sca * np.pi * a_cm**2) / m_g

    else:
        raise ValueError("func_of must be 'wavelength' or 'grainsize'")

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
def calculate_pol_deg(lambda_micron, scatt_angles_deg, a_max_micron, a_min_micron=0.001):
    
    # Grain size array (microns)
    grain_sizes_micron = np.logspace(np.log10(a_min_micron), np.log10(a_max_micron), 100)
    weights = size_distribution_weights(grain_sizes_micron)

    # Empty arrays for Z11 and Z12
    Z11_total = np.zeros_like(scatt_angles_deg, dtype=float)
    Z12_total = np.zeros_like(scatt_angles_deg, dtype=float)

    # Refractive index from Lin et al. (2023)
    m = complex(1.65, 1e-5)

    for j, a in enumerate(grain_sizes_micron):
        weight = weights[j]
        x = size_parameter(a, lambda_micron)

        for i, theta_deg in enumerate(scatt_angles_deg):
            theta_rad = np.deg2rad(theta_deg)
            Z = miepython.phase_matrix(m, x, theta_rad)

            Z11_total[i] += Z[0][0] * weight
            Z12_total[i] += Z[0][1] * weight

    # Normalize (optional, for meaningful plotting)
    Z11_total /= np.sum(Z11_total)
    Z12_total /= np.sum(Z11_total)

    pol_deg = -Z12_total / Z11_total

    return pol_deg

# --------------------------------------------------------------------------------------