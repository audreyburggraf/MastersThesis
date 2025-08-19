import numpy as np 

# Modified Bessel function of the first kind of real order 
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.iv.html
from scipy.special import iv  

# Local minimization of the scalar function of one variable
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize_scalar.html
from scipy.optimize import minimize_scalar # this is the new tnmin

def bessel_for_debias(P_obs, sigma_PI):

    def neg_pdf(P_I):
        x = (P_obs * P_I) / sigma_PI**2
        
        # (v, z): v is order, z is argument 
        I_0_x = iv(0, x)
        
        PDF = (P_obs / sigma_PI**2) * I_0_x * np.exp(-(P_obs**2 + P_I**2) / (2 * sigma_PI**2))
        
        return - PDF

    # (function, bounds, method)
    result = minimize_scalar(neg_pdf, bounds=(0, P_obs), method='bounded')
    
    return result.x if result.success else 0.0 # return 0 if the minimization was unsuccessful






from tqdm import tqdm

def debias_func(POLI_map, POLI_err_map, StokesQ_map, StokesU_map, band):
    ny, nx = POLI_map.shape
    POLI_debiased = np.zeros_like(POLI_map)

    for y in tqdm(range(ny), desc=f"Debiasing {band}"):
        for x in range(nx):
            sigma_PI = POLI_err_map[y, x]
            P_obs = POLI_map[y, x]
            StokesQ = StokesQ_map[y, x]
            StokesU = StokesU_map[y, x]
            S_N = P_obs / sigma_PI

            if S_N > 9:
                debiased = np.sqrt(StokesQ**2 + StokesU**2 - sigma_PI**2)
            else:
                debiased = bessel_for_debias(P_obs, sigma_PI)
               

            POLI_debiased[y, x] = debiased

    return POLI_debiased
