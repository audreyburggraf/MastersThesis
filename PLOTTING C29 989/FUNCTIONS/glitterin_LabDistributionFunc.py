#Note: These are not my functions and these are from Daniel
import numpy as np




class SizeContainer:
    """
    I will use r to be the radius of the particles. r here is in units of micron

    The number of particles within the range r to r + dr is defined as 
        n(r) dr
    In this case, n(r) has units of counts per micron. 

    Let y := log(r)
    The experiment measures the number of particles within y to y + dy
    That value is N_log_r which has units of counts

    So these two are different quantities. Be careful
    """

    QUANTITY = ['log_r', 'N_log_r', 'S_log_r', 'V_log_r']

    def __init__(self, name=None, log_r=None, N_log_r=None, S_log_r=None, V_log_r=None):
        self.name = name
        self.log_r = log_r
        self.N_log_r = N_log_r
        self.S_log_r = S_log_r
        self.V_log_r = V_log_r

    def x(self, w_micron):
        """
        the size parameter if we know the wavelength

        Note that r in this object is in micron
        """
        return 2 * np.pi * self.r / w_micron

    def resample(self, n=2):
        """
        for purposes of size averaging opacities, it may be useful to increase the sampling, especially for Mie theory with fine oscillations

        Returns
        -------
            a SizeContainer object
        """
        new_n = int(np.ceil(n * len(self.log_r)))

        new_grid = np.arange(new_n) / (new_n-1) * (len(self.log_r)-1)
        log_r = np.interp(new_grid, np.arange(len(self.log_r)), self.log_r)

        N_log_r = np.interp(log_r, self.log_r, self.N_log_r)

        S_log_r = np.interp(log_r, self.log_r, self.S_log_r)

        V_log_r = np.interp(log_r, self.log_r, self.V_log_r)

        out = SizeContainer(name=self.name,
            log_r=log_r,
            N_log_r=N_log_r,
            S_log_r=S_log_r,
            V_log_r=V_log_r
        )

        return out

    @property
    def r(self):
        """
        instead of the logarithmic quantity, this is the actual size in microns
        """
        return 10**self.log_r

    @property
    def d_log_r(self):
        """
        calculate the bin size
        We only know log_r defined on cells, so I need to do a bit of extrapolation
        """
        walls = np.zeros([len(self.log_r)+1])

        walls[1:-1] = 0.5 * (self.log_r[1:] + self.log_r[:-1])

        # extrapolate to first point
        walls[0] = self.log_r[0] - (walls[1] - self.log_r[0])

        # extrapolate to last point
        walls[-1] = self.log_r[-1] + (self.log_r[-1] - walls[-2])

        # now I can get the cell sizes
        d_log_r = np.diff(walls)

        return d_log_r

    @property
    def n_r(self):
        """
        This is the quantity related to the size distribution 

        The log r of the database is log base 10

        n(r) dr = N(log r) d log r
            = N(log r) / r / ln(10) dr

        in units of counts per unit length
        """
        return self.N_log_r / self.r / np.log(10.)

    @property
    def n_pwl(self):
        """
        The powerlow of the size density distribution
        """
        return np.log(self.n_r[1:] / self.n_r[:-1]) / np.log(self.r[1:] / self.r[:-1])

    def __repr__(self):
        txt = 'SizeContainer from gals_database.py'

        return txt

    def __str__(self):
        txt = 'SizeContainer from gals_database.py'
        if self.name is not None:
            txt += 'name: {}\n'.format(self.name)

        try:
            txt += 'log_r range: [{}, ]\n'.format(self.log_r[0], self.log_r[-1])
        except:
            pass

        return txt
    
    
    
    
def read_gals_size_table(fname, skiprows=0, out_type='table'):
    """
    in most cases, we only have
    log r, N(log r), S(log r), V(log r)

    This will be the default file reader
    """
    # read the entire file
    with open(fname, 'r') as f:
        tab = np.loadtxt(fname, skiprows=skiprows)

    if out_type == 'table':
        return tab

    elif out_type == 'container':
        sizes = SizeContainer()
        column_to_attribute = {0:'log_r', 1:'N_log_r', 2:'S_log_r', 3:'V_log_r'}
        for icolumn in column_to_attribute:
            ikey = column_to_attribute[icolumn]
            setattr(sizes, ikey, tab[:,icolumn])

        return sizes
    
    
    
    
    
def average_over_size(C, N_log_r, d_log_r):
    """
    average a cross-sectional quantity over the size distribution

    Let the number of particles in the range [r, r+dr] be n(r)dr
    
    The contribution to the intensity using Cabs as an example: 
        d I_nu = int n(r) Cabs dr B_nu
    The averaged Cabs is the representative value for the number of particles that will produce the same level of d I_nu which means
        avg(Cabs) int n(r) dr = int n(r) Cabs dr

    This is the definition of avg(Cabs). It's equivalent for Cext and Csca. 

    The scattering matrix quantities are in cross-section per steradian. 
        d I_nu(Omega) = int n(r) Z11(theta) dr I_0 dOemga
    This won't affect how we integrate. 

    The experiments usually give n(log r) at locations in log r instead of the n(r) in r. 
    Instead of converting to n(r), we should use this logarithmic sampling which helps with integration anyways. 

    The database has a detailed explanation of what their measurements of N_log_r means. 
    https://scattering.iaa.es/size

    The total number of particles from r1 to r2 is: 

        int_{r1}^{r2} n(r) dr = int_{log r1}^{log r2} N(log r) d log r

    Their N_log_r is precisely N(log r) here. Apparently, their n(r) is normalized such that the integral is 1. 

    The average quantity we need to assess is 
        avg(Cabs) = int_{r1}^{r2} n(r) Cabs dr 
            = int_{log r1}^{log r2} N(log r) Cabs d log r

    Numerically, we can express this as 
        avg(Cabs) = sum_i N(log_r)_i Cabs_i d log r_i

    Parameters
    ----------
    C : ndarray
        The first dimension is the size dimension to average over

    Returns
    -------

    """
    Ntot = np.sum(N_log_r * d_log_r)

    dim = C.shape
    if len(dim) == 1:
        # we only have a 1d ndarray
        Ctot = np.sum(C * N_log_r * d_log_r)
    else:
        expand = tuple(range(1, len(dim)))
        Ctot = np.sum(C * np.expand_dims(N_log_r * d_log_r, axis=expand), axis=0)

    avg = Ctot / Ntot

    return avg
# ------------------------------------------------------------
# Import the necessary packages
# ------------------------------------------------------------
import glitterin
# ------------------------------------------------------------


class SizeContainer:
    """
    I will use r to be the radius of the particles. r here is in units of micron

    The number of particles within the range r to r + dr is defined as 
        n(r) dr
    In this case, n(r) has units of counts per micron. 

    Let y := log(r)
    The experiment measures the number of particles within y to y + dy
    That value is N_log_r which has units of counts

    So these two are different quantities. Be careful
    """

    QUANTITY = ['log_r', 'N_log_r', 'S_log_r', 'V_log_r']

    def __init__(self, name=None, log_r=None, N_log_r=None, S_log_r=None, V_log_r=None):
        self.name = name
        self.log_r = log_r
        self.N_log_r = N_log_r
        self.S_log_r = S_log_r
        self.V_log_r = V_log_r

    def x(self, w_micron):
        """
        the size parameter if we know the wavelength

        Note that r in this object is in micron
        """
        return 2 * np.pi * self.r / w_micron

    def resample(self, n=2):
        """
        for purposes of size averaging opacities, it may be useful to increase the sampling, especially for Mie theory with fine oscillations

        Returns
        -------
            a SizeContainer object
        """
        new_n = int(np.ceil(n * len(self.log_r)))

        new_grid = np.arange(new_n) / (new_n-1) * (len(self.log_r)-1)
        log_r = np.interp(new_grid, np.arange(len(self.log_r)), self.log_r)

        N_log_r = np.interp(log_r, self.log_r, self.N_log_r)

        S_log_r = np.interp(log_r, self.log_r, self.S_log_r)

        V_log_r = np.interp(log_r, self.log_r, self.V_log_r)

        out = SizeContainer(name=self.name,
            log_r=log_r,
            N_log_r=N_log_r,
            S_log_r=S_log_r,
            V_log_r=V_log_r
        )

        return out

    @property
    def r(self):
        """
        instead of the logarithmic quantity, this is the actual size in microns
        """
        return 10**self.log_r

    @property
    def d_log_r(self):
        """
        calculate the bin size
        We only know log_r defined on cells, so I need to do a bit of extrapolation
        """
        walls = np.zeros([len(self.log_r)+1])

        walls[1:-1] = 0.5 * (self.log_r[1:] + self.log_r[:-1])

        # extrapolate to first point
        walls[0] = self.log_r[0] - (walls[1] - self.log_r[0])

        # extrapolate to last point
        walls[-1] = self.log_r[-1] + (self.log_r[-1] - walls[-2])

        # now I can get the cell sizes
        d_log_r = np.diff(walls)

        return d_log_r

    @property
    def n_r(self):
        """
        This is the quantity related to the size distribution 

        The log r of the database is log base 10

        n(r) dr = N(log r) d log r
            = N(log r) / r / ln(10) dr

        in units of counts per unit length
        """
        return self.N_log_r / self.r / np.log(10.)

    @property
    def n_pwl(self):
        """
        The powerlow of the size density distribution
        """
        return np.log(self.n_r[1:] / self.n_r[:-1]) / np.log(self.r[1:] / self.r[:-1])

    def __repr__(self):
        txt = 'SizeContainer from gals_database.py'

        return txt

    def __str__(self):
        txt = 'SizeContainer from gals_database.py'
        if self.name is not None:
            txt += 'name: {}\n'.format(self.name)

        try:
            txt += 'log_r range: [{}, ]\n'.format(self.log_r[0], self.log_r[-1])
        except:
            pass

        return txt
    
    
    
    
def read_gals_size_table(fname, skiprows=0, out_type='table'):
    """
    in most cases, we only have
    log r, N(log r), S(log r), V(log r)

    This will be the default file reader
    """
    # read the entire file
    with open(fname, 'r') as f:
        tab = np.loadtxt(fname, skiprows=skiprows)

    if out_type == 'table':
        return tab

    elif out_type == 'container':
        sizes = SizeContainer()
        column_to_attribute = {0:'log_r', 1:'N_log_r', 2:'S_log_r', 3:'V_log_r'}
        for icolumn in column_to_attribute:
            ikey = column_to_attribute[icolumn]
            setattr(sizes, ikey, tab[:,icolumn])

        return sizes
    
    
    
    
    
def average_over_size(C, N_log_r, d_log_r):
    """
    average a cross-sectional quantity over the size distribution

    Let the number of particles in the range [r, r+dr] be n(r)dr
    
    The contribution to the intensity using Cabs as an example: 
        d I_nu = int n(r) Cabs dr B_nu
    The averaged Cabs is the representative value for the number of particles that will produce the same level of d I_nu which means
        avg(Cabs) int n(r) dr = int n(r) Cabs dr

    This is the definition of avg(Cabs). It's equivalent for Cext and Csca. 

    The scattering matrix quantities are in cross-section per steradian. 
        d I_nu(Omega) = int n(r) Z11(theta) dr I_0 dOemga
    This won't affect how we integrate. 

    The experiments usually give n(log r) at locations in log r instead of the n(r) in r. 
    Instead of converting to n(r), we should use this logarithmic sampling which helps with integration anyways. 

    The database has a detailed explanation of what their measurements of N_log_r means. 
    https://scattering.iaa.es/size

    The total number of particles from r1 to r2 is: 

        int_{r1}^{r2} n(r) dr = int_{log r1}^{log r2} N(log r) d log r

    Their N_log_r is precisely N(log r) here. Apparently, their n(r) is normalized such that the integral is 1. 

    The average quantity we need to assess is 
        avg(Cabs) = int_{r1}^{r2} n(r) Cabs dr 
            = int_{log r1}^{log r2} N(log r) Cabs d log r

    Numerically, we can express this as 
        avg(Cabs) = sum_i N(log_r)_i Cabs_i d log r_i

    Parameters
    ----------
    C : ndarray
        The first dimension is the size dimension to average over

    Returns
    -------

    """
    Ntot = np.sum(N_log_r * d_log_r)

    dim = C.shape
    if len(dim) == 1:
        # we only have a 1d ndarray
        Ctot = np.sum(C * N_log_r * d_log_r)
    else:
        expand = tuple(range(1, len(dim)))
        Ctot = np.sum(C * np.expand_dims(N_log_r * d_log_r, axis=expand), axis=0)

    avg = Ctot / Ntot

    return avg
