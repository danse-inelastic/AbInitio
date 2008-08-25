__doc__ = """This modules provides functions to calculate electronic heat capacities and entropy from an electron DOS."""

import numpy as np

kB = 8.6173857e-5 # in eV/K

def Lorentzian(x, center=0.0, gamma=1.0):
    """Lorentzian function:
    L(x) = 1/PI * gamma/2 / ( (x-x0)^2 + (gamma/2)^2)
    of center x0 and FWHM gamma.
    The Lorentzian is normalized to a unit integral."""
    return gamma/2.0 / (np.pi * ((x-center)**2.0 + (gamma/2.0)**2.0))

def convolvePSF(x, y, kernel):
    """This function calculates the convolution product of the function y(x),
    passed as two 1D arrays of same dimensions (x and y),
    with the kernel function kernel (an actual python function
    returning scalars, and which should be defined for all x values).
    The returned results are normalized to the integral of the input array."""
    if (len(x) != len(y)): raise ValueError, "x and y should have same dimension"
    size = len(x)
    res = np.zeros(size, 'd')
    dx = x[1]-x[0]
    #norm = np.sum(y)*dx
    for n in range(size):
        res[n] = 0
        for j in range(size):
            res[n] = res[n] + y[j]*kernel(x[n]-x[j])
    #res /= norm
    res *= dx
    return res

def convolveSymArr(x, y, ker):
    """convolves the y array representing a function y(x),
    with a kernel function k(x), represented by an array 'kernel'.
    The kernel function is supposed to be even around zero,
    and the array kernel represents only the positive side of the kernel."""
    dx = x[1]-x[0]
    n=len(x)
    kk = np.zeros((n, n), 'd')
    for j in range(n):
        kk[j:n,j] = ker[0:n-j]
    diag = np.diag(kk)
    kk = kk + np.transpose(kk)
    for i in range(n):
        kk[i,i] = diag[i]
    return np.dot(kk, y)*dx
    


def integrate(x,y):
    """Returns the integral of the function y(x),
    defined by the arrays x and y."""
    if (len(x) != len(y)): raise ValueError, "x and y should have same dimension"
    size = len(x)
    res = np.zeros(size, 'd')
    for n in range(size-1):
        res[n+1] = res[n] + (x[n+1]-x[n])*y[n]
    return res

def integral(x,y):
    """return the numerical integral of y(x)"""
    if (len(x) != len(y)): raise ValueError, "x and y should have same dimension"
    size = len(x)
    res = 0
    for n in range(size-1):
        res = res + (x[n+1]-x[n])*y[n]
    return res
    

def findFillEnergy(e, intd, Nel):
    """this function solves for the chemical potential.
    It takes the integrated electron dos defined by two 1D arrays
    for the energy and the integrated dos, and the number of electrons
    for which to find the chemical potential."""
    close_bin = np.argmin( (intd - Nel)**2 )
    # the integrated dos always has positive slope,
    # since the derivative is the dos...
    # so we ahve to distinguish based on wether
    # the minum found is larger or smaller than Nel:
    if (intd[close_bin]) < Nel:
        # the zero is between (approx_zero_bin) and (approx_zero_bin + 1)
        slope = (intd[close_bin+1]-intd[close_bin]) / (e[close_bin+1]-e[close_bin])
        mu = e[close_bin] - (intd[close_bin]-Nel) / slope
    else:
        # the zero is between (approx_zero_bin-1) and (approx_zero_bin)
        slope = (intd[close_bin]-intd[close_bin-1]) / (e[close_bin]-e[close_bin-1])
        mu = e[close_bin] - (intd[close_bin]-Nel) / slope
    return mu

def intdosfermi(e, d, mu, T):
    """calculates the integral of the electron DOS times the Fermi distribution,
    for a given chemical potential mu (eV) and a temperature T (K).
    The electron DOS is d in states/eV and energy is e (eV)."""
    fet = np.array([FermiMu(x, mu, T) for x in e ])
    weightedDos = np.array(d * fet)
    return integral(e, weightedDos)

def findmu(e,d, mus, T, Nel):
    """Calculate the chemical potential at temperature T (K),
    for an electron dOS d(e) and a number of electrons Nel,
    using a discretized array of energy-shifts 'mus'."""
    nels = np.array([intdosfermi(e, d, mu, T) for mu in mus])
    return findFillEnergy(mus, nels, Nel)
    

def findNEF(e, d, ef):
    """This finds the (interpolated) density at the Fermi level,
    from an array 'e' for the energies,
    and an array 'd' for the electron DOS.
    'ef' is the Fermi energy."""
    close_bin = np.argmin( (e-ef)**2 )
    # we ahve to distinguish based on wether
    # the mininum found is larger or smaller than 0:
    if (e[close_bin] <= ef):
        if ((e[close_bin]-ef) * (e[close_bin+1]-ef) > 0.0):
            raise ValueError, 'E_Fermi is not in the energy range'
        # Efermi is between (approx_zero_bin) and (approx_zero_bin + 1)
        slope = (d[close_bin+1]-d[close_bin]) / (e[close_bin+1]-e[close_bin])
        dfermi = d[close_bin] + slope * (ef - e[close_bin])
    elif (e[close_bin] > ef ):
        if ((e[close_bin]-ef) * (e[close_bin+1]-ef) > 0.0):
            raise ValueError, 'E_Fermi is not in the energy range'
        # the zero is between (approx_zero_bin-1) and (approx_zero_bin)
        slope = (d[close_bin] - d[close_bin-1]) / (e[close_bin]-e[close_bin-1])
        dfermi = d[close_bin] - slope * (e[close_bin] - ef)
    else:
        raise ValueError, 'E_Fermi could not be found.'
    return dfermi
    

def gamma(lam, T):
    """calculates the smearing width for electron-phonon
    lorentzian broadening of electronic levels:
    gamma = PI * lambda * kB * T.
    T is expected in Kelvin and gamma returned in eV."""
    return 8.617343e-5 * np.pi * lam * T

def mutemp(e, d, Nel, lam, Ts):
    """calculates the chemical potential as function of temperature (array Ts) in K,
    for an electron DOS given by d (states/ev/atom) at energies e (in eV)."""
    muT = np.zeros(len(Ts))
    iT = 0
    for T in Ts:
        gamT = gamma(lam, T)
        #print "%s K: gamma = %s \n" % (T, gamT)
        de = e[1]-e[0]
        kerT = np.array([Lorentzian(x*de, 0.0, gamT) for x in range(len(e))])
        #pl.clf()
        #pl.plot(kerT)
        #blah = raw_input()
        depT = convolveSymArr(e, d, kerT)
        intdepT = integrate(e, depT)
        #pl.clf()
        #pl.plot(e,d,e,depT, e, intdepT)
        #blah = raw_input()
        muT[iT] = findmu(e, intdepT, Nel)
        #print "mu = %s \n" % muT[iT]
        iT += 1
    return muT

def derivate(x, y):
    """Calculates the numerical derivative of y(x),
    with x and y given by 1d arrays."""
    n = len(x)
    res = np.zeros(n-1, 'd')
    for i in range(n-1):
        res[i] = (y[i+1]-y[i]) / (x[i+1]-x[i])
    return res
    

def antiderivate(x, y):
    """Returns the antiderivative of the function y(x),
    defined by the arrays x and y."""
    if (len(x) != len(y)): raise ValueError, "x and y should have same dimension"
    size = len(x)
    res = np.zeros(size, 'd')
    for n in range(size-1):
        res[n+1] = res[n] + (x[n+1]-x[n])*y[n]
    return res


def Fermi(e, T):
    """calculate the Fermi distribution for energy e in eV (wrt to chemical potential mu),
    at temperature T (in K)."""
    return 1.0 / (np.exp(e/(8.6173857e-5*T)) + 1.0)

def FermiMu(e, mu, T):
    """calculate the Fermi distribution for energy e in eV,
    with a chemical potential mu in eV,
    and a temperature T in K."""
    return 1.0 / (np.exp((e-mu)/(8.6173857e-5*T)) + 1.0)

def derivFermi(e,T):
    """Returns the opposite of the derivative of the Fermi distribution function:
    -df/dE ."""
    fet = Fermi(e,T)
    return fet * (1-fet) / (8.6173857e-5 * T)

def xlnx(x):
    """returns x*Ln(x), and its continuation of x*Ln(x) in zero."""
    if (x>0):
        return x* np.log(x)
    elif (x==0):
        return 0.0
    else:
        return None
    

def thermalFactor(e, T):
    """calculates the thermal factor for electronic entropy."""
    return -(xlnx(Fermi(e, T)) + xlnx(1.0 - Fermi(e,T)))

def Selbare(e, d, T):
    """Calculates the bare electronic entropy from a DOS."""
    thmwgt = np.array([thermalFactor(x, T) for x in e])
    integrand = thmwgt * np.array(d)
    antideriv = antiderivate(e, integrand)
    return antideriv[-1]


def Selph(e, d, Nel, lam, T):
    """Returns the electronic entropy for the density of states d,
    broadened by a Lorentzian of width Gamma(T),
    representing the el-ph coupling effect.
    T is the temperature in K, e the energy array in eV."""
    #calculate the kernel
    de = e[1]-e[0]
    kerT = np.array([Lorentzian(i*de, 0.0, gamma(lam, T)) for i in range(len(e))])
    #calculate the smeared electron DOS
    depT = convolveSymArr(e, d, kerT)
    intdepT = antiderivate(e, depT)
    #calculate the chemical potential mu(T)
    muT = findmu(e, intdepT, Nel)
    #calculate the entropy
    Sel = Selbare( (e-muT), depT, T)
    return Sel


    

def Celbare(e, d, T):
    """Returns the bare electron heat capacity."""
    pass
    return
