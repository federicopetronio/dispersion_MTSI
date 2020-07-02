

from numerical_z.plasma import PlasmaParameters
from numerical_z.numerical_z import NumericZ

from numerical_z.evdfs import maxw_f
from plasmapy.mathematics import plasma_dispersion_func

import numpy as np
import scipy as sc

"Normalized permittivity"


def eps_IAW(omg, kx, ky, kz, prt=PlasmaParameters(),
            Ze=plasma_dispersion_func,
            Zi=plasma_dispersion_func):
    def Zep(x):
        return -2 * (1 + x * Ze(x))

    def Zip(x):
        return -2 * (1 + x * Zi(x))

    kp = np.sqrt(ky ** 2 + kz ** 2)
    eBeta = kp ** 2 * prt.ope ** 2 / prt.oce ** 2
    k2 = kx ** 2 + kp ** 2

    if k2 == 0:
        raise RuntimeError("The wave vector is Zero !!")

    iChi = (omg) / (ky * prt.vti)
    eChi = (omg - ky * prt.vde) / (ky * prt.vte)

    iEps = 1 / (ky**2 * prt.vti ** 2) * Zip(iChi)
    eEps = prt.ope ** 2 / (ky**2 * prt.vte ** 2) * Zep(eChi)
    #eEps = -1 / (ky**2 )
    #iEps =  (1 + 3/2*ky**2/omg**2 * prt.vti**2)/omg**2
    #iEps =  1/omg**2


    return 1 - iEps - eEps

def analytic_IAW(k, ti=0):
    return np.sqrt(1/(1 + 1/(k**2))*(1 + 3*ti*(1+k**2)))


def analytic_IAW_simple(k, ti=0):
    return np.sqrt(1/(1 + 1/(k**2)))
