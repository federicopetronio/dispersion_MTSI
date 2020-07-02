

from numerical_z.plasma import PlasmaParameters
from numerical_z.numerical_z import NumericZ

from numerical_z.evdfs import maxw_f
from plasmapy.mathematics import plasma_dispersion_func

import numpy as np
import scipy as sc

"Normalized permittivity"


def eps_ECDI(omg, kx, ky, kz, prt=PlasmaParameters(), Zlocal=plasma_dispersion_func, Nmax = 5):
    def Zp(x):
        return -2 * (1 + x * Zlocal(x))

    kp = np.sqrt(ky ** 2 + kz ** 2)
    eBeta = kp ** 2 * prt.ope ** 2 / prt.oce ** 2
    k2 = kx ** 2 + kp ** 2

    if k2 == 0:
        raise RuntimeError("The wave vector is Zero !!")

    iChi = (omg - kz * prt.vdi) / (np.sqrt(k2) * prt.vti)
    eChi = (omg - ky * prt.vde) / (kx * prt.vte)
    bessSum = np.array([sc.special.iv(n, eBeta) * Zlocal(eChi - n * prt.oce / (kx * prt.vte)) for n in range(-Nmax, Nmax)])

    iEps = 1 / (k2 * prt.vti ** 2) * Zp(iChi)
    eEps = 1 / k2 * (1. + np.exp(-eBeta) * eChi * sum(bessSum))

    return 1 - iEps + eEps
