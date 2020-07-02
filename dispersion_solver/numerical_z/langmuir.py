

from numerical_z.plasma import PlasmaParameters
from numerical_z.numerical_z import NumericZ

from numerical_z.evdfs import maxw_f
from plasmapy.mathematics import plasma_dispersion_func

import numpy as np


"Normalized permittivity"

def epsLangmuir_Maxw(omg, kx, ky, kz):
    eChi = omg / (np.sqrt(2) * ky)

    Z_echi = plasma_dispersion_func(eChi)
    Z_echip = -2 * (1 + eChi * Z_echi)

    eEps = 1 / (2 * ky ** 2) * Z_echip

    return (1 - eEps)



def epsLangmuirGene(omg, kx, ky, kz, Z=plasma_dispersion_func):
    eChi = omg / (np.sqrt(2) * ky)

    Z_echi = Z(eChi)
    Z_echip = -2 * (1 + eChi * Z_echi)

    eEps = 1 / (2 * ky ** 2) * Z_echip

    return (1 - eEps)

def langmuir_analytic(k, prt):
    return np.sqrt( prt.ope**2 + 3 * k * prt.vte)
