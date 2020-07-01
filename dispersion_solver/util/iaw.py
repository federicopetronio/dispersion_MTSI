

from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
from .parameters import PlasmaParameters
import numpy as np

"Normalized permittivity"


def eps_IAW(omg, kx, ky, kz, prt=PlasmaParameters(),
            Ze=plasma_dispersion_func,
            Zi=plasma_dispersion_func):
    def Zep(x):
        return -2 * (1 + x * Ze(x))

    def Zip(x):
        return -2 * (1 + x * Zi(x))

    kp = np.sqrt(ky ** 2 + kz ** 2)
    k2 = kx ** 2 + kp ** 2

    if k2 == 0:
        raise RuntimeError("The wave vector is Zero !!")

    iChi = (omg) / (ky * prt.vti/prt.BohmSpeed)
    eChi = (omg - ky * prt.driftSpeed/prt.BohmSpeed) / (ky * prt.vte/prt.BohmSpeed)

    iEps = 1 / (ky**2 * (prt.vti/prt.BohmSpeed) ** 2) * Zip(iChi)
    eEps = (prt.electronPlasmaFrequency/prt.ionPlasmaFrequency) ** 2 / (ky**2 * (prt.vte/prt.BohmSpeed) ** 2) * Zep(eChi)

    return 1 - iEps - eEps

def analytic_IAW(k, ti=0):
    return np.sqrt(1/(1 + 1/(k**2))*(1 + 3*ti*(1+k**2)))


def analytic_IAW_simple(k, ti=0):
    return np.sqrt(1/(1 + 1/(k**2)))
