
from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
from .parameters import PlasmaParameters
import numpy as np
from astropy import units as u

"Normalized permittivity"


def eps_MTSI(omg, kx, ky, kz, prt=PlasmaParameters()):
    """Plamsa permitivity for the MTSI, correspondr to the function to solve to find the dispertion relation
    :param omg: (complex) frequency $\omega$, normalized by the ion plasma frequency
    :param kx: wave number in the $x$ direction, normalised by the Debye length. Is supposed to be 0.
    :param ky: wave number in the $y$ direction, normalised by the Debye length. $y$ is the direction of the $ExB$ drift (azymuthal direction in Hall Effect thrusters)
    :param kz: wave number in the $y$ direction, normalised by the Debye length. $z$ is the direction parallel to the magnetif field lines (radial direction in Hall effect thrusters)
    :param prt: Plasma parameters,
    :return:
    """

    k2 = kz ** 2 + ky ** 2

    if k2 == 0:
        raise RuntimeError("The wave vector is Zero !!")

    iEps = 1/omg**2
    eEpsz = prt.mi_over_me * ( kz**2 ) / ( (omg - ky * prt.driftSpeed/prt.BohmSpeed)**2 * k2 )
    eEpsy = prt.mi_over_me * ( ky**2 ) / ( ((omg - ky * prt.driftSpeed/prt.BohmSpeed)**2 - prt.electronCyclotronFrequency**2/ (prt.ionPlasmaFrequency/u.rad)**2)* k2 )

    return 1 - iEps - eEpsz - eEpsy
