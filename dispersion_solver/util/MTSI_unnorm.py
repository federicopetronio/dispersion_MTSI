
from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
from .parameters import PlasmaParameters
import numpy as np
from astropy import units as u

"Normalized permittivity"




def eps_MTSI_unnorm(omg, kx, ky, kz, prt=PlasmaParameters()):
    """Plamsa permitivity for the MTSI, correspondr to the function to solve to find the dispertion relation
    :param omg: (complex) frequency $\omega$, normalized by the ion plasma frequency
    :param kx: wave number in the $x$ direction, normalised by the Debye length. Is supposed to be 0.
    :param ky: wave number in the $y$ direction, normalised by the Debye length. $y$ is the direction of the $ExB$ drift (azymuthal direction in Hall Effect thrusters)
    :param kz: wave number in the $y$ direction, normalised by the Debye length. $z$ is the direction parallel to the magnetif field lines (radial direction in Hall effect thrusters)
    :param prt: Plasma parameters,
    :return:
    """

    # ky = ky*u.m**(-1)
    # kz = kz*u.m**(-1)
    # omg = omg*u.rad*u.s**(-1)

    k2 = kz ** 2 + ky ** 2
    # print(prt.electronCyclotronFrequency) #1/s
    # print(prt.ionPlasmaFrequency) #rad/s


    if k2 == 0:
        raise RuntimeError("The wave vector is Zero !!")

    iEps = prt.ionPlasmaFrequency**2/omg**2
    eEpsz = prt.electronPlasmaFrequency**2 * ( kz**2 ) / ( (omg - ky * prt.driftSpeed*u.rad)**2 * k2 )
    eEpsy = prt.electronPlasmaFrequency**2 * ( ky**2 ) / ( ((omg - ky * prt.driftSpeed*u.rad)**2 - (prt.electronCyclotronFrequency*u.rad)**2)* k2 )
    # print(eEpsz)
    return 1 - iEps - eEpsz - eEpsy
