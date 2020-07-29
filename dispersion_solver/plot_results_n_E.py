from astropy import units as u
from plasmapy.formulary.parameters import plasma_frequency, Debye_length

import directsolver
import os

from functools import partial
from importlib import reload

reload(directsolver)
import util
reload(util)
from util.MTSI  import eps_MTSI
from util.iaw import eps_IAW, analytic_IAW, analytic_IAW_simple,first_guess,first_guess_1,first_guess_mod
from util.iaw import precedent_guess,precedent_guess_mod
from util.tools_dispersion import precedent_openfile
from directsolver import solvekys
from scipy import optimize

import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import m_e, m_p
me = m_e.value
mi = 131*m_p.value

# "PROBLEM IN OPENING THE FILES WHERE THE KY ARE NOT EXACTLY THE ONES EXPECTED"
#~~~~~~~~~~~~~~~~~~~~~~

from util.parameters import PlasmaParameters
Te = 10*u.eV
plasmaDensity=5e16 *u.m**(-3)
pp = PlasmaParameters(plasmaDensity=plasmaDensity, electronTemperature=Te)

#~~~~~~~~~~~~~~~~~~~~~~
from datetime import date, datetime

#calculate kz

sentierino = os.getcwd()

kymin = 0.001
kymax = 0.20

pas = 0.00023803827751196175

Nkys = (kymax-kymin)/pas
Nkys = int(Nkys)

kys = np.arange(kymin,kymax,pas)

kz0 = 0.0258
fig = plt.figure(figsize=(6,5))
plt.title("kz0={:5.4f}".format(kz0))
plt.grid()

path1 = sentierino + "/dispersion_data/change_E_Field/{:}/".format(10000.0)
path2 = sentierino + "/dispersion_data/change_n_E/20000.0_2e+17/"
libello = ["standard", "AT"]

prt_base=PlasmaParameters(plasmaDensity=5e16*u.m**(-3),
                     electronTemperature=10*u.eV,
                     magneticField=0.02*u.T,
                     electricField=1e4*u.V/u.m,
                     ionTemperature=0.5*u.eV)

prt_AT=PlasmaParameters(plasmaDensity=2e17*u.m**(-3),
                     electronTemperature=10*u.eV,
                     magneticField=0.02*u.T,
                     electricField=2e4*u.V/u.m,
                     ionTemperature=0.5*u.eV)

particella = [prt_base,prt_AT]

for ind,path in enumerate([path1,path2]):
    # path = sentierino + "/dispersion_data/change_E_Field/change_n_E/20000.0_2e+17/".format(den)

    prt = particella[ind]
    kz = (kz0/prt_base.Debye_length)*prt.Debye_length

    omega1, gamma1 = precedent_openfile(kz,Nkys,path)


    for index in range(len(gamma1)):
        if gamma1[index]>100:
            gamma1[index] = 1e-12
    gamma1 = gamma1*prt.ionPlasmaFrequency
    kys_denorm = kys/prt.Debye_length
    # plt.plot(kysref1, dispersion[i,1,:], "green", label="$\omega_r$ solver")
    # plt.plot(kysref1, dispersion[i,2,:], "magenta", label="$\gamma$ solver")
    plt.plot(kys_denorm,abs(gamma1),label = libello[ind])
    plt.xlabel("Azimuthal wave number $k_{\\theta} [1/m]$")
    plt.ylabel("Pulsations  $\\gamma [rad/s] $")
    plt.tight_layout()

    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
plt.legend()
plt.tight_layout()
plt.savefig(sentierino   + "/images_dispersion/dispersion_kz={:5.4f}_gamma_AT.png".format(kz))
# plt.show()
plt.close()
