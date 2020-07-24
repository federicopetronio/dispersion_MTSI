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
from util.iaw import precedent_openfile, precedent_guess,precedent_guess_mod
from directsolver import solvekys
from scipy import optimize

import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import m_e, m_p,eps0,e
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

kz0 = 0.0190
fig = plt.figure(figsize=(6,5))
plt.title("kz={:5.4f}".format(kz0))
plt.grid()

density = [5e16,1e17,2e17,3e17]

prt_base=PlasmaParameters(plasmaDensity=5e16*u.m**(-3),
                     electronTemperature=10*u.eV,
                     magneticField=0.02*u.T,
                     electricField=1e4*u.V/u.m,
                     ionTemperature=0.5*u.eV)

for ind,den in enumerate(density):
    path = sentierino + "/dispersion_data/change_n/{:}/".format(den)


    prt=PlasmaParameters(plasmaDensity=den*u.m**(-3),
                         electronTemperature=10*u.eV,
                         magneticField=0.02*u.T,
                         electricField=1e4*u.V/u.m,
                         ionTemperature=0.5*u.eV)


    # kz = kz0*(5e16/den)**0.5
    kz = (kz0/prt_base.Debye_length)*prt.Debye_length
    print("kz : {:5.4f}".format(kz))
    omega1, gamma1 = precedent_openfile(kz,Nkys,path)
    print("omega_tilda ",omega1[np.argmax(gamma1)])
    print("gamma_tilda ",gamma1[np.argmax(gamma1)])
    print("ky_tilda",kys[np.argmax(gamma1)])

    # omega_pi = (den*(1.6e-19)**2/(eps0*mi))**0.5
    omega_pi = prt.ionPlasmaFrequency
    gamma_1 = gamma1 * omega_pi
    print(omega1[np.argmax(gamma1)]* omega_pi)
    print(gamma1[np.argmax(gamma1)]* omega_pi)
    print(kys[np.argmax(gamma1)]/prt.Debye_length,kz/prt.Debye_length)


    kys_denorm = kys/prt.Debye_length
    plt.plot(kys_denorm, abs(gamma_1),"--",label = "n = {:}".format(den)+" m-3")
    plt.xlabel("Azimuthal wave number $k_{\\theta} [1/m]$")
    plt.ylabel("Pulsations  $\\gamma [rad/s]$")
    plt.tight_layout()
    break

plt.legend()
plt.tight_layout()
# plt.savefig(sentierino   + "/images_dispersion/dispersion_kz={:5.4f}_gamma.png".format(kz))
plt.show()
plt.close()
