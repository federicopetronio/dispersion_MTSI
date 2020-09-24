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

density = 1e17

# print(len(kys))
# print(Nkys)
prt_base=PlasmaParameters(plasmaDensity=density*u.m**(-3),
electronTemperature=10*u.eV,
magneticField=0.02*u.T,
electricField=1e4*u.V/u.m,
ionTemperature=0.5*u.eV)


Lr = 0.01*u.m
Ltheta = 0.006*u.m

kz0 = 2*np.pi/Lr
ktheta0 = 2*np.pi/Ltheta

print("kz0 = ",kz0)
print("ktheta0 = ",ktheta0)

fig = plt.figure(figsize=(6,5))
plt.title("n = {:}, kz ={:5.0f}, $L_r$ = {:.4f}, $L_t$ = {:.4f}".format(density,kz0,Lr,Ltheta) )
plt.grid()

path1 = sentierino + "/dispersion_data/change_n/{:}/".format(density)
path1 = sentierino + "/dispersion_data/change_n_E/20000.0_2e+17/"

print(path1)

kz = kz0*prt_base.Debye_length
print("kz_tilde  = ", kz0*prt_base.Debye_length)

kymin = 0.001
kymax = 0.22001
pas = 0.0002383025027203481
if kz > 0.0516 :
    kymax = 0.44001
    pas = 0.0002381025027203481

kys = np.arange(kymin,kymax,pas)
Nkys = len(kys)

print(Nkys)


omega1, gamma1, kz = precedent_openfile(kz,Nkys,path1)

gamma1 = gamma1*prt_base.ionPlasmaFrequency
kys_denorm = kys/prt_base.Debye_length

    # plt.plot(kys_denorm[:-5],abs(gamma1[:-5]),label = libello[ind])
plt.plot(kys_denorm[:],abs(gamma1[:]),label = "$\gamma$")

plt.xlabel("Azimuthal wave number, $k_{\\theta}$ 1/m")
plt.ylabel("Growth rate, $\\gamma$ rad/s")
plt.axvline(x = 1*ktheta0*u.m, linestyle='dashed', label="1 period theta")
plt.axvline(x = 2*ktheta0*u.m, linestyle='dashed', label="2 periods theta",color="green")
plt.axvline(x = 3*ktheta0*u.m, linestyle='dashed', label="3 periods theta",color='magenta')
plt.legend()
plt.tight_layout()
plt.savefig(sentierino   + "/images_dispersion/dispersion_kz={:.0f}_Lr={:.0f}_Lt={:.0f}.png".format(kz*1000,Lr*1000/u.m,Ltheta*1000/u.m))
plt.show()
plt.close()
