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

"""COMPARE ANTOINE SOLUTION WITH THE ONE IN THE BENCHMARK"""
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
kymax = 0.22001

pas = 0.0002383025027203481
kys = np.arange(kymin,kymax,pas)
Nkys = len(kys)

# print(len(kys))
# print(Nkys)
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

print("lambda_D_AT: ",prt_AT.Debye_length)

Lr=0.02*u.m
Ltheta=Lr/1000*511

Lr = 0.01*u.m
Ltheta = 0.0026*u.m

Lr = 0.02*u.m
Ltheta = 0.005*u.m

real_Lr = Lr - 0*prt_AT.Debye_length

# Lr=0.0128*u.m
# Ltheta=0.0128*u.m

# 1 period in z direction
kz0 = 2*np.pi/(real_Lr)
# 1/2 period in z direction
# kz0 = np.pi/Lr

ktheta0 = 2*np.pi/Ltheta

print("kz0 = ",kz0)
fig = plt.figure(figsize=(6,5))
plt.title("kz ={:5.0f}, $L_r$ = {:.4f}, $L_t$ = {:.4f}".format(kz0,Lr,Ltheta) )
plt.grid()

path1 = sentierino + "/dispersion_data/change_E_Field/{:}/".format(10000.0)
path2 = sentierino + "/dispersion_data/change_n_E/20000.0_2e+17/"
libello = ["Benchmark", "AT"]


particella = [prt_base,prt_AT]

for ind,path in enumerate([path1,path2]):
    # path = sentierino + "/dispersion_data/change_E_Field/change_n_E/20000.0_2e+17/".format(den)

    prt = particella[ind]
    kz = kz0*prt.Debye_length
    print("kz_tilde  = ", kz0*prt.Debye_length)

    omega1, gamma1, kz = precedent_openfile(kz,Nkys,path)

    # for index in range(len(gamma1)):
    #     if gamma1[index]>100:
    #         gamma1[index] = 1e-12
    gamma1 = gamma1*prt.ionPlasmaFrequency
    kys_denorm = kys/prt.Debye_length

    plt.plot(kys_denorm[:-5],abs(gamma1[:-5]),label = libello[ind])
plt.xlabel("Azimuthal wave number, $k_{\\theta}$ 1/m")
plt.ylabel("Growth rate, $\\gamma$ rad/s")
plt.axvline(x = 1*ktheta0*u.m, linestyle='dashed', label="1 period theta")
plt.axvline(x = 2*ktheta0*u.m, linestyle='dashed', label="2 periods theta",color="green")
plt.axvline(x = 3*ktheta0*u.m, linestyle='dashed', label="3 periods theta",color='magenta')
plt.legend()
plt.tight_layout()
plt.savefig(sentierino   + "/images_dispersion/dispersion_kz={:.0f}_Lr={:.0f}Lt={:.0f}_gamma_AT.png".format(kz*1000,Lr*1000/u.m,Ltheta*1000/u.m))
plt.show()
plt.close()
