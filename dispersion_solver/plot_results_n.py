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
import pylab
from matplotlib import rcParams

import argparse

import rcparams

from astropy.constants import m_e, m_p
me = m_e.value
mi = 131*m_p.value

# "PROBLEM IN OPENING THE FILES WHERE THE KY ARE NOT EXACTLY THE ONES EXPECTED"
#~~~~~~~~~~~~~~~~~~~~~~

from util.parameters import PlasmaParameters
#~~~~~~~~~~~~~~~~~~~~~~
from datetime import date, datetime
sentierino = os.getcwd()

denormalize = True

parser = argparse.ArgumentParser()
parser.add_argument('--L_T', type=float)
parser.add_argument('--L_r', type=float)
# parser.add_argument('--gammap', type=float)
# parser.add_argument('--omegap', type=float)

args = parser.parse_args()

# gamma_pic = 0
L_theta = float(args.L_T)*0.01*u.m
Lr = float(args.L_r)*0.01*u.m
# gamma_pic = float(args.gammap)
# omega_pic = float(args.omegap)


density = 5e16
# density = 2e17


prt_base=PlasmaParameters(plasmaDensity=density*u.m**(-3),
                    electronTemperature=10*u.eV,
                    magneticField=0.02*u.T,
                    electricField=3e4*u.V/u.m,
                    ionTemperature=0.5*u.eV)

ky1 = 2*np.pi/L_theta*prt_base.Debye_length
kz0 = 2*np.pi/Lr/2
ktheta0 = 2*np.pi/L_theta

print("kz0 = ",kz0)
print("ktheta0 = ",ktheta0)

fig = plt.figure(figsize=(6,5))
# plt.title("n = {:}, kz ={:5.0f}, $L_r$ = {:.4f}, $L_t$ = {:.4f}".format(density,kz0,Lr,Ltheta) )
plt.grid()

path1 = sentierino + "/dispersion_data/change_n/{:}/".format(density)
# path1 = sentierino + "/dispersion_data/change_n_E/30000.0_2e+17/"

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

# kys = kys/prt.Debye_length

# plt.plot(kys[:],abs(gamma1[:]),color="magenta",label = "$\gamma$")
# plt.plot(kys[:],abs(omega1[:]),color="blue",label = "$\omega_r$")
# plt.plot(kys[:],kys[:]*50000/prt_base.Debye_length/prt_base.ionPlasmaFrequency,color="blue",label = "ky*50000")


if denormalize :
    gamma1 = gamma1*prt_base.ionPlasmaFrequency
    omega1 = omega1*prt_base.ionPlasmaFrequency
    kys_denorm = kys/prt_base.Debye_length
    kz = kz/prt_base.Debye_length
    ky1=ky1/prt_base.Debye_length*u.m
    plt.plot(kys_denorm[:],abs(gamma1[:]),label = "$\gamma$")
    plt.plot(kys_denorm[:],abs(omega1[:]),label = "$\omega_r$")
    plt.plot(kys_denorm[:],abs(omega1[:]+1j*gamma1),label = "abs_value$")
    plt.plot(kys_denorm[:],kys_denorm[:]*50000,color="blue",label = "ky*ub")
    print("{:e},{:e}".format(np.amax(omega1),gamma1[np.argmax(omega1)]))



plt.title("Radial normalized wave number, $k_r = {:.4f}$".format(kz))
plt.xlabel("Azimuthal wave number, $k_{\\theta}$")
plt.ylabel("($\omega$, $\gamma$) / $\omega_{pi}$")
# plt.axvline(x = 1*ktheta0*u.m, linestyle='dashed', label="1 period theta")
[plt.axvline(x=xfct, linestyle=(3,(3,6)),color="tab:red") for xfct in [ky1,ky1*2,ky1*3]]
# plt.axhline(y=gamma_pic, linestyle=(3,(3,6)),color="tab:blue",label="$\\gamma$ pic")
# plt.axhline(y=omega_pic, linestyle=(3,(3,6)),color="tab:green",label="$\\omega$ pic")


# [plt.axvline(x=xfct, linestyle=(3,(3,6)),color="tab:red") for xfct in [ky1,ky1*2,ky1*3,ky1*4]]

# plt.axvline(x = 3*ktheta0*u.m, linestyle='dashed', label="3 periods theta",color='magenta')
plt.legend()
plt.grid(color='gray', linestyle='dotted')
plt.tight_layout()
# plt.savefig(sentierino   + "/images_keep/dispersion_kz={:.0f}_Lr={:.0f}_Lt={:.0f}.png".format(kz*1000,Lr*10000/u.m,L_theta*10000/u.m))
# plt.savefig('/home/petronio/Nextcloud/theseLPP/reports/MTSI/images/'+ 'strongerE.png')
plt.savefig('/home/petronio/Downloads/comparsion_with_v.png')
plt.show()
plt.close()
