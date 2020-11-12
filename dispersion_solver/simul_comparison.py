from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import matplotlib.pyplot as plt
from util.iaw import precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion,precedent_openfile
from util.parameters import PlasmaParameters
from astropy.constants import m_e, m_p
from astropy import units as u
import os
import rcparams

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--caso',   help='1,2,3,4,5,6,7,8')
case = parser.parse_args()

caso = int(case.caso)

print(caso)
density = 5e16
L_theta = 1.28e-2*u.m
prt_base=PlasmaParameters(plasmaDensity=density*u.m**(-3),
                    electronTemperature=10*u.eV,
                    magneticField=0.02*u.T,
                    electricField=1e4*u.V/u.m,
                    ionTemperature=0.5*u.eV)

prt_AT=PlasmaParameters(plasmaDensity=2e17*u.m**(-3),
                    electronTemperature=10*u.eV,
                    magneticField=0.02*u.T,
                    electricField=2e4*u.V/u.m,
                    ionTemperature=0.5*u.eV)

if caso == 1:
    ky1 = 0.1554
    ky2 = 0.1554
    kz = 0.0516*2
    gamma_exp = 0.015
    omega_exp = 0.0586
    Lr = 0.64e-2*u.m
    kz = np.pi/Lr*prt_base.Debye_length/2
    ky1 = 2*np.pi/L_theta*prt_base.Debye_length

if caso == 2:
    ky1 = 0.1554
    ky2 = 0.1554
    kz = 0.0516/2
    gamma_exp = 0.04
    # omega_exp = 0.609
    omega_exp = 0.36 # only during the growth
    Lr = 1.28e-2*u.m
    ky1 = 2*np.pi/L_theta*prt_base.Debye_length * 3
    omega_err = omega_exp/4 /2

if caso == 3:
    ky1 = 0.104
    ky2 = 0.1554
    kz = 0.0344/2
    gamma_exp = 0.064
    # omega_exp = 0.4538
    omega_exp = 0.23 # only during the growth
    ky1 = 2*np.pi/L_theta*prt_base.Debye_length * 2
    omega_err = omega_exp /2

    Lr = 1.92e-2*u.m

if caso == 4:
    ky1 = 0.104
    ky2 = 0.1554
    kz = 0.0258
    # kz = 0.0258*3/2
    gamma_exp = 0.16
    # omega_exp = 0.4709
    Lr = 2.56e-2*u.m
    omega_exp = 0.48 # only during the growth
    ky1 = 2*np.pi/L_theta*prt_base.Debye_length * 2
    kz = 2*np.pi/(2*Lr) * prt_base.Debye_length
    omega_err = omega_exp/2 /2

    Lr = 2.56e-2*u.m

if caso == 5:
    ky1 = 0.1554
    ky2 = 0.1554
    kz = 0.0516*2
    gamma_exp = 0.022
    omega_exp = 0.0902
    Lr = 0.64e-2*u.m
    kz = np.pi/Lr*prt_base.Debye_length/2
    ky1 = 2*np.pi/L_theta*prt_base.Debye_length

if caso == 6:
    ky1 = 0.155
    ky2 = 0.2072
    kz = 0.0516/2
    gamma_exp = 0.036
    # omega_exp = 0.4778
    omega_exp = 0.4061 # only during the growth
    Lr = 1.28e-2*u.m
    ky1 = 2*np.pi/L_theta*prt_base.Debye_length * 3
    omega_err = omega_exp/3 /2

if caso == 7:
    ky1 = 0.104
    ky2 = 0.1554
    kz = 0.0344/2
    gamma_exp = 0.079
    # omega_exp = 0.2958
    omega_exp = 0.27 # only during the growth
    Lr = 1.92e-2*u.m
    ky1 = 2*np.pi/L_theta*prt_base.Debye_length * 2
    omega_err = omega_exp /2

if caso == 8:
    ky1 = 0.104
    ky2 = 0.1554
    kz = 0.0258
    gamma_exp = 0.091
    # omega_exp = 0.34
    omega_exp = 0.40 # only during the growth
    Lr = 2.56e-2*u.m
    ky1 = 2*np.pi/L_theta*prt_base.Debye_length * 2
    kz = 2*np.pi/(2*Lr) * prt_base.Debye_length
    omega_err = omega_exp/2 /2

if caso == 21:
    ky1 = 0.104
    kz = 0.0172
    gamma_exp = 0.068
    omega_exp = 0.2457
    Lr = 3.84e-2*u.m



if caso == 0:
    # ky1 = 2*np.pi/(0.01*u.m)*prt_AT.Debye_length
    ky1 = 2*np.pi/(0.005*u.m)*prt_AT.Debye_length

    # ky1 = 0.226
    kz = 2*np.pi/(0.03*u.m)*prt_AT.Debye_length
    print("ktheta = {:}".format(ky1/prt_AT.Debye_length))
    gamma_exp = 0.06*prt_base.ionPlasmaFrequency/prt_AT.ionPlasmaFrequency
    omega_exp = 0.5538
    density = 2e17


kymin = 0.001
kymax = 0.22
pas = 0.0002383025027203481
kappa = np.arange(0.001,0.2200,0.0002383025027203481)
Nkys = len(kappa)
print(Nkys)

kymin = 0.001
kymax = 0.22001
pas = 0.0002383025027203481
if kz > 0.0517 :
    kymax = 0.44001
    pas = 0.0002381025027203481

kappa = np.arange(kymin,kymax,pas)
# ky1 = 492.79*prt_base.Debye_length
print("kz = ",kz)
# kz = kz/2
# print("kz_after = ",kz)

# Lz = 0.02*u.m
# kz = 2*np.pi/Lz*prt_base.Debye_length
# print(kz)
# ky = np.pi*100/u.m*prt_base.Debye_length*1000/511
# print(ky)
# kys = ky*[1,2,3,4]*u.m
# print(kys)


gammax = np.ones(Nkys)
omegax = np.ones(Nkys)

current = os.getcwd()
path = current + "/dispersion_data/change_n/{:}/".format(density)
# path = current + "/dispersion_data/change_n_E/20000.0_2e+17/"

omegax, gammax, kz = precedent_openfile(kz, Nkys=Nkys, path=path)
print(kz)

gamma = gammax[:len(kappa)]
omega = omegax[:len(kappa)]

plt.figure(figsize=(5,4))
plt.plot(kappa,abs(gamma),color="magenta",label="$\\gamma$")
plt.plot(ky1,gamma_exp,"v",color="magenta")
plt.plot(kappa,abs(omega),color='blue',label="$\\omega$")
plt.plot(ky1,omega_exp,"^",color='blue')
plt.errorbar(ky1,omega_exp,yerr=omega_err,ecolor = 'blue', zorder = 1, capsize = 2)

plt.legend()
plt.title("Radial wave number, $k_r\cdot\\lambda_D = ${:.4f}".format(kz))
plt.xlabel("Azimuthal wave number, $k_{\\theta} \cdot \\lambda_D$")
plt.ylabel("$(\\gamma, \\omega_r) / \\omega_{pi}$ ")

ky1 = 2*np.pi/L_theta*prt_base.Debye_length
print(ky1,prt_base.Debye_length)
[plt.axvline(x=xfct, linestyle=(3,(3,6)),color="tab:red") for xfct in [ky1,ky1*2,ky1*3,ky1*4]]
[plt.axvline(x=xfct, linestyle=(3,(3,6)),color="tab:red") for xfct in [ky1,ky1*2]]


# [plt.axvline(x=xfct, linestyle='dashed') for xfct in kys/u.m]
plt.grid(True)
plt.tight_layout()
# plt.savefig(current + "/images_dispersion/" + "{:}_case_JAN.png".format(caso,ky1,kz))
# plt.savefig('/home/petronio/Nextcloud/theseLPP/reports_vari/MTSI/images/'+ "{:}_case_JAN.png".format(caso,ky1,kz))
# plt.savefig(current + "/images_dispersion/" + "{:}_case_ky={:.3f}_kz={:.3f}.png".format(caso,ky1,kz))
plt.show()
plt.close()
