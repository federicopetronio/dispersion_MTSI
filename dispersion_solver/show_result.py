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
parser.add_argument('--Lt')
parser.add_argument('--Lr')
parser.add_argument('--gam')
parser.add_argument('--ome')
parser.add_argument('--per')



case = parser.parse_args()
Lr = float(case.Lr)/100*u.m
Lt = float(case.Lt)/100*u.m

omega_exp = float(case.ome)
gamma_exp = float(case.gam)
per = int(case.per)
omega_err = omega_exp/2 /per


print(Lr,Lt)
density = 5e16
prt_base=PlasmaParameters(plasmaDensity=density*u.m**(-3),
                    electronTemperature=10*u.eV,
                    magneticField=0.02*u.T,
                    electricField=1e4*u.V/u.m,
                    ionTemperature=0.5*u.eV)

kz = np.pi/Lr*prt_base.Debye_length
kz = 2*np.pi/Lr*prt_base.Debye_length

ky1 = 2*np.pi/Lt*prt_base.Debye_length

kymin = 0.001
kymax = 0.22
pas = 0.0002383025027203481
kappa = np.arange(0.001,0.2200,0.0002383025027203481)
Nkys = len(kappa)

if kz > 0.0517 :
    kymax = 0.44001
    pas = 0.0002381025027203481

kappa = np.arange(kymin,kymax,pas)

gammax = np.ones(Nkys)
omegax = np.ones(Nkys)

current = os.getcwd()
path = current + "/dispersion_data/change_n/{:}/".format(density)

omegax, gammax, kz = precedent_openfile(kz, Nkys=Nkys, path=path)

gamma = gammax[:len(kappa)]
omega = omegax[:len(kappa)]

plt.figure(figsize=(5,4))
plt.plot(kappa,abs(gamma),color="magenta",label="$\\gamma$")
plt.plot(kappa,abs(omega),color='blue',label="$\\omega$")

plt.legend()
plt.title("Radial wave number, $k_r\cdot\\lambda_D = ${:.4f}".format(kz) + ", Lr = {:.2f}".format(Lr*100))
plt.xlabel("Azimuthal wave number, $k_{\\theta} \cdot \\lambda_D$")
plt.ylabel("$(\\gamma, \\omega_r) / \\omega_{pi}$ ")

[plt.axvline(x=xfct, linestyle=(3,(3,6)),color="tab:red") for xfct in [ky1,ky1*2,ky1*3]]
[plt.axvline(x=xfct, linestyle=(3,(3,6)),color="tab:red") for xfct in [ky1,ky1*2,ky1*3,ky1*4]]

# [plt.axvline(x=xfct, linestyle=(3,(3,6)),color="tab:red") for xfct in [ky1,ky1*2]]
# [plt.axvline(x=xfct, linestyle=(3,(3,6)),color="tab:red") for xfct in [ky1]]

# plt.plot(ky1,gamma_exp,"v",color="magenta")
# plt.plot(ky1,omega_exp,"^",color='blue')
plt.errorbar(ky1,omega_exp,yerr=omega_err,ecolor = 'blue', zorder = 1, capsize = 2)



plt.grid(True)
plt.tight_layout()
# plt.savefig('/home/petronio/Nextcloud/theseLPP/reports_vari/MTSI/images_DR/'+ "{:}_{:}.png".format(Lr,Lt))
plt.show()
