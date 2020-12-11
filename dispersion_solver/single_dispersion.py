from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import matplotlib.pyplot as plt
import sys
from util.iaw import precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion,precedent_openfile
from util.parameters import PlasmaParameters
from astropy.constants import m_e, m_p
from astropy import units as u
import os
import rcparams

import argparse

""" TRACE THE DR WITH IN IMPUT THE DIMENSIONS IN CM"""

Lr = float(sys.argv[1])*0.01*u.m
L_theta = float(sys.argv[2])*0.01*u.m

density = 5e16
prt_base=PlasmaParameters(plasmaDensity=density*u.m**(-3),
                    electronTemperature=10*u.eV,
                    magneticField=0.02*u.T,
                    electricField=1e4*u.V/u.m,
                    ionTemperature=0.5*u.eV)

ky1 = 2*np.pi/L_theta*prt_base.Debye_length
kz = 2*np.pi/Lr*prt_base.Debye_length
print("ktheta = {:}".format(ky1/prt_base.Debye_length))
print("kr = {:}".format(kz/prt_base.Debye_length))



kymin = 0.001
kymax = 0.22001
pas = 0.0002383025027203481
if kz > 0.0517 :
    kymax = 0.44001
    pas = 0.0002381025027203481

kappa = np.arange(kymin,kymax,pas)
Nkys = len(kappa)
gammax = np.ones(Nkys)
omegax = np.ones(Nkys)

current = os.getcwd()
path = current + "/dispersion_data/change_n/{:}/".format(density)
try:
    os.mkdir(current+"/graphs_report/")
except:
    print("existing folder Tio")
# path = current + "/dispersion_data/change_n_E/20000.0_2e+17/"

omegax, gammax, kz = precedent_openfile(kz, Nkys=Nkys, path=path)
print(kz)

gamma = gammax[:len(kappa)]
omega = omegax[:len(kappa)]

plt.figure(figsize=(6,5))
plt.plot(kappa,abs(gamma),color="magenta",label="$\\gamma$")
# plt.plot(ky1,gamma_exp,"v",color="magenta")
plt.plot(kappa,abs(omega),color='blue',label="$\\omega$")
# plt.plot(ky1,omega_exp,"^",color='blue')

plt.legend()
plt.title("Lt={:.2f} cm, Lr={:.2f} cm,".format(L_theta*100/u.m,Lr*100/u.m)+" $k_r\cdot\\lambda_D = ${:.4f}".format(kz))
plt.xlabel("Azimuthal wave number, $k_{\\theta} \cdot \\lambda_D$")
plt.ylabel("$(\\gamma, \\omega_r) / \\omega_{pi}$ ")

# ky1 = 2*np.pi/L_theta*prt_base.Debye_length
# print(ky1,prt_base.Debye_length)
# [plt.axvline(x=xfct, linestyle=(3,(3,6)),color="tab:red") for xfct in [ky1,ky1*2,ky1*3]]
# [plt.axvline(x=xfct, linestyle=(3,(3,6)),color="tab:red") for xfct in [ky1,ky1*2]]
[plt.axvline(x=xfct, linestyle=(3,(3,6)),color="tab:red") for xfct in [ky1]]



# [plt.axvline(x=xfct, linestyle='dashed') for xfct in kys/u.m]
plt.grid(True)
plt.savefig(current + "/graphs_report/" + "Lt={:.2f}_Lr={:.2f}.png".format(L_theta*100,Lr*100))
plt.show()
plt.close()
