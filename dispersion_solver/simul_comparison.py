from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "STIXGeneral"
from util.iaw import precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion,precedent_openfile
from util.parameters import PlasmaParameters
from astropy.constants import m_e, m_p
from astropy import units as u
import os

caso = 8
ky = 0.104
kz = 0.0258
gamma_exp = 0.091
omega_exp = 0.3813

density = 5e16
kymin = 0.001
kymax = 0.22
pas = 0.0002383025027203481
kappa = np.arange(0.001,0.2200,0.0002383025027203481)
Nkys = len(kappa)
print(Nkys)

prt_base=PlasmaParameters(plasmaDensity=density*u.m**(-3),
                    electronTemperature=10*u.eV,
                    magneticField=0.02*u.T,
                    electricField=1e4*u.V/u.m,
                    ionTemperature=0.5*u.eV)

# Lz = 0.02*u.m
# kz = 2*np.pi/Lz*prt_base.Debye_length
# print(kz)
# ky = np.pi*100/u.m*prt_base.Debye_length*1000/511
# print(ky)
# kys = ky*[1,2,3,4]*u.m
# print(kys)


gamma = np.ones(Nkys)
omega = np.ones(Nkys)

current = os.getcwd()
path = current + "/dispersion_data/change_n/{:}/".format(density)
# path = current + "/dispersion_data/change_n_E/20000.0_2e+17/"

omega, gamma, kz = precedent_openfile(kz, Nkys=Nkys, path=path)

plt.figure(figsize=(8,6))
plt.plot(kappa,abs(gamma),color="magenta",label="$\\gamma$")
plt.plot(ky,gamma_exp,"v",color="magenta")
plt.plot(kappa,abs(omega),color='blue',label="$\\omega$")
plt.plot(ky,omega_exp,"^",color='blue')

plt.legend()
plt.title("Radial wave number: {:.4f}".format(kz))
plt.xlabel("Azimuthal wave number $k_{\\theta} \cdot \\lambda_D$")
plt.ylabel("Growth rate $\\gamma / \\omega_{pi}$ and frequency $ \\omega / \\omega_{pi}$ ")

[plt.axvline(x=xfct*prt_base.Debye_length, linestyle='dashed') for xfct in [492.79884762,  985.59769524, 1478.39654287, 1971.19539049]/u.m]
# [plt.axvline(x=xfct, linestyle='dashed') for xfct in kys/u.m]
plt.grid(True)
plt.savefig(current + "/images_dispersion/" + "{:}_case_ky={:.3f}_kz={:.3f}.png".format(caso,ky,kz))
plt.show()
plt.close()
