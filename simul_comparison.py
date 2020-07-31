from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import matplotlib.pyplot as plt
from util.iaw import precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion,precedent_openfile
from util.parameters import PlasmaParameters
from astropy.constants import m_e, m_p
from astropy import units as u
import os


kz = 0.104
kz = 0.0258
density = 5e16
# max_pos = verification_dispersion(kz, density=density,unnorm=True)

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

# densities = [5e16,1e17,2e17,3e17]
densities = [5e16,2e17]

kappa = np.ones((len(densities),Nkys))
gamma = np.ones((len(densities),Nkys))
omega = np.ones((len(densities),Nkys))

current = os.getcwd()
path = current + "/dispersion_data/change_n/{:}/".format(density)

gamma, omega, kz = precedent_openfile(kz_z, Nkys=Nkys, path=path)

plt.figure(figsize=(8,6))
plt.plot(kappa,abs(gamma),"--")
plt.legend()
plt.xlabel("Azimuthal wave number $k_{\\theta}$")
plt.ylabel("Growth rate  $\\gamma$ ")

plt.figure(figsize=(8,6))
plt.plot(kappa,abs(omega),"--")
plt.legend()
plt.xlabel("Azimuthal wave number $k_{\\theta}$")
plt.ylabel("Growth rate  $\\omega$ ")

plt.grid(True)
currentdir = os.getcwd()
# plt.savefig(currentdir + "/images_dispersion/" + "invariance_density.png")
plt.show()
plt.close()
