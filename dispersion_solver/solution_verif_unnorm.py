from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import matplotlib.pyplot as plt
from util.iaw import precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion,precedent_openfile
from util.parameters import PlasmaParameters
from astropy.constants import m_e, m_p
from astropy import units as u


kz = 0.0370
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
Lr=0.014*u.m
kz = 2*np.pi/Lr

# densities = [5e16,1e17,2e17,3e17]
densities = [5e16,2e17]
kappa = np.ones((len(densities),Nkys))
gamma = np.ones((len(densities),Nkys))
omega = np.ones((len(densities),Nkys))


for index,dindondensity in enumerate(densities):
    prtd=PlasmaParameters(plasmaDensity=dindondensity*u.m**(-3),
                        electronTemperature=10*u.eV,
                        magneticField=0.02*u.T,
                        electricField=1e4*u.V/u.m,
                        ionTemperature=0.5*u.eV)
    kz_z = kz*prtd.Debye_length
    print("kz_z",kz_z)
    # kz_z = kz/prt_base.Debye_length*prtd.Debye_length
    # print(kz_z)
    kappa[index,:], gamma[index,:], omega[index,:] = verification_dispersion(kz_z, density=dindondensity,unnorm=True)

plt.figure(figsize=(8,6))
for index,dindondensity in enumerate(densities):
    plt.plot(kappa[index,:],abs(gamma[index,:]),"--",label=str(dindondensity))
plt.legend()
plt.xlabel("Azimuthal wave number $k_{\\theta}$")
plt.ylabel("Growth rate  $\\gamma$ ")

plt.figure(figsize=(8,6))
for index,dindondensity in enumerate(densities):
    plt.plot(kappa[index,:],abs(omega[index,:]),"--",label=str(dindondensity))
plt.legend()
plt.xlabel("Azimuthal wave number $k_{\\theta}$")
plt.ylabel("Growth rate  $\\omega$ ")

plt.grid(True)
plt.show()
plt.cl
