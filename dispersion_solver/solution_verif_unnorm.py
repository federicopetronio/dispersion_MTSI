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


kz = 0.0370/u.m
density = 5e16
electricField=1e4
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
                    electricField=electricField*u.V/u.m,
                    ionTemperature=0.5*u.eV)
Lr=0.0128*u.m
L_theta = 0.01345*u.m
kz = 2*np.pi/Lr
kyy = 2*np.pi/L_theta


densities = [5e16,1e17,2e17,3e17]
densities = [5e16]
# densities = [2e17,5e16]
# densities = [1e17,1e17,1e17]

kappa = np.ones((len(densities),Nkys))
gamma = np.ones((len(densities),Nkys))
omega = np.ones((len(densities),Nkys))

for index,dindondensity in enumerate(densities):
    prtd=PlasmaParameters(plasmaDensity=dindondensity*u.m**(-3),
                        electronTemperature=10*u.eV,
                        magneticField=0.02*u.T,
                        electricField=electricField*u.V/u.m,
                        ionTemperature=0.5*u.eV)
    kz_z = kz*prtd.Debye_length
    # use this to verify the invariance with respect to the density

    print("kz_z",kz_z)
    # kz_z = 0.005
    kappa[index,:], gamma[index,:], omega[index,:] = verification_dispersion(kz_z, density=dindondensity,EF=electricField,unnorm=True)

plt.close("all")
import mpltex

plt.figure(figsize=(8,4))

plt.subplot(1,2,1)
linestyles = mpltex.linestyle_generator(lines=['-'], hollow_styles=[])
for index,dindondensity in enumerate(densities):
    # plt.plot(kappa[index,:],abs(gamma[index,:]),linestyle=(0, (3, 3)),linewidth=1.5,label=str(dindondensity*u.m**(-3)),alpha=0.4)
    plt.plot(kappa[index,:],abs(gamma[index,:]),linewidth=1.5,label=str(dindondensity)+" $m^{-3}$",**next(linestyles),alpha=0.6,markevery=50)

# plt.legend()
plt.grid(True)
plt.text(kappa[3,10],np.amax(gamma),"(a)")
plt.xlabel("Azimuthal wave number $k_{\\theta}$  [1/m]")
plt.ylabel("Growth rate  $\\gamma$ [rad/s]")

# plt.figure(figsize=(8,6))

plt.subplot(1,2,2)
linestyles = mpltex.linestyle_generator(lines=['-'], hollow_styles=[])
for index,dindondensity in enumerate(densities):
    plt.plot(kappa[index,:],abs(omega[index,:]),label=str(dindondensity)+" $m^{-3}$",**next(linestyles),alpha=0.6,markevery=50)
plt.legend()
plt.text(kappa[3,10],np.amax(omega),"(b)")
plt.xlabel("Azimuthal wave number $k_{\\theta}$ [1/m]")
plt.ylabel("Frequency  $\\omega$ [rad/s]")
plt.tight_layout()

plt.grid(True)
currentdir = os.getcwd()

plt.show()
plt.close()
