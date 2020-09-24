from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import matplotlib.pyplot as plt
from util.iaw import  precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion,precedent_openfile
from astropy import units as u
from util.parameters import PlasmaParameters
import os


Lr=0.0128*u.m
L_theta = 0.01345*u.m
kz = 2*np.pi/Lr

kymin = 0.001
kymax = 0.22
pas = 0.0002383025027203481
kappa = np.arange(0.001,0.2200,0.0002383025027203481)
Nkys = len(kappa)

densities = [5e16,2e17]
densities = [5e16,1e17,2e17,3e17]

kappa = np.ones((len(densities),Nkys))
gamma = np.ones((len(densities),Nkys))
omega = np.ones((len(densities),Nkys))

for index,den in enumerate(densities):
    prtd=PlasmaParameters(plasmaDensity=den*u.m**(-3),
                        electronTemperature=10*u.eV,
                        magneticField=0.02*u.T,
                        electricField=1e4*u.V/u.m,
                        ionTemperature=0.5*u.eV)
    kz_z = kz*prtd.Debye_length
    print(kz_z)
    kappa[index,:], gamma[index,:], omega[index,:] = verification_dispersion(kz_z, density=den,unnorm=False)

plt.figure(figsize=(8,4))
plt.subplot(1,2,1)
for index,dindondensity in enumerate(densities):
    prtd=PlasmaParameters(plasmaDensity=dindondensity*u.m**(-3),
                        electronTemperature=10*u.eV,
                        magneticField=0.02*u.T,
                        electricField=1e4*u.V/u.m,
                        ionTemperature=0.5*u.eV)
    plt.plot(kappa[index,:],abs(gamma[index,:])*prtd.ionPlasmaFrequency,linestyle=(0, (3, 3)),linewidth=1.5,label=str(dindondensity*u.m**(-3)))
# plt.legend()
plt.grid(True)
plt.text(kappa[3,-10],np.amax(gamma)*prtd.ionPlasmaFrequency*u.s/u.rad,"(a)")
plt.text(kappa[3,-10],2e7,"(a)")
plt.xlabel("Azimuthal wave number $k_{\\theta}\lambda_D$")
plt.ylabel("Growth rate  $\\gamma$ 1/m ")

# plt.figure(figsize=(8,6))

plt.subplot(1,2,2)
for index,dindondensity in enumerate(densities):
    prtd=PlasmaParameters(plasmaDensity=dindondensity*u.m**(-3),
                        electronTemperature=10*u.eV,
                        magneticField=0.02*u.T,
                        electricField=1e4*u.V/u.m,
                        ionTemperature=0.5*u.eV)
    plt.plot(kappa[index,:],abs(omega[index,:])*prtd.ionPlasmaFrequency,linestyle=(0, (3, 3)),label=str(dindondensity*u.m**(-3)))
plt.legend()
# plt.text(kappa[3,10],np.amax(omega)*prtd.ionPlasmaFrequency*u.s/u.rad,"(b)")
plt.text(kappa[3,-10],4e7,"(b)")
plt.xlabel("Azimuthal wave number $k_{\\theta}\lambda_D$")
plt.ylabel("Frequency  $\\omega$ 1/m")

plt.grid(True)
currentdir = os.getcwd()
plt.savefig(currentdir + "/images_dispersion/" + "invariance_density_norm.png")
# plt.show()

plt.close()
