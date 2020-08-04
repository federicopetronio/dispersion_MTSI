from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import matplotlib.pyplot as plt
from util.iaw import precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion,precedent_openfile
from util.parameters import PlasmaParameters
from astropy.constants import m_e, m_p
from astropy import units as u
from functools import partial
import os

kz = 0.0370/u.m
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
Lr=0.0128*u.m
kz = 2*np.pi/Lr

densities = [5e16,1e17,2e17,3e17]
densities = [5e16,6e16,7e16,8e16]
# densities = [2e17,5e16]

kappa = np.ones((len(densities),Nkys))
gamma = np.ones((len(densities),Nkys))
omega = np.ones((len(densities),Nkys))

primo =True
kz_z = kz*prt_base.Debye_length

for index,dindondensity in enumerate(densities):
    prt=PlasmaParameters(plasmaDensity=dindondensity*u.m**(-3),
                        electronTemperature=10*u.eV,
                        magneticField=0.02*u.T,
                        electricField=1e4*u.V/u.m,
                        ionTemperature=0.5*u.eV)
    # use this to verify the invariance with respect to the density

    # kz_z = kz_zz[index]
    print("kz_z",kz_z)
    # kz_z = kz/prt_base.Debye_length*prtd.Debye_length
    # print(kz_z)

    current = os.getcwd()
    path = current + "/dispersion_data/change_n/{:}/".format(dindondensity)

    kx = 0.0

    if primo:
        primo = False
        kappa = np.arange(0.001,0.2200,0.0002383025027203481)
        omega,gamma,kz = precedent_openfile(kz_z,path=path)

        kappa = kappa/prt.Debye_length
        omega = omega*prt.ionPlasmaFrequency
        gamma = gamma*prt.ionPlasmaFrequency

        ky,ome,gam = find_max_gamma(kz,path=path)

        print("ky = ",ky,ky/prt.Debye_length)
        gioco_omega = np.arange(0.95*ome,1.05*ome,0.005)
        gioco_gamma = np.arange(0.7*gam,1.05*gam,0.001)

        ky = ky/prt.Debye_length
        kz = kz/prt.Debye_length
        gam = gam*prt.ionPlasmaFrequency
        ome = ome*prt.ionPlasmaFrequency
        gioco_omega = np.arange(0.9*ome/(u.rad/u.s),1.1*ome/(u.rad/u.s),2e4)*u.rad/u.s
        gioco_gamma = np.arange(0.9*gam/(u.rad/u.s),1.1*gam/(u.rad/u.s),2e4)*u.rad/u.s

    print("kz: ",kz)
    print("max_ome: ", ome)
    print("max_gam: ", gam)

    solution_real = np.zeros((len(gioco_omega),len(gioco_gamma)))
    solution_imag = np.zeros((len(gioco_omega),len(gioco_gamma)))

    from util.MTSI_unnorm  import eps_MTSI_unnorm

    plasmaEps = partial(eps_MTSI_unnorm, prt=prt) #assign to the function eps_MTSI the value of prt from now on


    for i,omega_1 in enumerate(gioco_omega) :
        for j,gamma_1 in enumerate(gioco_gamma) :

            zia = 1/plasmaEps(omg=omega_1+1j*gamma_1,kx=0.0,kz=kz,ky=ky)

            solution_real[i,j] = zia.real
            solution_imag[i,j] = zia.imag

    abs_sol=abs(solution_real+1j*solution_imag)

    plt.figure()
    plt.title("invers of susceptibility, "+"density: {:}, n: {:.4f}".format(dindondensity, kz))
    plt.pcolor(gioco_gamma*u.s/u.rad,gioco_omega*u.s/u.rad, abs(solution_real+1j*solution_imag))
    # plt.pcolor(gioco_gamma,gioco_omega, abs(solution_real+1j*solution_imag))
    plt.plot(gam,ome,'*')
    plt.ylabel("$\omega$")
    plt.xlabel("$\gamma$")
    plt.colorbar()

    # plt.figure(figsize=(6,5))
    # plt.title("density {:}".format(plasmaDensity))
    # plt.plot(kappa,abs(gamma), label="solver solution")
    # plt.plot(kyons[0],gioco_gamma[int(max_pos[0,1])],'o',color='blue',label = "computed solution")
    # plt.plot(ky, gam, '*',color='magenta')
    # plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
    # plt.ylabel("Growth rate  $\\gamma/\\omega_{pi}$ ")
    #
    # if unnorm:
    #     plt.xlabel("Azimuthal wave number $k_{\\theta}$")
    #     plt.ylabel("Growth rate  $\\gamma$ ")
    # plt.legend()
    #
    #
    # plt.figure(figsize=(6,5))
    # plt.plot(kappa,omega, label="solver solution")
    # plt.title("density {:}".format(plasmaDensity))
    # plt.plot(kyons[0],gioco_omega[int(max_pos[0,0])],'o',color='blue',label = "computed solution")
    # plt.plot(ky, ome, '*',color='magenta')
    # # plt.plot(kyons[1],gioco_omega[int(max_pos[1,0])],'*',color='blue')
    # # plt.plot(kyons[2],gioco_omega[int(max_pos[2,0])],'*',color='blue')
    # plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
    # plt.ylabel("Pulsations  $\\omega/\\omega_{pi}$ ")
    # if unnorm:
    #     plt.xlabel("Azimuthal wave number $k_{\\theta}$")
    #     plt.ylabel("Pulsations  $\\omega$ ")
    # plt.legend()
plt.show()
