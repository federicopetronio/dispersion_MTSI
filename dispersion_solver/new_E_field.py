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
from astropy.constants import m_e, m_p
me = m_e.value
mi = 131*m_p.value

# "PROBLEM IN OPENING THE FILES WHERE THE KY ARE NOT EXACTLY THE ONES EXPECTED"
#~~~~~~~~~~~~~~~~~~~~~~

from util.parameters import PlasmaParameters
Te = 10*u.eV
plasmaDensity=5e16 *u.m**(-3)
electricField = 3e4*u.V/u.m
pp = PlasmaParameters(plasmaDensity=plasmaDensity, electronTemperature=Te)

#~~~~~~~~~~~~~~~~~~~~~~
from datetime import date, datetime

#calculate kz

sentierino = os.getcwd()
try:
    os.mkdir(sentierino + "/dispersion_data")
except:
    print("already existing folder dispersion_data")
try:
    os.mkdir(sentierino + "/dispersion_data/change_E_Field")
except:
    print("already existing folder dispersion_data/change_E_Field")
try:
    os.mkdir(sentierino + "/dispersion_data/change_E_Field/{:}/".format(electricField*u.m**(1)/u.V))
except:
    print("already existing folder dispersion_data/change_E_Field/{:}/".format(electricField*u.m**(1)/u.V))


# path = sentierino + "/dispersion_data/all_data/"
path2 = sentierino + "/dispersion_data/change_E_Field/{:}/".format(electricField*u.m**(1)/u.V)

# path = "/home/petronio/Nextcloud/theseLPP/runs/runs_benchmark/MTSI/dispersion_MTSI/dispersion_solver/dispersion_data/"
print(path2)

kx = 0.0
prt=PlasmaParameters(plasmaDensity=plasmaDensity,
                     electronTemperature=10*u.eV,
                     magneticField=0.02*u.T,
                     electricField=electricField,
                     ionTemperature=0.5*u.eV)

Lr = 0.0128*u.m
kz = 2*np.pi*prt.Debye_length/Lr

kzetas = np.arange(0.490,0.0530,0.002)


try:
    os.mkdir(path2 + "/images_dispersion")
except:
    print("t'appost")


f = open(path2   + "/parameters.txt","w+")
f.write("plasma density: " + str(prt.plasmaDensity)+"\r\n"+ "\r\n")
f.write("Te: " + str(prt.electronTemperature)+"\r\n"+ "\r\n")
f.write("Ti: " + str(prt.ionTemperature)+"\r\n"+ "\r\n")
f.write("E field: " + str(prt.electricField)+"\r\n"+ "\r\n")
f.write("B field: " + str(prt.magneticField)+"\r\n"+ "\r\n")
f.close()

# kz = 0.005550
# kz = 0.001
print("kz * lambda_d = ",kz)
kymin = 0.001
kymax = 0.20
pas = 0.00023803827751196175


Nkys = (kymax-kymin)/pas
Nkys = int(Nkys)

plasmaEps = partial(eps_MTSI, prt=prt) #assign to the function eps_MTSI the value of prt from now on
primo = True



dispersion = np.zeros((len(kzetas),4,Nkys))
dispersion_clean = np.zeros((len(kzetas),4,Nkys))



for i,kz in enumerate(kzetas):
    print("kz * lambda_d = ",kz)

    omega_1, gamma_1 = precedent_openfile(kz=kz,Nkys=Nkys)
    ky_1 = np.arange(kymin,kymax,pas)
    # primo = False

    wrfunct = lambda k: precedent_guess_mod(k=k,ky=ky_1,ome=omega_1,gam=gamma_1)

    kysref1,  xsref1 = solvekys(plasmaEps, kx=kx, kz=kz, kymin=kymin, kymax=kymax,
                             parall=False, wrfunct=wrfunct,Nkys=Nkys,already_solved=True,root=False)

    dispersion[i,0,:] = kysref1
    dispersion[i,1,:] = xsref1[:,0]
    dispersion[i,2,:] = xsref1[:,1]
    dispersion[i,3,:] = xsref1[:,2]

    f = open(path2   + "/kz={:5.4f}".format(kz) + "_omega_r.txt","w+")
    for ind in range(len(kysref1)):
        f.write(str(xsref1[ind,0]) + "  ")
    f.close()

    f = open(path2   + "/kz={:5.4f}".format(kz) + "_gamma.txt","w+")
    for ind in range(len(kysref1)):
        f.write(str(xsref1[ind,1]) + "  ")
    f.close()
    # eliminates the points with a too large error
    for j in np.arange(0,Nkys):
        if dispersion[i,3,j]>0.5 :
            dispersion[i,1,j] = 1e-12
            dispersion[i,2,j] = 1e-12




    ky_1 = dispersion[i,0,:]
    omega_1 = dispersion[i,1,:]
    gamma_1 = dispersion[i,2,:]
    argmax_1 = np.argmax(gamma_1)
    # print("max : ", argmax_1, gamma_1[argmax_1])
    cont = 0
    for j in np.arange(0,Nkys):
        if j<argmax_1:
            if gamma_1[j-cont]<1e-8 :
                ky_1 = np.delete(ky_1,j-cont)
                omega_1 = np.delete(omega_1,j-cont)
                gamma_1 = np.delete(gamma_1,j-cont)
                cont = cont+1
        else:
            if gamma_1[j-cont]<0 :
                gamma_1[j-cont] = -gamma_1[j-cont]

    # plot
    fig = plt.figure(figsize=(6,5))
    plt.title("kz={:5.4f}".format(kz))
    plt.grid()
    # plt.plot(kysref1, dispersion[i,1,:], "green", label="$\omega_r$ solver")
    plt.plot(kysref1, abs(dispersion[i,2,:]), "magenta", label="$\gamma$ solver")
    plt.plot(ky_1,gamma_1,"--")
    plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
    plt.ylabel("Pulsations  $\\gamma/\\omega_{pi} $")
    plt.tight_layout()
    for j in np.arange(0,Nkys):
        if dispersion[i,3,j]>0.5 :
            plt.plot(kysref1[j], dispersion[i,2,j], "*",color="blue")

            #plt.xlim(left=0)
            #plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path2   + "/images_dispersion/dispersion_kz={:5.4f}_gamma.png".format(kz))
    plt.close()



    fig = plt.figure(figsize=(6,5))
    plt.title("kz={:5.4f}".format(kz))
    plt.grid()
    plt.plot(kysref1, dispersion[i,1,:], "green", label="$\omega$ solver")
    plt.plot(ky_1,omega_1,"--")
    plt.plot()
    plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
    plt.ylabel("Pulsations $ \\omega_r/\\omega_{pi} $")
    plt.tight_layout()
    for j in np.arange(0,Nkys):
        if dispersion[i,3,j]>0.5 :
            plt.plot(kysref1[j], dispersion[i,2,j], "*",color="blue")

    plt.tight_layout()
    plt.savefig(path2   + "/images_dispersion/dispersion_kz={:5.4f}_omega.png".format(kz))
    plt.close()

    plt.figure(figsize=(6,5))
    plt.semilogy(kysref1, xsref1[:,2], "magenta")
    # plt.semilogy(kysref1, xsref1[:,2], "+")
    plt.title("Error in the solution")
    plt.grid(True)
    plt.ylabel("Pulsations  $\\gamma/\\omega_{pi}$ and $ \\omega_r/\\omega_{pi} $")
    for ii in np.arange(0,Nkys,20):
        plt.plot(kysref1[ii],abs(plasmaEps(omg=xsref1[ii,0]+1j * xsref1[ii,1],kx=0.0,kz=kz,ky=kysref1[ii])),'o',color='blue')
    plt.savefig(path2   + "/images_dispersion/error_kz={:5.4f}".format(kz) + ".png")
    plt.close()
f = open(path2   + "/ky.txt","w+")
for i in range(len(kysref1)):
    f.write(str(kysref1[i]) + "  ")
f.close()

# plt.plot(kysref, prt.driftSpeed/prt.BohmSpeed* np.sqrt(me/mi)*np.sqrt(np.pi/8)*kysref/(1 + kysref**2)**(3/2), "r--", label="$\gamma$ IAW  Analytic")
# plt.plot(kysref, analytic_IAW(kysref, ti=prt.ionTemperature/ prt.electronTemperature), "b--", label="$\omega_r$ IAW Analytic")

fig = plt.figure(figsize=(12,9))
plt.grid()

for i,kz in enumerate(kzetas):
    gamma_1 = dispersion[i,2,:]
    ky_1 = dispersion[i,0,:]

    argmax_1 = np.argmax(gamma_1)
    print("max : ", argmax_1, gamma_1[argmax_1])
    cont = 0
    for j in np.arange(0,argmax_1):
        if gamma_1[j-cont]<1e-8 :
            gamma_1 = np.delete(gamma_1,j-cont)
            ky_1 = np.delete(ky_1,j-cont)
            cont = cont+1

    # plt.plot(kysref1, dispersion[i,1,:], "green", label="$\omega_r$ solver")
    plt.plot(ky_1, gamma_1, label="$\gamma$ solver, kz = %5.4f"%kz)
    plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
    plt.ylabel("Pulsations  $\\gamma/\\omega_{pi}$ ")
    plt.title("Solution of the dispersion relation")
    plt.legend()
    plt.tight_layout()
plt.savefig(path2   + "/images_dispersion/dispersion_kz_multiple.png")


# plt.show()
