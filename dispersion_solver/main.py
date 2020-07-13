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
from util.iaw import eps_IAW, analytic_IAW, analytic_IAW_simple,first_guess,first_guess_1,precedent_guess
from directsolver import solvekys
from scipy import optimize

import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import m_e, m_p
me = m_e.value
mi = 131*m_p.value

#~~~~~~~~~~~~~~~~~~~~~~

from util.parameters import PlasmaParameters
Te = 10*u.eV
plasmaDensity=1e16 *u.m**(-3)
pp = PlasmaParameters(plasmaDensity=plasmaDensity, electronTemperature=Te)

#~~~~~~~~~~~~~~~~~~~~~~
from datetime import date, datetime

today = date.today()
now = datetime.now()
current = today.strftime("%y%m%d_")+now.strftime("%H:%M:%S")
print("Today's date:", current)
#calculate kz

path = "/home/petronio/Nextcloud/theseLPP/runs/runs_benchmark/MTSI/dispersion_MTSI/dispersion_solver/dispersion_data/"

kx = 0.0
prt=PlasmaParameters(plasmaDensity=5e16/u.m**3,
                     electronTemperature=10*u.eV,
                     magneticField=0.02*u.T,
                     electricField=1e4*u.V/u.m,
                     ionTemperature=0.5*u.eV)

Lr = 0.0128*u.m
kz = 2*np.pi*prt.Debye_length/Lr

os.mkdir(path + current)
os.mkdir(path + current + "/images_dispersion")

f = open(path + current + "/parameters.txt","w+")
f.write("plasma density: " + str(prt.plasmaDensity)+"\r\n")
f.write("Te: " + str(prt.electronTemperature)+"\r\n")
f.write("Ti: " + str(prt.ionTemperature)+"\r\n"+ "\r\n")
f.write("E field: " + str(prt.electricField)+"\r\n"+ "\r\n")
f.write("B field: " + str(prt.magneticField)+"\r\n"+ "\r\n")
f.close()

# kz = 0.005550
# kz = 0.001
print("kz * lambda_d = ",kz)
Nkys = 200

plasmaEps = partial(eps_MTSI, prt=prt) #assign to the function eps_MTSI the value of prt from now on
primo = True
dispersion = np.zeros((11,3,Nkys))
kzetas = np.arange(0.001,0.05,0.005)
# we still use the IAW as a first guess
# wrfunct = lambda k: analytic_IAW(k, ti=prt.ionTemperature/ prt.electronTemperature)

for i,kz in enumerate(kzetas):
    print("kz * lambda_d = ",kz)
    # different first guess
    # USE ONLY WITH kymin=0.01, kymax=2, Nkys=200
    if primo :
        wrfunct = lambda k: first_guess_1(k)
        primo = False
    else :
        wrfunct = lambda k: precedent_guess(k=k,ky=kysref1,ome=dispersion[i-1,1,:],gam=dispersion[i-1,2,:])

    kysref1,  xsref1 = solvekys(plasmaEps, kx=kx, kz=kz, kymin=0.01, kymax=0.5,
                             parall=False, wrfunct=wrfunct,Nkys=Nkys,already_solved=True)

    dispersion[i,0,:] = kysref1
    dispersion[i,1,:] = xsref1[:,0]
    dispersion[i,2,:] = xsref1[:,1]


    f = open(path + current + "/kz=" + str(kz) + "_omega_r.txt","w+")
    for i in range(len(kysref1)):
        f.write(str(xsref1[i,0]) + "  ")
    f.close()

    f = open(path + current + "/kz=" + str(kz) + "_gamma.txt","w+")
    for i in range(len(kysref1)):
        f.write(str(xsref1[i,1]) + "  ")
    f.close()


    # plot
    fig = plt.figure(figsize=(6,5))
    # plt.plot(kysref, xsref[:,0], "b", label="$\omega_r$ solver old")
    # plt.plot(kysref, xsref[:,1], "red", label="$\gamma$   solver old")
    # plt.plot(kysref1, xsref1[:,0], "green", label="$\omega_r$ solver")
    # plt.plot(kysref1, xsref1[:,1], "magenta", label="$\gamma$   solver")
    plt.title("kz = " + str(kz))
    plt.grid()
    plt.plot(kysref1, xsref1[:,0], "green", label="$\omega_r$ solver")
    plt.plot(kysref1, xsref1[:,1], "magenta", label="$\gamma$ solver")
    plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
    plt.ylabel("Pulsations  $\\gamma/\\omega_{pi}$ and $ \\omega_r/\\omega_{pi} $")
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path + current + "/images_dispersion/dispersion_kz=" + str(kz) + ".png")

    plt.figure(figsize=(6,5))
    plt.semilogy(kysref1, xsref1[:,2], "magenta")
    # plt.semilogy(kysref1, xsref1[:,2], "+")
    plt.title("Error in the solution")
    plt.grid(True)
    plt.ylabel("Pulsations  $\\gamma/\\omega_{pi}$ and $ \\omega_r/\\omega_{pi} $")
    for i in np.arange(0,200,20):
        plt.plot(kysref1[i],abs(plasmaEps(omg=xsref1[i,0]+1j * xsref1[i,1],kx=0.0,kz=kz,ky=kysref1[i])),'o',color='blue')
    plt.savefig(path + current + "/images_dispersion/error_kz=" + str(kz) + ".png")

f = open(path + current + "/ky.txt","w+")
for i in range(len(kysref1)):
    f.write(str(kysref1[i]) + "  ")
f.close()

# plt.plot(kysref, prt.driftSpeed/prt.BohmSpeed* np.sqrt(me/mi)*np.sqrt(np.pi/8)*kysref/(1 + kysref**2)**(3/2), "r--", label="$\gamma$ IAW  Analytic")
# plt.plot(kysref, analytic_IAW(kysref, ti=prt.ionTemperature/ prt.electronTemperature), "b--", label="$\omega_r$ IAW Analytic")



plt.show()
