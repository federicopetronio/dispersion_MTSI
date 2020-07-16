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
plasmaDensity=5e16 *u.m**(-3)
pp = PlasmaParameters(plasmaDensity=plasmaDensity, electronTemperature=Te)

#~~~~~~~~~~~~~~~~~~~~~~
from datetime import date, datetime

today = date.today()
now = datetime.now()
current = today.strftime("%y%m%d_")+now.strftime("%H:%M:%S")+"map2D"
sentierino = os.getcwd()
while True:
    try:
        os.mkdir(sentierino + "/dispersion_data")
    except:
        print("already existing folder dispersion_data")
    try:
        os.mkdir(sentierino + "/dispersion_data/2Dmap")
    except:
        print("existing folder 2Dmap")
        break
path = sentierino + "/dispersion_data/2Dmap/"
# path = "/home/petronio/Nextcloud/theseLPP/runs/runs_benchmark/MTSI/dispersion_MTSI/dispersion_solver/dispersion_data/"
print(path)

kx = 0.0
prt=PlasmaParameters(plasmaDensity=5e16/u.m**3,
                     electronTemperature=10*u.eV,
                     magneticField=0.02*u.T,
                     electricField=1e4*u.V/u.m,
                     ionTemperature=0.5*u.eV)

kz = 0.001
ky = 0.0268

gioco_omega = np.arange(0.02,0.5,0.005)
gioco_gamma = np.arange(0.01,1,0.005)
# gioco_gamma = np.arange(0.02,3,0.01)
# gioco_omega = np.arange(0.01,3,0.01)


solution_real = np.zeros((len(gioco_omega),len(gioco_gamma)))
solution_imag = np.zeros((len(gioco_omega),len(gioco_gamma)))
plasmaEps = partial(eps_MTSI, prt=prt) #assign to the function eps_MTSI the value of prt from now on

for i,omega in enumerate(gioco_omega) :
    for j,gamma in enumerate(gioco_gamma) :
        # solution[i,j] = plasmaEps(omg=omega+1j*gamma,kx=0.0,kz=kz,ky=ky)
        zia = 1/plasmaEps(omg=omega+1j*gamma,kx=0.0,kz=kz,ky=ky)
        # print(zia)
        solution_real[i,j] = zia.real
        solution_imag[i,j] = zia.imag

a=abs(solution_real+1j*solution_imag)
max_pos = np.unravel_index(abs(a).argmax(), a.shape)
print(max_pos)
# plt.figure()
# plt.pcolor(gioco_gamma,gioco_omega, solution_real)
# plt.colorbar()
#
# plt.figure()
# plt.pcolor(gioco_gamma,gioco_omega, solution_imag)
# plt.colorbar()

plt.figure()
plt.title("invers of susceptibility ")
plt.pcolor(gioco_gamma,gioco_omega, abs(solution_real+1j*solution_imag))
plt.xlabel("$\gamma/\omega_{pi}$")
plt.ylabel("$\omega/\omega_{pi}$")
plt.text(x=gioco_gamma[-70],y=gioco_omega[-25],s="kz = %5.4f \n"%kz + "ky = %5.4f \n"%ky+
        "$\omega_{max}$ = %6.4f \n"%gioco_omega[max_pos[0]] + "$\gamma_{max}$ = %6.4f"%gioco_gamma[max_pos[1]],color='red')
# plt.text(x=gioco_gamma[-70],y=gioco_omega[-15],s="$\omega_{max}$ = %6.4f \n"%gioco_omega[max_pos[0]] + "$\gamma_{max}$ = %6.4f"%gioco_gamma[max_pos[1]],color='red')
plt.savefig(path + "kz=%5.4f"%kz + "_ky+%5.4f"%ky+".png")
plt.colorbar()


plt.show()
