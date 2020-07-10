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

kx = 0.0
prt=PlasmaParameters(plasmaDensity=5e16/u.m**3,
                     electronTemperature=10*u.eV,
                     magneticField=0.02*u.T,
                     electricField=1e4*u.V/u.m,
                     ionTemperature=0.5*u.eV)

kz = 0.05
ky = 0.05

# gioco_x = np.arange(0.01,0.025,0.0005)
# gioco_y = np.arange(0.255,0.275,0.0005)
gioco_x = np.arange(0.01,3,0.01)
gioco_y = np.arange(0.02,3,0.01)


solution_real = np.zeros((len(gioco_x),len(gioco_y)))
solution_imag = np.zeros((len(gioco_x),len(gioco_y)))
plasmaEps = partial(eps_MTSI, prt=prt) #assign to the function eps_MTSI the value of prt from now on

for i,omega in enumerate(gioco_x) :
    for j,gamma in enumerate(gioco_y) :
        # solution[i,j] = plasmaEps(omg=omega+1j*gamma,kx=0.0,kz=kz,ky=ky)
        zia = 1/plasmaEps(omg=omega+1j*gamma,kx=0.0,kz=kz,ky=ky)
        # print(zia)
        solution_real[i,j] = zia.real
        solution_imag[i,j] = zia.imag

# plt.figure()
# plt.pcolor(gioco_y,gioco_x, solution_real)
# plt.colorbar()
#
# plt.figure()
# plt.pcolor(gioco_y,gioco_x, solution_imag)
# plt.colorbar()

plt.figure()
plt.pcolor(gioco_y,gioco_x, abs(solution_real+1j*solution_imag))
plt.xlabel("k_r")
plt.ylabel("k_z")
plt.colorbar()


plt.show()
