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

path = '/home/petronio/Nextcloud/theseLPP/runs/runs_benchmark/MTSI/dispersion_MTSI/dispersion_solver/dispersion_data/general_results/'
# kappa = np.genfromtxt(path + "ky.txt", delimiter="  ")
kz = 0.0405

if kz < 0.0099:
    start = 0.001
    stop = 0.1
    steps = 500
    kappa = np.arange(start,stop,(stop-start)/steps)
else:
    start = 0.001
    stop = 0.20
    steps = 668+84+84+84
    kappa = np.arange(start,stop,(stop-start)/steps)

# print(kappa[:10])
omega = np.genfromtxt(path + "kz={:5.4f}".format(kz) + "_omega_r.txt", delimiter="  ", unpack=False)
# print(omega[:10])
gamma = np.genfromtxt(path + "kz={:5.4f}".format(kz) + "_gamma.txt", delimiter="  ", unpack=False)
# print(gamma[:10])
solution_prec = np.ones(len(kappa))*omega[-1]+1j*gamma[-1]
for k in kappa :
    for index,ka in enumerate(kappa) :
        if ka >= k :
            if index < len(omega):
                solution_prec[index] = complex(omega[index],gamma[index])

plt.figure()
plt.plot(kappa, solution_prec.imag)
plt.plot(kappa,solution_prec.real)
# # plt.show()

for index, k in enumerate(kappa) :
    solution_prec[index] = precedent_openfile(k,kz)

plt.plot(kappa, solution_prec.imag,"--",color="magenta")
plt.plot(kappa,solution_prec.real,"--")
plt.show()
