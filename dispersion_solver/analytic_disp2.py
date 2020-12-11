from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import scipy as spy
import matplotlib.pyplot as plt
# from util.iaw import precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion,precedent_openfile
from util.parameters import PlasmaParameters
from util.MTSI_unnorm  import eps_MTSI_unnorm
from functools import partial

from astropy.constants import m_e, m_p
from astropy import units as u
import os
import rcparams

verobse =  False
unnorm = True

"""SOLVE THE APPROX SOLUTION"""


import cubic_solver
plasmaDensity=5e16*u.m**-3
prt=PlasmaParameters(plasmaDensity=plasmaDensity,
                     electronTemperature=10*u.eV,
                     magneticField=0.02*u.T,
                     electricField=1e4*u.V/u.m,
                     ionTemperature=0.5*u.eV)

# ky = np.arange(10,8000,10)
ky = np.arange(0.001,0.2200,0.0002383025027203481)/prt.Debye_length*u.m

kz = 0.011*np.ones(len(ky))/prt.Debye_length*u.m
# kz = 6.28/0.0456*np.ones(len(ky))
v0 = 5e4
Nkys = 920


print("kz = ",kz[0],"ky = ",ky[0])
if verobse: print("kyv0 = {:e} ".format(ky[0]*v0))

k = (kz**2 + ky**2)**0.5

omega_r = np.zeros((len(kz),3), dtype=complex)
omega_tot = np.zeros((len(kz),3), dtype=complex)
gamma = np.zeros((len(kz),3), dtype=complex)

omega_ce = prt.electronCyclotronFrequency*u.s
omega_pe = prt.electronPlasmaFrequency*u.s/u.rad
omega_pi = prt.ionPlasmaFrequency*u.s/u.rad

kappa_solver = np.ones(Nkys)
gamma_solver = np.ones(Nkys)
omega_solver = np.ones(Nkys)

kz_norm = kz[0]*prt.Debye_length/u.m
kappa_solver, gamma_solver, omega_solver = verification_dispersion(kz_norm, density=plasmaDensity*u.m**3,unnorm=True)
omega_cpx = omega_solver + 1j*gamma_solver
omega_cpx = omega_cpx*u.s/u.rad
true_solution = 1 - omega_pi**2/omega_cpx**2 - (kz**2 * omega_pe**2)/(omega_cpx-ky * v0)**2/k**2 - (ky**2 * omega_pe**2)/((omega_cpx-ky * v0)**2 -omega_ce**2)/k**2

plt.figure()
plt.semilogy(ky,abs(true_solution.real))
plt.semilogy(ky,abs(true_solution.imag))

alpha = (1 - (kz**2 * omega_pe**2)/(k**2 * ky**2 * v0**2) + (ky**2 * omega_pe**2)/(k**2 * omega_ce**2) * (1+ky**2 * v0**2 / omega_ce**2))
beta = -(2 * (kz**2 * omega_pe**2)/(k**2 * ky**3 * v0**3) + 2 * (ky**3 * omega_pe**2 * v0)/(k**2 * omega_ce**4))
appr_solution = (alpha*omega_cpx**2 + beta * omega_cpx**3 - omega_pi**2)/omega_cpx**2

appr_derivative_omega = -alpha/(6*beta)
solver_derivative_omega = (omega_solver[1:]-omega_solver[:-1])/(kappa_solver[1]-kappa_solver[0])

plt.figure()
plt.plot(appr_derivative_omega)
plt.plot(solver_derivative_omega)

plt.figure()
plt.semilogy(ky,abs(appr_solution.real))
plt.semilogy(ky,abs(appr_solution.imag))
plt.show()
