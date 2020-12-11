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

import rcparams

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
pp = PlasmaParameters(plasmaDensity=plasmaDensity, electronTemperature=Te)

#~~~~~~~~~~~~~~~~~~~~~~
from datetime import date, datetime

#calculate kz

sentierino = os.getcwd()

kymin = 0.001
kymax = 0.22

pas = 0.00023803827751196175

Nkys = (kymax-kymin)/pas
Nkys = int(Nkys)

kys = np.arange(kymin,kymax,pas)

kz = 0.0344
fig = plt.figure(figsize=(6,5))
plt.title("$k_r\cdot \\lambda_D$ = {:5.4f}".format(kz))
plt.grid(False)

E_field = [10000.0,15000.0,20000.0,30000.0]
colors = ["black","red","blue","green"]

for ind,den in enumerate(E_field):
    path = sentierino + "/dispersion_data/change_E_Field/{:}/".format(den)

    omega1, gamma1, ksa= precedent_openfile(kz,Nkys+1,path)

    # plt.plot(kysref1, dispersion[i,1,:], "green", label="$\omega_r$ solver")
    # plt.plot(kysref1, dispersion[i,2,:], "magenta", label="$\gamma$ solver")
    plt.plot(kys,abs(gamma1),label = "E = {:.0f}".format(den/1000)+" kV/m",color=colors[ind])
    plt.xlabel("$k_{\\theta}\cdot \\lambda_{De}$")
    plt.ylabel("$\\gamma/\\omega_{pi} $")
    print("okkkk")
    # plt.tight_layout()

    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
plt.legend()
plt.tight_layout()
# plt.savefig(sentierino   + "/images_dispersion/dispersion_kz={:5.4f}_gamma_E.png".format(kz))
# plt.savefig('/home/petronio/Nextcloud/theseLPP/reports/MTSI/images/'+ 'differentE.png')

plt.show()
plt.close()
