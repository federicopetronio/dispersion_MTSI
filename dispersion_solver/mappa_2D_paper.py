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

import argparse

density = 5e16
L_theta = 1.28e-2*u.m
prt_base=PlasmaParameters(plasmaDensity=density*u.m**(-3),
                    electronTemperature=10*u.eV,
                    magneticField=0.02*u.T,
                    electricField=1e4*u.V/u.m,
                    ionTemperature=0.5*u.eV)

kz = 0.001

kymin = 0.001
kymax = 0.22
pas = 0.0002383025027203481

if kz > 0.0517 :
    kymax = 0.44001
    pas = 0.0002381025027203481

kappa = np.arange(kymin,kymax,pas)
Nkys = len(kappa)

kzetas = np.arange(0.0010,0.05,0.0001)
# kzetas = np.arange(0.0010,0.005,0.0001)

current = os.getcwd()
path = current + "/dispersion_data/change_n/{:}/".format(density)

omega_plot = np.zeros((len(kzetas),len(kappa)))
gamma_plot = np.zeros((len(kzetas),len(kappa)))
kz_plot = np.zeros(len(kzetas))

for ind, kz in enumerate(kzetas) :
    omega_plot[ind,:], gamma_plot[ind,:], kz_plot[ind] = precedent_openfile(kz, Nkys=Nkys, path=path)
omega_plot = abs(omega_plot)
gamma_plot = abs(gamma_plot)

# plt.figure()
# plt.pcolor(kappa,kz_plot,gamma_plot[:,:])
# plt.colorbar()
# plt.xlabel("Azimuthal wave number, $k_{\\theta} \cdot \\lambda_D$")
# plt.ylabel("Radial wave number, $k_{r} \cdot \\lambda_D$")
#
# plt.figure()
# plt.pcolor(kappa,kz_plot,omega_plot[:,:])
# plt.colorbar()
# plt.xlabel("Azimuthal wave number, $k_{\\theta} \cdot \\lambda_D$")
# plt.ylabel("Radial wave number, $k_{r} \cdot \\lambda_D$")


index = 334
print(kzetas[index-5:index+5])
# index = 30

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.set_xlabel("Azimuthal wave number, $k_{\\theta} \cdot \\lambda_D$")
ax1.set_ylabel("Radial wave number, $k_{r} \cdot \\lambda_D$")
im = ax1.pcolor(kappa,kz_plot,gamma_plot[:,:])

cax = fig.add_axes([0.2, 0.9, 0.3, 0.05])
# cax = fig.add_axes([0.37, 0.9, 0.3, 0.05])

fig.colorbar(im,cax=cax, orientation='horizontal')
ax2.set_ylabel("$\\gamma/\\omega_{pi}$"+", kz = {:.3f}".format(kzetas[index]),color='red')
ax2.plot(kappa,gamma_plot[index,:],color='red')
fig.tight_layout()
plt.savefig("/home/petronio/Downloads/gamma2D.png")


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.set_xlabel("Azimuthal wave number, $k_{\\theta} \cdot \\lambda_D$")
ax1.set_ylabel("Radial wave number, $k_{r} \cdot \\lambda_D$")
im = ax1.pcolor(kappa,kz_plot,omega_plot[:,:])

cax = fig.add_axes([0.2, 0.9, 0.3, 0.05])
# cax = fig.add_axes([0.37, 0.9, 0.3, 0.05])

fig.colorbar(im,cax=cax, orientation='horizontal')
ax2.set_ylabel("$\\omega_r/\\omega_{pi}$"+", kz = {:.3f}".format(kzetas[index]),color='red')
ax2.plot(kappa,omega_plot[index,:],color='red')
fig.tight_layout()
plt.savefig("/home/petronio/Downloads/omega2D.png")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig, axes = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=False)
cmap = "viridis"
fig.set_figheight(4.1)
fig.set_figwidth(10.5)

color1="red"
ax1 = axes[0]
ax2 = ax1.twinx()
ax2.tick_params(colors=color1)


ax1.set_xlabel("Azimuthal wave number, $k_{\\theta} \cdot \\lambda_D$")
ax1.set_ylabel("Radial wave number, $k_{r} \cdot \\lambda_D$")
im = ax1.pcolor(kappa,kz_plot,gamma_plot[:,:],cmap=cmap)

cax = fig.add_axes([0.12, 0.9, 0.13, 0.05])
# cax = fig.add_axes([0.37, 0.9, 0.3, 0.05])

cbar = fig.colorbar(im,cax=cax, orientation='horizontal')
cbar.set_ticks([0,0.5])

ax2.set_ylabel("$\\gamma/\\omega_{pi}$"+", for $k_r \cdot \\lambda_D$ = {:.3f}".format(kzetas[index]),color=color1)
ax2.plot(kappa,gamma_plot[index,:],color=color1)
# fig.tight_layout()


# fig, ax1 = plt.subplots(1,2,2)
color2 = "blue"
ax3 = axes[1]

ax4 = ax3.twinx()
ax4.tick_params(colors=color2)
# ax4.major_ticklabels.set_color(color2)

ax3.set_xlabel("Azimuthal wave number, $k_{\\theta} \cdot \\lambda_D$")
ax3.set_ylabel("Radial wave number, $k_{r} \cdot \\lambda_D$")
im = ax3.pcolor(kappa,kz_plot,omega_plot[:,:],cmap=cmap)

cax = fig.add_axes([0.62, 0.9, 0.13, 0.05])


cbar = fig.colorbar(im,cax=cax, orientation='horizontal')
cbar.set_ticks([0,0.5,1])

ax4.set_ylabel("$\\omega_r/\\omega_{pi}$"+", for $k_r \cdot \\lambda_D$ = {:.3f}".format(kzetas[index]),color=color2)
ax4.plot(kappa,omega_plot[index,:],color=color2)
fig.tight_layout()


plt.tight_layout()
plt.savefig("/home/petronio/Downloads/gamma_omega2D.png")
plt.show()
plt.close()
