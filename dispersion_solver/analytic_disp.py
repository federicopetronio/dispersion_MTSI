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

import cubic_solver

"""SOLVE THE APPROX SOLUTION"""

verobse =  False
unnorm = True


plasmaDensity=5e16*u.m**-3
prt=PlasmaParameters(plasmaDensity=plasmaDensity,
                     electronTemperature=10*u.eV,
                     magneticField=0.02*u.T,
                     electricField=1e4*u.V/u.m,
                     ionTemperature=0.5*u.eV)




# kz = np.arange(0.00001,0.0001,0.0000001)/prt.Debye_length*u.m
# ky = np.arange(0.001,0.8,0.001)/prt.Debye_length*u.m
ky = np.arange(10,8000,10)
kz = 0.011*np.ones(len(ky))/prt.Debye_length*u.m
# kz = 6.28/0.0456*np.ones(len(ky))
v0 = 5e4


print("kz = ",kz[0],"ky = ",ky[0])
if verobse: print("kyv0 = {:e} ".format(ky[0]*v0))
# kz = 0.1/prt.Debye_length *u.m
# ky = 0.025/prt.Debye_length *u.m

k = (kz**2 + ky**2)**0.5

omega_r = np.zeros((len(kz),3), dtype=complex)
omega_tot = np.zeros((len(kz),3), dtype=complex)
gamma = np.zeros((len(kz),3), dtype=complex)

# print(ky, 2*np.pi/(2*0.0128))

omega_ce = prt.electronCyclotronFrequency*u.s
omega_pe = prt.electronPlasmaFrequency*u.s/u.rad
omega_pi = prt.ionPlasmaFrequency*u.s/u.rad

# print(kz,ky,k,v0,omega_ce,omega_pi,omega_pe)
# alpha = (1*u.rad**2 - (kz**2 * omega_pe**2)/(k**2 * ky**2 * v0**2) + (ky**2 * omega_pe**2)/(k**2 * omega_ce**2) * (1+ky**2 * v0**2 / omega_ce**2))

alpha = (1 - (kz**2 * omega_pe**2)/(k**2 * ky**2 * v0**2) + (ky**2 * omega_pe**2)/(k**2 * omega_ce**2) * (1+ky**2 * v0**2 / omega_ce**2))
beta = -(2 * (kz**2 * omega_pe**2)/(k**2 * ky**3 * v0**3) + 2 * (ky**3 * omega_pe**2 * v0)/(k**2 * omega_ce**4))
if verobse: print("From equation (A): alpha = {:e}, beta = {:e}".format(alpha[0],beta[0]))

# a = 8 * beta**2 *u.rad**-4 *u.s**-2
# b = 8 * alpha * beta *u.rad**-4 *u.s**-1
# c = 2 * alpha**2 *u.rad**-4
# d = beta * omega_pi**2 *u.rad**-4 *u.s**1

# print(np.roots([a.value,b.value,c.value,d.value]))
for ind in range(len(kz)):
    a = 8 * beta[ind]**2
    b = 8 * alpha[ind] * beta[ind]
    c = 2 * alpha[ind]**2
    d = beta[ind] * omega_pi**2
    if verobse: print("coefficients of equation (C): {:e},{:e},{:e},{:e}".format(a.value,b.value,c.value,d.value))
    omega_r[ind,:] = np.roots([a.value,b.value,c.value,d.value])
    gamma[ind,:] = ((-omega_pi**2 + omega_r[ind]**2 * alpha[ind] + omega_r[ind]**3 * beta[ind])/(alpha[ind] + 3 * omega_r[ind] * beta[ind]))**0.5
# print("unit",a.unit,b.unit,c.unit,d.unit)
# print("value",a.value,b.value,c.value,d.value)

omega_cpx = omega_r[:,2] + 1j*gamma[:,2]
if verobse: print("omega = {:e} + j{:e}".format(omega_cpx[0].real,omega_cpx[0].imag))
if verobse: print("|omega| = {:e}".format((omega_cpx[0].real**2+omega_cpx[0].imag**2)**0.5))

# print(omega_r)
# print(gamma)

# plt.plot(kz,abs(gamma[:,0]),label='0')
if unnorm == False:
    plt.title("kz = {:.4f}".format(kz[0]*prt.Debye_length/u.m))
else:
    plt.title("kz = {:.0f}".format(kz[0]))

if unnorm == False:
    plt.plot(ky*prt.Debye_length/u.m,abs(gamma[:,2])/omega_pi,label='gamma')
    plt.plot(ky*prt.Debye_length/u.m,abs(omega_r[:,2]/omega_pi),label='omega')
else:
    plt.plot(ky,abs(gamma[:,2]),label='gamma')
    plt.plot(ky,abs(omega_r[:,2]),label='omega')

# plt.plot(kz,abs(gamma[:,2]),label='2')
plt.legend()
plt.show()
# print("DEBYE", prt.Debye_length)

# print(a.value*omega_r[:,0]**3 + b.value*omega_r[:,0]**2 + c.value*omega_r[:,0] + d.value)
# print(a.value*omega_r[:,1]**3 + b.value*omega_r[:,1]**2 + c.value*omega_r[:,1] + d.value)
if verobse: print("(C) solved function for omega: ",a.value*omega_r[:,2]**3 + b.value*omega_r[:,2]**2 + c.value*omega_r[:,2] + d.value)
if verobse: print("(A) solved function for omega/gamma: ",alpha*omega_cpx**2 + beta * omega_cpx**3 - omega_pi**2)
#
# print(alpha*(omega_r[:,0] + 1j*gamma[:,0])**2 + beta * (omega_r[:,0] + 1j*gamma[:,0])**3 - omega_pi**2)
# print(alpha*(omega_r[:,1] + 1j*gamma[:,1])**2 + beta * (omega_r[:,1] + 1j*gamma[:,1])**3 - omega_pi**2)


# omega_cpx[0] = 1e6 + 1j*1e6
true_solution = 1 - omega_pi**2/omega_cpx**2 - (kz**2 * omega_pe**2)/(omega_cpx-ky * v0)**2/k**2 - (ky**2 * omega_pe**2)/((omega_cpx-ky * v0)**2 -omega_ce**2)/k**2
appr_solution = (alpha*omega_cpx**2 + beta * omega_cpx**3 - omega_pi**2)/omega_cpx**2

if verobse: print("true_solution with the calculated omega and gamma: ", true_solution)
if verobse: print("appr_solution with the calculated omega and gamma: ", appr_solution)

# true_solution_short_1 = -(kz**2 * omega_pe**2)/(omega_cpx-ky * v0)**2/k**2 - (ky**2 * omega_pe**2)/((omega_cpx-ky * v0)**2 -omega_ce**2)/k**2
# appr_solution_short_1 = -(kz**2 * omega_pe**2)/(ky**2 * v0**2)*(1+2*omega_cpx/(ky*v0))/k**2 - (ky**2 * omega_pe**2)/(omega_ce**2 * k**2) * (1+ky**2*v0**2/omega_ce**2 - 2*ky*v0*omega_cpx/omega_ce**2)

true_solution_short_1 = -(kz**2 * omega_pe**2)/(omega_cpx-ky * v0)**2/k**2
appr_solution_short_1 = -(kz**2 * omega_pe**2)/(ky**2 * v0**2)*(1+2*omega_cpx/(ky*v0))/k**2

true_solution_short_2 = - (ky**2 * omega_pe**2)/((omega_cpx-ky * v0)**2 -omega_ce**2)/k**2
appr_solution_short_2 =   (ky**2 * omega_pe**2)/(omega_ce**2 * k**2) * (1+ky**2*v0**2/omega_ce**2 - 2*ky*v0*omega_cpx/omega_ce**2)

# appr_solution_short = ((alpha-1)*omega_cpx**2 + beta * omega_cpx**3)/omega_cpx**2
# print("omega {:e} {:e}".format(omega_cpx[0].real,omega_cpx[0].imag))
if verobse: print("true_solution short (a) {:e} {:e}j, (b) {:e} {:e}j".format(true_solution_short_1[0].real,true_solution_short_1[0].imag,true_solution_short_2[0].real,true_solution_short_2[0].imag))
if verobse: print("appr_solution short (a) {:e} {:e}j, (b) {:e} {:e}j".format(appr_solution_short_1[0].real,appr_solution_short_1[0].imag,appr_solution_short_2[0].real,appr_solution_short_2[0].imag))
if verobse: print(true_solution_short_1+true_solution_short_2)
if verobse: print(appr_solution_short_1+appr_solution_short_2)
if verobse: print(1 - omega_pi**2/omega_cpx**2+true_solution_short_1+true_solution_short_2)
if verobse: print(1 - omega_pi**2/omega_cpx**2+appr_solution_short_1+appr_solution_short_2)

##################################### CHECK WITH A 2D MAP
ky_check = ky[250]/10**0.5/u.m
kappa = ky/10**0.5
kz=kz/u.m
ome = omega_r[250,2]*u.rad/u.s
gam = gamma[250,2]*u.rad/u.s
print(ome,gam)
# ome = ome.real
# gam = gam.real
# print(ome,gam)
# gioco_omega = np.arange(0.9*omega_check,1.1*omega_check,0.0005*prt.ionPlasmaFrequency)
# gioco_gamma = np.arange(0.9*gamma_check,1.1*gamma_check,0.001*prt.ionPlasmaFrequency)
gioco_omega = np.arange(0.1*ome/(u.rad/u.s),1.5*ome/(u.rad/u.s),ome/(u.rad/u.s)/50)*u.rad/u.s
gioco_gamma = np.arange(0.1*gam/(u.rad/u.s),1.5*gam/(u.rad/u.s),gam/(u.rad/u.s)/50)*u.rad/u.s

kyons = [ky_check*prt.Debye_length-0.0001,ky_check*prt.Debye_length+0.0001]/prt.Debye_length

solution_real = np.zeros((len(gioco_omega),len(gioco_gamma)))
solution_imag = np.zeros((len(gioco_omega),len(gioco_gamma)))
plasmaEps = partial(eps_MTSI_unnorm, prt=prt) #assign to the function eps_MTSI the value of prt from now on
max_pos = np.zeros((len(kyons),2))
for kk,kaps in enumerate(kyons):
    for i,omega_1 in enumerate(gioco_omega) :
        for j,gamma_1 in enumerate(gioco_gamma) :

            zia = 1/plasmaEps(omg=omega_1+1j*gamma_1,kx=0.0,kz=kz[0],ky=kaps)

            solution_real[i,j] = zia.real
            solution_imag[i,j] = zia.imag

    abs_sol=abs(solution_real+1j*solution_imag)
    max_pos[kk,:] = np.unravel_index(abs(abs_sol).argmax(), abs_sol.shape)
    break
print("min_calc = ", 1/abs(abs_sol).argmax())
# print("dist_ome ", ome/prt.ionPlasmaFrequency,gioco_omega[int(max_pos[0,0])]/prt.ionPlasmaFrequency)
# print("dist_ome ", ome/prt.ionPlasmaFrequency-gioco_omega[int(max_pos[0,0])]/prt.ionPlasmaFrequency)
# print("dist_gam ", gam/prt.ionPlasmaFrequency,gioco_gamma[int(max_pos[0,1])]/prt.ionPlasmaFrequency)
# print("dist_gam ", gam/prt.ionPlasmaFrequency-gioco_gamma[int(max_pos[0,1])]/prt.ionPlasmaFrequency)

current = os.getcwd()

plt.figure(figsize=(8,8))

plt.subplot(2,2,(1,2))
plt.title("Inverse of permittivity and calculated dispersion")
if unnorm == False:
    plt.pcolor(gioco_gamma.real*u.s/u.rad,gioco_omega.real*u.s/u.rad, abs(solution_real+1j*solution_imag),cmap='Reds')
else:
    plt.pcolor(gioco_gamma.real*u.s/u.rad,gioco_omega.real*u.s/u.rad, abs(solution_real+1j*solution_imag),cmap='Reds')
plt.plot(gam,ome,'o',color='blue',label='solver solution')
plt.xlabel("$\gamma/\omega_{pi}$")
plt.ylabel("$\omega/\omega_{pi}$")
if unnorm == False:
    plt.xlabel("$\gamma$ rad/s")
    plt.ylabel("$\omega$ rad/s")
    plt.text(gioco_gamma[5]*u.s/u.rad,gioco_omega[5]*u.s/u.rad,"n: {:}, ".format(plasmaDensity) + "$k_{\\theta}$: "+"{:.0f}, $k_r$: {:.0f}".format(ky*u.m,kz*u.m))
    plt.text(gioco_gamma[5]*u.s/u.rad,gioco_omega[-5]*u.s/u.rad,"(a)")
plt.colorbar().set_label("$1 / \epsilon$")
plt.legend()
plt.tight_layout()

# plt.figure(figsize=(6,5))
plt.subplot(2,2,3)
# plt.title("n: {:}, $k_r$: {:.0f}".format(plasmaDensity,kz))
plt.plot(kappa,abs(gamma[:,2]), label="solver dispersion",color='blue')
plt.plot(ky_check, gam, 'o',color='blue',label='solver solution')
plt.plot(kyons[0],gioco_gamma[int(max_pos[0,1])],'*',color='red',label = "computed solution")
plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
plt.ylabel("Growth rate  $\\gamma/\\omega_{pi}$ ")
plt.grid(True)

if unnorm== False:
    plt.xlabel("Azimuthal wave number $k_{\\theta}$ 1/m")
    plt.ylabel("Growth rate  $\\gamma$ rad/s")
    plt.text(kappa[5]*u.m,gam*u.s/u.rad,"(b)")
plt.legend()
plt.tight_layout()


# plt.figure(figsize=(6,5))
plt.subplot(2,2,4)
plt.plot(kappa,omega_r[:,2], label="solver dispersion",color='blue')
# plt.title("n: {:}, $k_r$: {:.0f}".format(plasmaDensity,kz))
plt.plot(ky_check, ome, 'o',color='blue',label="solver solution")
plt.plot(kyons[0],gioco_omega[int(max_pos[0,0])],'*',color='red',label = "computed solution")
# plt.plot(kyons[1],gioco_omega[int(max_pos[1,0])],'*',color='blue')
# plt.plot(kyons[2],gioco_omega[int(max_pos[2,0])],'*',color='blue')
plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
plt.ylabel("Pulsations  $\\omega/\\omega_{pi}$ ")
plt.grid(True)
# plt.show()
if unnorm == False:
    plt.xlabel("Azimuthal wave number $k_{\\theta}$ 1/m")
    plt.ylabel("Pulsations  $\\omega$ rad/s")
    plt.text(kappa[5]*u.m,np.amax(omega)*u.s/u.rad,"(c)")
plt.legend()
plt.tight_layout()


#################################################################################


kymin = 0.001
kymax = 0.22
pas = 0.0002383025027203481
kappa = np.arange(0.001,0.2200,0.0002383025027203481)
Nkys = len(kappa)
print("~~~~~~~~~~~~~~~~~~~~~~~Nkys",Nkys)

densities = [5e16]

kappa_solver = np.ones((len(densities),Nkys))
gamma_solver = np.ones((len(densities),Nkys))
omega_solver = np.ones((len(densities),Nkys))

for index,den in enumerate(densities):
    prtd=PlasmaParameters(plasmaDensity=den*u.m**(-3),
                        electronTemperature=10*u.eV,
                        magneticField=0.02*u.T,
                        electricField=1e4*u.V/u.m,
                        ionTemperature=0.5*u.eV)
    kz_z = kz[0]*prtd.Debye_length
    print("kz_z",kz_z)
    kappa_solver[index,:], gamma_solver[index,:], omega_solver[index,:] = verification_dispersion(kz_z, density=den,unnorm=True)
# plt.close()

plt.figure()
plt.plot(kappa_solver[0,:],gamma_solver[0,:])
plt.plot(ky/3.20,abs(gamma[:,2]),label='gamma_r')
plt.plot(kappa_solver[0,:],omega_solver[0,:])
plt.plot(ky/3.20,abs(omega_r[:,2]),label='omega_r')

# plt.plot(ky,abs(gamma[:,2]),label='gamma_r')
# plt.plot(ky,abs(omega_r[:,2]),label='omega_r')


plt.plot(kappa_solver[0,:],(gamma_solver[0,:]**2+omega_solver[0,:]**2)**0.5,"--")
plt.plot(ky/3.20,(abs(omega_r[:,2]**2)+abs(gamma[:,2]**2))**0.5,label='sum')


plt.legend()
# plt.xlim(0,1500)

# use the solution from the solver inside the real equation and the approximated equation
omega_cpx = np.zeros(len(kz))
omega_cpx[:len(kz)] = gamma_solver[0,:len(kz)]*1j+omega_solver[0,:len(kz)]
kz = kz*u.m
true_solution1 = 1 - omega_pi**2/omega_cpx**2 - (kz**2 * omega_pe**2)/(omega_cpx-ky * v0)**2/k**2 - (ky**2 * omega_pe**2)/((omega_cpx-ky * v0)**2 -omega_ce**2)/k**2
appr_solution1 = (alpha*omega_cpx**2 + beta * omega_cpx**3 - omega_pi**2)/omega_cpx**2

plt.figure()
plt.semilogy(ky,abs(true_solution1.real),color='r')
plt.semilogy(ky,abs(true_solution1.imag),color='b')
plt.show()

print(true_solution)
# print(appr_solution)
# plt.show()



# print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
# eps = kz**2*omega_pe**2/(k**2 * ky**2 * v0**2)
# for ind in range(len(kz)):
#     a = 2 * eps[ind]/(ky*v0)
#     b = 1 - eps[ind]
#     c = 0*u.rad
#     d = omega_pi**2
#
#     # a = 1*u.rad
#     # b = -3*u.rad
#     # c = 3*u.rad
#     # d = -1*u.rad
#     print(a.value,b.value,c.value,d.value)
#     omega_tot[ind,:] = np.roots([a.value,b.value,c.value,d.value])
#     print(omega_tot)
#
# print(a.value*omega_tot[0,0]**3 + b.value*omega_tot[0,0]**2 + c.value*omega_tot[0,0] + d.value)
# print(a.value*omega_tot[0,1]**3 + b.value*omega_tot[0,1]**2 + c.value*omega_tot[0,1] + d.value)
# print(a.value*omega_tot[0,2]**3 + b.value*omega_tot[0,2]**2 + c.value*omega_tot[0,2] + d.value)
#
# omega_cpx = omega_tot[0,0]
# true_solution = 1 - omega_pi**2/omega_cpx**2 - (kz**2 * omega_pe**2)/(omega_cpx-ky * v0)**2/k**2
# print("omega_tot {:e} {:e}".format(omega_tot[0,0].real,omega_tot[0,0].imag))
# print(true_solution)
#
# omega_cpx = omega_tot[0,1]
# true_solution = 1 - omega_pi**2/omega_cpx**2 - (kz**2 * omega_pe**2)/(omega_cpx-ky * v0)**2/k**2 - (ky**2 * omega_pe**2)/((omega_cpx-ky * v0)**2 -omega_ce**2)/k**2
# print("omega_tot {:e} {:e}".format(omega_tot[0,1].real,omega_tot[0,1].imag))
# print(true_solution)
#
# omega_cpx = omega_tot[0,2]
# true_solution = 1 - omega_pi**2/omega_cpx**2 - (kz**2 * omega_pe**2)/(omega_cpx-ky * v0)**2/k**2 - (ky**2 * omega_pe**2)/((omega_cpx-ky * v0)**2 -omega_ce**2)/k**2
# print("omega_tot {:e} {:e}".format(omega_tot[0,2].real,omega_tot[0,2].imag))
# print(true_solution)
#
# print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
# omega_cpx = 3e6-2.7j*1e6
# L_r = 0.0128
# L_t = 0.0128/3
# kz = 2*np.pi/L_r/2
# ky = 2*np.pi/L_t
# ky = 245
# k = (kz**2 + ky**2)**0.5
# # true_solution = 1 - omega_pi**2/omega_cpx**2 - (kz**2 * omega_pe**2)/(omega_cpx-ky * v0)**2/k**2 - (ky**2 * omega_pe**2)/((omega_cpx-ky * v0)**2 -omega_ce**2)/k**2
# uno = 1 - omega_pi**2/omega_cpx**2
# due = - (kz**2 * omega_pe**2)/(omega_cpx-ky * v0)**2/k**2
# due_approx = -(kz**2 * omega_pe**2)/(ky**2 * v0**2)*(1+2*omega_cpx/(ky*v0))/k**2
# tre = - (ky**2 * omega_pe**2)/((omega_cpx-ky * v0)**2 -omega_ce**2)/k**2
# tre_approx = (ky**2 * omega_pe**2)/(omega_ce**2 * k**2) * (1+ky**2*v0**2/omega_ce**2 - 2*ky*v0*omega_cpx/omega_ce**2)
# print("omega",omega_cpx,"kz",kz,"ky",ky)
# print("(1): {:e} {:e}j".format(uno.real,uno.imag))
# print("(2) true: ",(due))
# print("(2) appr: ",(due_approx))
# print("(3) true: ",(tre))
# print("(3) appr: ",(tre_approx))
#
#
# print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
# alpaca = (kz*omega_pe/(k*ky*v0))**2
# a = 1
# b = ky*v0/(2*alpaca)*(1-alpaca)
# c = 0
# d = ky*v0/(2*alpaca)*omega_pi**2
# print(alpaca,a,b,c,d)
# omega_tot[0,:] = np.roots([a,b.value,c,d.value])
# print(omega_tot)
# omega_tot[0,:] = cubic_solver.cubic_solver(a,b.value,c,d.value)
# print(omega_tot)
#
#
# print(a*omega_tot[0,0]**3 + b.value*omega_tot[0,0]**2 + c*omega_tot[0,0] + d.value)
# print(a*omega_tot[0,1]**3 + b.value*omega_tot[0,1]**2 + c*omega_tot[0,1] + d.value)
# print(a*omega_tot[0,2]**3 + b.value*omega_tot[0,2]**2 + c*omega_tot[0,2] + d.value)
#
#
#
#
# # print(np.roots([1,-3,12,1]))
# # print(cubic_solver.cubic_solver(1,-3,12,1))
#
#
#
#
#
#
# # plt.plot(kz*prt.Debye_length/u.m,gamma)
