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
# ky = np.arange(1,8000,1)
ky = np.arange(0.001,0.8200,0.0002383025027203481)/prt.Debye_length*u.m
kz = 0.019*np.ones(len(ky))/prt.Debye_length*u.m
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
omega_sol_compls = np.zeros((len(kz),3), dtype=complex)
evaluation = np.zeros((len(kz),3), dtype=complex)
evaluation2 = np.zeros((len(kz),3), dtype=complex)
evaluation3 = np.zeros((len(kz),3), dtype=complex)
evaluation4 = np.zeros((len(kz),3), dtype=complex)
evaluation5 = np.zeros((len(kz),3), dtype=complex)
evaluation6 = np.zeros((len(kz),3), dtype=complex)


# print(ky, 2*np.pi/(2*0.0128))

omega_ce = prt.electronCyclotronFrequency*u.s
omega_pe = prt.electronPlasmaFrequency*u.s/u.rad
omega_pi = prt.ionPlasmaFrequency*u.s/u.rad
lambda_D = prt.Debye_length/u.m

##### normalize

kz = kz*lambda_D
ky = ky*lambda_D
k = k*lambda_D
omega_pe = omega_pe/omega_pi
omega_ce = omega_ce/omega_pi
mass_ratio = 1836.152*131
v0 = 5e4/(lambda_D*omega_pi)
print(mass_ratio)

# print(kz,ky,k,v0,omega_ce,omega_pi,omega_pe)
# alpha = (1*u.rad**2 - (kz**2 * omega_pe**2)/(k**2 * ky**2 * v0**2) + (ky**2 * omega_pe**2)/(k**2 * omega_ce**2) * (1+ky**2 * v0**2 / omega_ce**2))

# alpha = (1 - (kz**2 * omega_pe**2)/(k**2 * ky**2 * v0**2) + (ky**2 * omega_pe**2)/(k**2 * omega_ce**2) * (1+ky**2 * v0**2 / omega_ce**2))
# beta = -(2 * (kz**2 * omega_pe**2)/(k**2 * ky**3 * v0**3) + 2 * (ky**3 * omega_pe**2 * v0)/(k**2 * omega_ce**4))

alpha = 1 - mass_ratio * kz**2 / (k**2 * ky**2 * v0**2) + mass_ratio * ky**2 / (k**2 * omega_ce**2) * (1 + ky**2 * v0**2 / omega_ce**2)
beta = - mass_ratio * kz **2 * 2 / (k**2 * ky**3 * v0**3) -2 * mass_ratio * ky**3 * v0 / (omega_ce**4 * k**2)

# print("alpha", alpha[0::10])
# print("beta", beta[0::10])

if verobse: print("From equation (A): alpha = {:e}, beta = {:e}".format(alpha[10],beta[10]))

# alpha = alpha.value
# beta = beta.value

for ind in range(len(kz)):
    # omega_sol_compls[ind,:] = np.roots([alpha[ind].value,beta[ind].value,0,-1])
    p = np.poly1d([alpha[ind].value,beta[ind].value,0,-1])
    omega_sol_compls[ind,:] = p.r
    evaluation[ind,:] = p(p.r)
    # evaluation2[ind,:] = alpha[ind].value*p.r**3 + beta[ind].value*p.r**2 + 0 - 1
    evaluation2[ind,:] = 1 - 1/p.r**2 - mass_ratio * kz[ind]**2/((p.r - ky[ind]*v0)**2 * k[ind]**2) - mass_ratio * ky[ind]**2/(((p.r-ky[ind]*v0)**2-omega_ce**2)*k[ind]**2)
    evaluation3[ind,:] = - mass_ratio * kz[ind]**2/((p.r - ky[ind]*v0)**2 * k[ind]**2)
    evaluation4[ind,:] = - mass_ratio * kz[ind]**2/(ky[ind]**2 * v0**2 * k[ind]**2) * (1+2*p.r/(ky[ind] * v0))
    # evaluation3[ind,:] = - 1/((p.r - ky[ind]*v0)**2)
    # evaluation4[ind,:] = - 1/(ky[ind]**2 * v0**2) * (1+2*p.r/(ky[ind] * v0))
    evaluation5[ind,:] = - mass_ratio * ky[ind]**2/(((p.r-ky[ind]*v0)**2-omega_ce**2)*k[ind]**2)
    evaluation6[ind,:] = mass_ratio * ky[ind]**2/(omega_ce**2*k[ind]**2) * (1 + ky[ind]**2 * v0**2 /omega_ce**2 - 2 *p.r*ky[ind]*v0/omega_ce**2)

# print(evolution2[])

plt.figure()
# plt.plot(ky,omega_sol_compls[:,0].real,label="0")
# plt.plot(ky,omega_sol_compls[:,0].imag,label="0")
plt.plot(ky,omega_sol_compls[:,1].real,label="1 real")
plt.plot(ky,omega_sol_compls[:,1].imag,label="1 imag")
plt.plot(ky,omega_sol_compls[:,2].real,"--",label="2 real")
plt.plot(ky,omega_sol_compls[:,2].imag,"--",label="2 imga")
plt.legend()

plt.figure()
for ind in range(len(kz)):
    a = 8 * beta[ind]**2
    b = 8 * alpha[ind] * beta[ind]
    c = 2 * alpha[ind]**2
    d = beta[ind]
    # if verobse: print("coefficients of equation (C): {:e},{:e},{:e},{:e}".format(a.value,b.value,c.value,d.value))
    # omega_r[ind,:] = np.roots([a.value,b.value,c.value,d.value])
    p = np.poly1d([a.value,b.value,c.value,d.value])
    omega_r[ind,:] = p.r
    gamma[ind,:] = ((-1**2 + omega_r[ind]**2 * alpha[ind] + omega_r[ind]**3 * beta[ind])/(alpha[ind] + 3 * omega_r[ind] * beta[ind]))**0.5
    pp = np.poly1d([alpha[ind].value,beta[ind].value,0,-1])
    soliz = pp.r
    # evaluation[ind,:] = p(soliz.real)
    # evaluation2[ind,:] = alpha[ind].value*soliz.real**3 + beta[ind].value*soliz.real**2 + 0 - 1





# print("unit",a.unit,b.unit,c.unit,d.unit)
# print("value",a.value,b.value,c.value,d.value)

omega_cpx = omega_r[:,1] + 1j*gamma[:,1]
if verobse: print("omega = {:e} + j{:e}".format(omega_cpx[10].real,omega_cpx[10].imag))
if verobse: print("|omega| = {:e}".format((omega_cpx[10].real**2+omega_cpx[10].imag**2)**0.5))
if verobse: print("(C) solved function for omega: ",a.value*omega_r[:,2]**3 + b.value*omega_r[:,2]**2 + c.value*omega_r[:,2] + d.value)
if verobse: print("(A) solved function for omega/gamma: ",alpha*omega_cpx**2 + beta * omega_cpx**3 - omega_pi**2)

# print(omega_r)
# print(gamma)

# plt.plot(kz,abs(gamma[:,0]),label='0')
plt.title("kz = {:.4f}".format(kz[0]))

# plt.plot(ky,abs(gamma[:,0])/omega_pi,label='0')
# plt.plot(ky,abs(omega_r[:,0]/omega_pi),label='0')
# plt.plot(ky,abs(gamma[:,1])/omega_pi,label='1')
# plt.plot(ky,abs(omega_r[:,1]/omega_pi),label='1')
ky = ky*3.2
# plt.plot(ky/3.2,(gamma[:,2]),label='gamma old')
# plt.plot(ky/3.2,(omega_r[:,2]),label='omega old')
# plt.plot(ky/3.2,(abs(omega_r[:,2]**2)+abs(gamma[:,2]**2))**0.5,label='somma old')

plt.semilogy(ky/3.2,(omega_sol_compls[:,1].real),label='omega')
plt.semilogy(ky/3.2,(omega_sol_compls[:,1].imag),label='gamma')
plt.semilogy(ky/3.2,((omega_sol_compls[:,1].real)**2+(omega_sol_compls[:,2].imag)**2)**0.5,label='somma')
ky = ky/3.2
plt.semilogy(ky,ky*v0)
plt.ylim(1e-4,10)




plt.legend()

current = os.getcwd()
path = current + "/dispersion_data/change_n/5e+16/"
kappa = np.arange(0.001,0.2200,0.0002383025027203481)
omega,gamma,kz = precedent_openfile(kz[0],path=path)
# plt.plot(kappa,omega,"--")
# plt.plot(kappa,gamma,"--")
# plt.plot(kappa,omega+gamma,"--")



plt.figure()
omega_r[:,2]/3.2
# plt.semilogy(ky,abs(beta.value*omega_sol_compls[:,2]**3 + alpha.value*omega_sol_compls[:,2]**2 + 0*omega_sol_compls[:,2] -1),"v",label = "complex solver")
# plt.semilogy(ky,abs(beta.value*omega_sol_compls[:,1]**3 + alpha.value*omega_sol_compls[:,1]**2 + 0*omega_sol_compls[:,1] -1),"v",label = "complex solver")
# plt.semilogy(ky,evaluation[:,0])
#
# plt.semilogy(ky,evaluation[:,1],'*',label="1")
# plt.semilogy(ky,evaluation[:,2],label="2")

# plt.semilogy(ky,evaluation2[:,0])
# plt.semilogy(ky,abs(evaluation2[:,1]),'*',label="1 bis")
# plt.semilogy(ky,abs(evaluation2[:,2]),label="2 bis")
plt.semilogy(ky,abs(evaluation3[:,1] - evaluation4[:,1]),label="diff 1o")
plt.semilogy(ky,abs(evaluation3[:,1]),label="evo3 1o")
plt.semilogy(ky,abs(evaluation4[:,1]),label="evo4 1o")

# plt.semilogy(ky,abs(evaluation5[:,1] - evaluation6[:,1]),"--",label="diff 2o")
# plt.semilogy(ky,abs(evaluation3[:,1] + evaluation5[:,1] - evaluation4[:,1] - evaluation6[:,1]),"--",label="diff 3o")



# plt.semilogy(abs(a.value*omega_r[:,1]**3 + b.value*omega_r[:,1]**2 + c.value*omega_r[:,1] + d.value),"^",label = "1")
# plt.semilogy(abs(a.value*omega_r[:,2]**3 + b.value*omega_r[:,2]**2 + c.value*omega_r[:,2] + d.value),"o",label = "2")
# plt.semilogy(abs(a.value*omega**3 + b.value*omega**2 + c.value*omega + d.value),"o",label = "solver")

plt.legend()
# plt.axhline(y=1)
# plt.plot(,a.value*omega_r[:,2]**3 + b.value*omega_r[:,2]**2 + c.value*omega_r[:,2] + d.value)



plt.show()
