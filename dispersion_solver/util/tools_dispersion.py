import numpy as np
from util.MTSI  import eps_MTSI
from functools import partial
import matplotlib.pyplot as plt





def open_disp_file(kz, path = None):
    if path == None:
        path = '/home/petronio/Nextcloud/theseLPP/runs/runs_benchmark/MTSI/dispersion_MTSI/dispersion_solver/dispersion_data/general_results/'
    kappa = np.genfromtxt(path + "ky.txt", delimiter="  ")
    omega = np.genfromtxt(path + "kz={:5.4f}".format(kz) + "_omega_r.txt", delimiter="  ", unpack=False)
    gamma = np.genfromtxt(path + "kz={:5.4f}".format(kz) + "_gamma.txt", delimiter="  ", unpack=False)
    return kappa,omega,gamma


def find_max_gamma(kz,path=None):
    kappa,omega,gamma = open_disp_file(kz,path)
    max_ind = np.argmax(gamma)
    # if kz < 0.0099:
    #     start = 0.001
    #     stop = 0.1
    #     steps = 500
    #     kapa = np.arange(start,stop,(stop-start)/steps)
    # else:
    start = 0.001
    stop = 0.22
    pas = 0.00023803827751196175

    kapa = np.arange(start,stop,pas)
    print(kz,kapa[max_ind],omega[max_ind],gamma[max_ind])

    return kapa[max_ind],omega[max_ind],gamma[max_ind]

def verification_dispersion(kz):
    from util.parameters import PlasmaParameters
    from astropy.constants import m_e, m_p
    from astropy import units as u

    path = "/home/petronio/Nextcloud/theseLPP/runs/runs_benchmark/MTSI/dispersion_MTSI/dispersion_solver/dispersion_data/change_n/5e+16/"
    Te = 10*u.eV
    plasmaDensity=5e16 *u.m**(-3)
    pp = PlasmaParameters(plasmaDensity=plasmaDensity, electronTemperature=Te)


    prt=PlasmaParameters(plasmaDensity=plasmaDensity,
                        electronTemperature=10*u.eV,
                        magneticField=0.02*u.T,
                        electricField=1e4*u.V/u.m,
                        ionTemperature=0.5*u.eV)

    kx = 0.0
    kappa,omega,gamma = open_disp_file(kz,path=path)
    ky,ome,gam = find_max_gamma(kz,path=path)

    gioco_omega = np.arange(0.8*ome,1.2*ome,0.005)
    gioco_gamma = np.arange(0.8*gam,1.2*gam,0.005)

    solution_real = np.zeros((len(gioco_omega),len(gioco_gamma)))
    solution_imag = np.zeros((len(gioco_omega),len(gioco_gamma)))
    plasmaEps = partial(eps_MTSI, prt=prt) #assign to the function eps_MTSI the value of prt from now on

    kyons = [ky-0.0001,ky,ky+0.0001]
    # kyons = ky
    max_pos = np.zeros((3,2))
    for kk,kaps in enumerate(kyons):
        for i,omega_1 in enumerate(gioco_omega) :
            for j,gamma_1 in enumerate(gioco_gamma) :
                # solution[i,j] = plasmaEps(omg=omega+1j*gamma,kx=0.0,kz=kz,ky=ky)
                zia = 1/plasmaEps(omg=omega_1+1j*gamma_1,kx=0.0,kz=kz,ky=kaps)
                # print(zia)
                solution_real[i,j] = zia.real
                solution_imag[i,j] = zia.imag

        abs_sol=abs(solution_real+1j*solution_imag)
        max_pos[kk,:] = np.unravel_index(abs(abs_sol).argmax(), abs_sol.shape)

    # print(gioco_gamma[int(max_pos[0,1])])
    plt.figure()
    plt.title("invers of susceptibility ")
    plt.pcolor(gioco_gamma,gioco_omega, abs(solution_real+1j*solution_imag))
    plt.xlabel("$\gamma/\omega_{pi}$")
    plt.ylabel("$\omega/\omega_{pi}$")
    # plt.text(x=gioco_gamma[-70],y=gioco_omega[-25],s="kz = %5.4f \n"%kz + "ky = %5.4f \n"%ky+
    #         "$\omega_{max}$ = %6.4f \n"%gioco_omega[max_pos[0]] + "$\gamma_{max}$ = %6.4f"%gioco_gamma[max_pos[1]],color='red')
    plt.colorbar()
    #
    plt.figure(figsize=(6,5))
    plt.plot(kappa,gamma, label="solver solution")
    plt.plot(kyons[0],gioco_gamma[int(max_pos[0,1])],'*',color='blue',label = "computed solution")
    plt.plot(kyons[1],gioco_gamma[int(max_pos[1,1])],'*',color='blue')
    plt.plot(kyons[2],gioco_gamma[int(max_pos[2,1])],'*',color='blue')
    plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
    plt.ylabel("Growth rate  $\\gamma/\\omega_{pi}$ ")
    plt.legend()


    plt.figure(figsize=(6,5))
    plt.plot(kappa,omega, label="solver solution")
    plt.plot(kyons[0],gioco_omega[int(max_pos[0,0])],'*',color='blue',label = "computed solution")
    plt.plot(kyons[1],gioco_omega[int(max_pos[1,0])],'*',color='blue')
    plt.plot(kyons[2],gioco_omega[int(max_pos[2,0])],'*',color='blue')
    plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
    plt.ylabel("Pulsations  $\\omega/\\omega_{pi}$ ")
    plt.legend()
    plt.show()
    return max_pos
