import numpy as np
from util.MTSI  import eps_MTSI
from util.MTSI_unnorm  import eps_MTSI_unnorm
from functools import partial
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "STIXGeneral"
import os




def open_disp_file(kz, path = None):
    if path == None:
        path = '/home/petronio/Nextcloud/theseLPP/runs/runs_benchmark/MTSI/dispersion_MTSI/dispersion_solver/dispersion_data/general_results/'
    kymin = 0.001
    kymax = 0.20
    pas = 0.00023803827751196175
    Nkys = (kymax-kymin)/pas
    Nkys = int(Nkys)
    kappa = np.arange(kymin,kymax,pas)

    omega = np.genfromtxt(path + "kz={:5.4f}".format(kz) + "_omega_r.txt", delimiter="  ", unpack=False)
    gamma = np.genfromtxt(path + "kz={:5.4f}".format(kz) + "_gamma.txt", delimiter="  ", unpack=False)
    return kappa,omega,gamma

def precedent_openfile(kz,Nkys=920,path=None):
    if path == None:
        path = '/home/petronio/Nextcloud/theseLPP/runs/runs_benchmark/MTSI/dispersion_MTSI/dispersion_solver/dispersion_data/general_results/'

    # kappa = np.genfromtxt(path + "ky.txt", delimiter="  ")
    # if kz < 0.0099:
    #     start = 0.001
    #     stop = 0.1
    #     steps = 500
    #     kappa = np.arange(start,stop,(stop-start)/steps)
    # else:
    print("kz_preopen : {:.4f}".format(kz) )
    loop = 1
    while True:
        try :
            omega_read = np.genfromtxt(path + "kz={:5.4f}".format(kz) + "_omega_r.txt", delimiter="  ", unpack=False)
            print("kz_open : {:.4f}".format(kz) )
            break
        except :
            print(path + "kz={:5.4f}".format(kz) + "_omega_r.txt")
            kz = kz + 0.0001*loop*(-1)**loop
            loop = loop+1
            # break

    gamma_read = np.genfromtxt(path + "kz={:5.4f}".format(kz) + "_gamma.txt", delimiter="  ", unpack=False)

    # print("Len",len(omega_read))

    if kz > 0.0516 :
        Nkys = 1844

    omega = np.ones(Nkys)*1e-12
    gamma = np.ones(Nkys)*1e-12

    omega[:len(omega_read)] = omega_read
    gamma[:len(gamma_read)] = gamma_read

    return omega, gamma, kz


def find_max_gamma(kz,path=None):
    if path == None:
        path = '/home/petronio/Nextcloud/theseLPP/runs/runs_benchmark/MTSI/dispersion_MTSI/dispersion_solver/dispersion_data/general_results/'

    omega,gamma,kz = precedent_openfile(kz,path=path)
    max_ind = np.argmax(gamma)

    start = 0.001
    stop = 0.22
    pas = 0.0002383025027203481

    kapa = np.arange(start,stop,pas)
    print(kz,kapa[max_ind],omega[max_ind],gamma[max_ind])

    return kapa[max_ind],omega[max_ind],gamma[max_ind]




def verification_dispersion(kz,density=5e16,unnorm = False,EF=1e4):
    from util.parameters import PlasmaParameters
    from astropy.constants import m_e, m_p
    from astropy import units as u

    current = os.getcwd()
    path = current + "/dispersion_data/change_n/{:}/".format(density)
    # path = current + "/dispersion_data/change_n_E/20000.0_1e+17/"


    Te = 10*u.eV
    plasmaDensity=density*u.m**(-3)
    EF = EF*u.V/u.m
    prt=PlasmaParameters(plasmaDensity=plasmaDensity,
                        electronTemperature=10*u.eV,
                        magneticField=0.02*u.T,
                        electricField=EF,
                        ionTemperature=0.5*u.eV)

    kx = 0.0

    kappa = np.arange(0.001,0.2200,0.0002383025027203481)
    # kymin = 0.001
    # kymax = 0.22001
    # pas = 0.0002383025027203481
    # if kz > 0.0516 :
    #     kymax = 0.44001
    #     pas = 0.0002381025027203481
    # kappa = np.arange(kymin,kymax,pas)

    omega,gamma,kz = precedent_openfile(kz,path=path)
    omega = omega[:len(kappa)]
    gamma = gamma[:len(kappa)]

    # print("kz_opened : {:.4f}".format(kz) )


    ky,ome,gam = find_max_gamma(kz,path=path)

    ind_ky = list(kappa).index(ky)
    # ind_ky = int(ind_ky*1.5)
    ky = kappa[ind_ky]
    ome = omega[ind_ky]
    gam = gamma[ind_ky]


    if unnorm:
        kappa = kappa/prt.Debye_length
        omega = omega*prt.ionPlasmaFrequency
        gamma = gamma*prt.ionPlasmaFrequency

    # plt.figure()
    # plt.plot(kappa,gamma)
    # plt.plot(ky/prt.Debye_length,gam*prt.ionPlasmaFrequency,'o')
    print("ky = ",ky,ky/prt.Debye_length)
    gioco_omega = np.arange(0.9*ome,1.1*ome,0.0005)
    gioco_gamma = np.arange(0.9*gam,1.1*gam,0.001)
    kyons = [ky-0.0001,ky+0.0001]

    if unnorm:
        ky = ky/prt.Debye_length
        kz = kz/prt.Debye_length
        gam = gam*prt.ionPlasmaFrequency
        ome = ome*prt.ionPlasmaFrequency

        # ome = 13920507.92438155 *u.rad / u.s
        # gam = 23312202.01594126 *u.rad / u.s
        # ky =  1793.4063823828753 / u.m
        # kz =  490.87385212340513 / u.m

        gioco_omega = np.arange(0.5*ome/(u.rad/u.s),1.5*ome/(u.rad/u.s),ome/(u.rad/u.s)/100)*u.rad/u.s
        gioco_gamma = np.arange(0.5*gam/(u.rad/u.s),1.5*gam/(u.rad/u.s),gam/(u.rad/u.s)/100)*u.rad/u.s
        kyons = kyons/prt.Debye_length
        # print(ky)
        # print("gioco_omega",gioco_omega,"\n gioco_gamma",gioco_gamma)
        # print("ome", ome, gioco_omega[int(len(gioco_omega)/2)])
        # print("gam", gam, gioco_gamma[int(len(gioco_gamma)/2)])
    print("density: ", plasmaDensity)
    print("max_ome: ", ome)
    print("max_gam: ", gam)
    print("ky = ",ky)
    print("kz = ",kz)

    solution_real = np.zeros((len(gioco_omega),len(gioco_gamma)))
    solution_imag = np.zeros((len(gioco_omega),len(gioco_gamma)))
    if unnorm:
        plasmaEps = partial(eps_MTSI_unnorm, prt=prt) #assign to the function eps_MTSI the value of prt from now on
    else:
        plasmaEps = partial(eps_MTSI, prt=prt) #assign to the function eps_MTSI the value of prt from now on

    max_pos = np.zeros((len(kyons),2))
    for kk,kaps in enumerate(kyons):
        for i,omega_1 in enumerate(gioco_omega) :
            for j,gamma_1 in enumerate(gioco_gamma) :

                zia = 1/plasmaEps(omg=omega_1+1j*gamma_1,kx=0.0,kz=kz,ky=kaps)

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
    if unnorm:
        plt.pcolor(gioco_gamma*u.s/u.rad,gioco_omega*u.s/u.rad, abs(solution_real+1j*solution_imag),cmap='Reds')
    else:
        plt.pcolor(gioco_gamma,gioco_omega, abs(solution_real+1j*solution_imag),cmap='Reds')
    plt.plot(gam,ome,'o',color='blue',label='solver solution')
    plt.xlabel("$\gamma/\omega_{pi}$")
    plt.ylabel("$\omega/\omega_{pi}$")
    if unnorm:
        plt.xlabel("$\gamma$ rad/s")
        plt.ylabel("$\omega$ rad/s")
        plt.text(gioco_gamma[5]*u.s/u.rad,gioco_omega[5]*u.s/u.rad,"n: {:}, ".format(plasmaDensity) + "$k_{\\theta}$: "+"{:.0f}, $k_r$: {:.0f}".format(ky,kz))
        plt.text(gioco_gamma[5]*u.s/u.rad,gioco_omega[-5]*u.s/u.rad,"(a)")
    plt.colorbar().set_label("$1 / \epsilon$")
    plt.legend()
    plt.tight_layout()

    # plt.figure(figsize=(6,5))
    plt.subplot(2,2,3)
    # plt.title("n: {:}, $k_r$: {:.0f}".format(plasmaDensity,kz))
    plt.plot(kappa,abs(gamma), label="solver dispersion",color='blue')
    plt.plot(ky, gam, 'o',color='blue',label='solver solution')
    plt.plot(kyons[0],gioco_gamma[int(max_pos[0,1])],'*',color='red',label = "computed solution")
    plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
    plt.ylabel("Growth rate  $\\gamma/\\omega_{pi}$ ")
    plt.grid(True)

    if unnorm:
        plt.xlabel("Azimuthal wave number $k_{\\theta}$ 1/m")
        plt.ylabel("Growth rate  $\\gamma$ rad/s")
        plt.text(kappa[5]*u.m,gam*u.s/u.rad,"(b)")
    plt.legend()
    plt.tight_layout()


    # plt.figure(figsize=(6,5))
    plt.subplot(2,2,4)
    plt.plot(kappa,omega, label="solver dispersion",color='blue')
    # plt.title("n: {:}, $k_r$: {:.0f}".format(plasmaDensity,kz))
    plt.plot(ky, ome, 'o',color='blue',label="solver solution")
    plt.plot(kyons[0],gioco_omega[int(max_pos[0,0])],'*',color='red',label = "computed solution")
    # plt.plot(kyons[1],gioco_omega[int(max_pos[1,0])],'*',color='blue')
    # plt.plot(kyons[2],gioco_omega[int(max_pos[2,0])],'*',color='blue')
    plt.xlabel("Azimuthal wave number $k_{\\theta} \\lambda_{De}$")
    plt.ylabel("Pulsations  $\\omega/\\omega_{pi}$ ")
    plt.grid(True)

    if unnorm:
        plt.xlabel("Azimuthal wave number $k_{\\theta}$ 1/m")
        plt.ylabel("Pulsations  $\\omega$ rad/s")
        plt.text(kappa[5]*u.m,np.amax(omega)*u.s/u.rad,"(c)")
    plt.legend()
    plt.tight_layout()
    # plt.savefig(current + "/images_dispersion/" + "solution_verif{:}.png".format(density))
    # zia = eps_MTSI_unnorm(omg=ome+1j*gam, kx=0.0, ky=1793.4/u.m, kz=kz, prt=prt,impr=True)
    # plt.close('all')
    # plt.show()
    return kappa,gamma,omega
