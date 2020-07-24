from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import matplotlib.pyplot as plt
from util.iaw import precedent_openfile, precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion
from util.parameters import PlasmaParameters
from astropy.constants import m_e, m_p
from astropy import units as u


kz = 0.0450
density = 5e16
# max_pos = verification_dispersion(kz, density=density,unnorm=True)


prt_base=PlasmaParameters(plasmaDensity=density*u.m**(-3),
                    electronTemperature=10*u.eV,
                    magneticField=0.02*u.T,
                    electricField=1e4*u.V/u.m,
                    ionTemperature=0.5*u.eV)

densities = [5e16,1e17,2e17,3e17]
for dindondensity in densities:
    prtd=PlasmaParameters(plasmaDensity=dindondensity*u.m**(-3),
                        electronTemperature=10*u.eV,
                        magneticField=0.02*u.T,
                        electricField=1e4*u.V/u.m,
                        ionTemperature=0.5*u.eV)

    kz_z = kz/prt_base.Debye_length*prtd.Debye_length
    print(kz_z)
    max_pos = verification_dispersion(kz_z, density=dindondensity,unnorm=True)
plt.show()
