from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import matplotlib.pyplot as plt
from util.iaw import  precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion,precedent_openfile

kz = 0.030
densities = [5e16,2e17]
for den in densities:
    kappa, gamma, omega = verification_dispersion(kz, density=den,unnorm=False)

plt.show()
