from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import matplotlib.pyplot as plt
from util.iaw import precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion,precedent_openfile
import os


### max gamma evolution
sentierino = os.getcwd()

fig, ax1 = plt.subplots(figsize=(6,5))
ax1.set_xlabel("kz $\lambda_D$")
ax1.set_ylabel('ky $\lambda_D$ max', color='blue')
ax2 = ax1.twinx()
ax2.set_ylabel('$\gamma / \omega_{pi}$', color='red')
kzetas = np.arange(0.0011,0.055,0.0004)
for kz in kzetas:
    kap,ome,gam = find_max_gamma(kz)
    ax1.plot(kz,kap,'*',color="blue")
    ax2.plot(kz,gam,'v',color="red")
plt.grid(True)
plt.savefig(sentierino   + "/images_dispersion/max_gamma.png")
