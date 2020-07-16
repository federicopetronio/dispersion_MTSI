from plasmapy.formulary.dispersionfunction import plasma_dispersion_func
import numpy as np
import matplotlib.pyplot as plt
from util.iaw import precedent_openfile, precedent_guess,precedent_guess_mod
from util.tools_dispersion import open_disp_file,find_max_gamma,verification_dispersion


### max gamma evolution

fig, ax1 = plt.subplots(figsize=(6,5))
ax1.set_xlabel("kz")
ax1.set_ylabel('ky max', color='blue')
ax2 = ax1.twinx()
ax2.set_ylabel('ky max', color='red')
kzetas = np.arange(0.0011,0.006,0.0001)
for kz in kzetas:
    kap,ome,gam = find_max_gamma(kz)
    ax1.plot(kz,kap,'*',color="blue")
    ax2.plot(kz,gam,'v',color="red")
plt.show()
