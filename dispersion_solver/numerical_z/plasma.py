import numpy as np

class PlasmaParameters:
    """Parameters needed for the Dispersion Relation """

    def __init__(self, ne=1e17, te=50, ti=0.1, bf=0.02, ef=4e4, ni=131):
        """sets the parameters of the Plasma"""
        self.qe = 1.6e-19  # "Electric charge (C)"
        self.c = 3e8  # "Speed of light (m/s)"
        self.me = 9.1e-31  # "Electron mass (Kg)"
        self.eps0 = 8.85e-12  # "Vacuum permittivity(F/m)"
        self.u = 1.6605e-27

        self.mi = ni*self.u  # "Ion mass (Kg)"

        "Parameters of the system"
        self.ne = ne  # "Plasma density (m^-3)"
        self.te = te  # "Electron temperature (eV)"
        self.ti = ti  # "Ion temperature (eV)"
        self.bf = bf  # "Magnetic field (T)"
        self.ef = ef  # "Electric field (V/m)"

        self.opi = np.sqrt(self.ne * self.qe ** 2 / (self.eps0 * self.mi))  # "Ion plasma frequency"
        self.ld = np.sqrt(self.eps0 * self.te / (self.ne * self.qe))  # "Debye length"
        self.cs = self.ld * self.opi  # "Sound speed"
        "Derived quantities (normalized)"

        self.ope = np.sqrt(self.ne * self.qe ** 2 / (self.eps0 * self.me)) / self.opi  # "Ion plasma frequency"
        self.oce = self.qe * self.bf / (self.me * self.opi)  # "Electron cyclotron frquency"
        self.vte = np.sqrt(2 * self.qe * self.te / self.me) / self.cs  # "Electron thermal velocity"
        self.vde = self.ef / (self.bf * self.cs)  # "Electron drift velocity"
        # vde=ef/bf
        # vde = 4.5*vte

        self.vti = np.sqrt(2 * self.qe * self.ti / self.mi) / self.cs  # "Ion thermal velocity"
        self.rL = self.vte / self.oce  # "Larmor radius"
        self.eps = 0.01
        self.vtb = 0.5 * self.vte
        self.vdi = 0  # 10000 / self.cs
