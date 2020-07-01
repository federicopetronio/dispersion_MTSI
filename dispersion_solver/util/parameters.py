from astropy import units as u
from plasmapy.formulary.parameters import plasma_frequency, Debye_length, thermal_speed
from plasmapy.formulary.drifts import ExB_drift
from astropy.constants import e, m_e


class PlasmaParameters():
    """this class "pack" the plasma parameters to be used for normalisation"""

    def __init__(self,
                 plasmaDensity=1e16 *u.m**(-3),
                 ionMass=131*u.u,
                 electronTemperature=10*u.eV,
                 ionTemperature=0.1*u.eV,
                 magneticField=2e-4*u.T,
                 electricField=1e4*u.V/u.m,
                 ionName="Xe+"

                 ):
        """ all the arguments should have a dimension, except of the ion name which is a `str` """

        self.plasmaDensity = plasmaDensity
        assert self.plasmaDensity.unit == u.m**(-3), f"The density unit is not similare to 'per cubic meters' m$^{-3}$, instead you provided {self.plasmaDensity.unit}"

        self.ionMass = ionMass
        assert self.ionMass.unit.physical_type == (1*u.kg).unit.physical_type, f"The ionmass is not similare to 'kg', instead you provided {self.ionMass.unit.physical_type}"
        self.particle = ionName

        self.electronTemperature = electronTemperature
        self.ionTemperature = ionTemperature

        self.magneticField = magneticField
        self.electricField = electricField


        self.compute_properties()

    def compute_properties(self):
        """the properties of the plasma should not evolves, so they can all be computed"""

        Te_k = self.electronTemperature.to(u.K, equivalencies=u.temperature_energy())

        self.Debye_length =  Debye_length(Te_k, self.plasmaDensity)

        self.ionPlasmaFrequency = plasma_frequency(self.plasmaDensity, self.particle)

        self.electronPlasmaFrequency = plasma_frequency(self.plasmaDensity)

        self.BohmSpeed = self.Debye_length * self.ionPlasmaFrequency/u.rad

        Ti_k = self.ionTemperature.to(u.K, equivalencies=u.temperature_energy())

        self.vti = thermal_speed(Ti_k, self.particle)

        self.vte =  thermal_speed(Te_k)

        self.driftSpeed = (self.electricField/self.magneticField).to(u.m/u.s)