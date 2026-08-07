import cantera as ct
import numpy as np

class SurfaceCubicEaData(ct.ExtensibleRateData):
    __slots__ = ("T","site_density")

    def __init__(self):
        self.T = None
        self.site_density = None

    def update(self, thermo):
        T = thermo.T
        site_density = thermo.site_density
        if self.T != T or self.site_density != site_density:
            self.T = T
            self.site_density = site_density
            return True
        return False

@ct.extension(name="surface-cubic-Ea", data=SurfaceCubicEaData)
class SurfaceCubicEaRate(ct.ExtensibleRate):
    def set_parameters(self, params, units):
        self.surface_order = float(params["surface-order"])
        self.gas_order = float(params["gas-order"])
        self.standard_pressure = float(params.get("standard-pressure", 100000))
        # Cubic coefficients for Ea(T)
        self.Ea_coeffs = [
            float(params["Ea0"]),      # a0
            float(params["Ea1"]),      # a1
            float(params["Ea2"]),      # a2
            float(params["Ea3"])       # a3
        ]

    def eval(self, data):
        Gamma = data.site_density
        T = data.T
        k_tst = ct.boltzmann * T / ct.planck
        
        T_adj = T - 1000
        Ea_T = (
            self.Ea_coeffs[0]
            + self.Ea_coeffs[1] * T_adj
            + self.Ea_coeffs[2] * T_adj**2
            + self.Ea_coeffs[3] * T_adj**3
        )
        if Ea_T < 0.00:
            Ea_T = 0.00
        r = k_tst * (Gamma ** (1.0 - self.surface_order)) * ((8314.32 * T / self.standard_pressure)**(self.gas_order)) * np.exp(-Ea_T / (8.31432 * T))
        
        return r
