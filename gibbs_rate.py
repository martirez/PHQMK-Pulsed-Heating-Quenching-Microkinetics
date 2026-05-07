
import cantera as ct
import numpy as np

class SurfaceCubicEaData(ct.ExtensibleRateData):
    __slots__ = ("T",)

    def __init__(self):
        self.T = None

    def update(self, thermo):
        T = thermo.T
        if self.T != T:
            self.T = T
            return True
        return False

@ct.extension(name="surface-cubic-Ea", data=SurfaceCubicEaData)
class SurfaceCubicEaRate(ct.ExtensibleRate):
    def set_parameters(self, params, units):
        self.A = params.convert_rate_coeff("A", units)
        self.b = params["b"]
        # Cubic coefficients for Ea(T)
        self.Ea_coeffs = [
            float(params["Ea0"]),      # a0
            float(params["Ea1"]),      # a1
            float(params["Ea2"]),      # a2
            float(params["Ea3"])       # a3
        ]

    def eval(self, data):
        T_adj = data.T - 1000
        Ea_T = (
            self.Ea_coeffs[0]
            + self.Ea_coeffs[1] * T_adj
            + self.Ea_coeffs[2] * T_adj**2
            + self.Ea_coeffs[3] * T_adj**3
        )
        if Ea_T < 0.00:
            Ea_T = 0.00
        r = self.A * data.T**self.b * np.exp(-Ea_T / (8.31432 * data.T))
        return r
