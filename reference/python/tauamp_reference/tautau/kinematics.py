from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class FourMomentum:
    px: float
    py: float
    pz: float
    energy: float

    def vector(self):
        return np.array([self.energy, self.px, self.py, self.pz], dtype=float)

    def spatial(self):
        return np.array([self.px, self.py, self.pz], dtype=float)

    def mass_squared(self):
        return self.energy ** 2 - self.px ** 2 - self.py ** 2 - self.pz ** 2

    def boosted(self, beta):
        beta = np.asarray(beta, dtype=float)
        beta_squared = float(beta @ beta)
        if beta_squared >= 1.0:
            raise ValueError('boost velocity must be subluminal')
        if beta_squared < 1.0e-18:
            return self
        gamma = 1.0 / math.sqrt(1.0 - beta_squared)
        momentum = self.spatial()
        beta_dot_p = float(beta @ momentum)
        factor = (gamma - 1.0) * beta_dot_p / beta_squared + gamma * self.energy
        shifted = momentum + factor * beta
        return FourMomentum(*shifted, gamma * (self.energy + beta_dot_p))

    def __add__(self, other):
        return FourMomentum(self.px + other.px, self.py + other.py, self.pz + other.pz, self.energy + other.energy)

    def __sub__(self, other):
        return FourMomentum(self.px - other.px, self.py - other.py, self.pz - other.pz, self.energy - other.energy)


@dataclass(frozen=True)
class BeamState:
    electron: FourMomentum
    positron: FourMomentum
    electron_polarization: float
    positron_polarization: float

    def sqrt_s(self):
        return math.sqrt((self.electron + self.positron).mass_squared())


@dataclass(frozen=True)
class TauPairKinematicPoint:
    beams: BeamState
    tau_minus: FourMomentum
    tau_plus: FourMomentum
    tau_mass: float
