from dataclasses import dataclass
from enum import Enum
import math

import numpy as np

from .kinematics import BeamState, FourMomentum, TauPairKinematicPoint


_TOLERANCE = 1.0e-9


@dataclass(frozen=True)
class HadronicTauPairObservation:
    beams: BeamState
    visible_minus: FourMomentum
    visible_plus: FourMomentum
    tau_mass: float


class HadronicTauPairSolutionStatus(str, Enum):
    NO_SOLUTION = 'no_solution'
    ONE_SOLUTION = 'one_solution'
    TWO_SOLUTIONS = 'two_solutions'
    NON_UNIQUE = 'non_unique'


@dataclass(frozen=True)
class HadronicTauPairKinematicSolution:
    point: TauPairKinematicPoint
    neutrino_minus: FourMomentum
    neutrino_plus: FourMomentum


@dataclass(frozen=True)
class HadronicTauPairSolutionSet:
    status: HadronicTauPairSolutionStatus
    solutions: tuple[HadronicTauPairKinematicSolution, ...]


class HadronicTauPairKinematicSolver:
    @staticmethod
    def solve(observation):
        tau_mass = observation.tau_mass
        sqrt_s = observation.beams.sqrt_s()
        if tau_mass <= 0.0 or sqrt_s <= 2.0 * tau_mass:
            return _no_solution()

        beam_total = observation.beams.electron + observation.beams.positron
        if beam_total.energy <= 0.0:
            return _no_solution()
        to_pair_cm = -beam_total.spatial() / beam_total.energy
        to_lab = -to_pair_cm
        visible_minus_cm = observation.visible_minus.boosted(to_pair_cm)
        visible_plus_cm = observation.visible_plus.boosted(to_pair_cm)
        minus_spatial = visible_minus_cm.spatial()
        plus_spatial = visible_plus_cm.spatial()
        minus_momentum = float(np.linalg.norm(minus_spatial))
        plus_momentum = float(np.linalg.norm(plus_spatial))
        if minus_momentum < _TOLERANCE or plus_momentum < _TOLERANCE:
            return _no_solution()

        tau_energy = sqrt_s / 2.0
        tau_momentum_squared = tau_energy ** 2 - tau_mass ** 2
        if tau_momentum_squared <= _TOLERANCE ** 2:
            return _no_solution()
        tau_momentum = math.sqrt(tau_momentum_squared)
        minus_mass_squared = visible_minus_cm.mass_squared()
        plus_mass_squared = visible_plus_cm.mass_squared()
        if minus_mass_squared < -_TOLERANCE or plus_mass_squared < -_TOLERANCE:
            return _no_solution()

        minus_direction = minus_spatial / minus_momentum
        plus_direction = plus_spatial / plus_momentum
        c_minus = (2.0 * tau_energy * visible_minus_cm.energy - tau_mass ** 2 - minus_mass_squared) / \
                  (2.0 * tau_momentum * minus_momentum)
        c_plus = (tau_mass ** 2 + plus_mass_squared - 2.0 * tau_energy * visible_plus_cm.energy) / \
                 (2.0 * tau_momentum * plus_momentum)
        if c_minus < -1.0 - _TOLERANCE or c_minus > 1.0 + _TOLERANCE or \
           c_plus < -1.0 - _TOLERANCE or c_plus > 1.0 + _TOLERANCE:
            return _no_solution()
        c_minus = _bound_unit(c_minus)
        c_plus = _bound_unit(c_plus)
        dot = float(minus_direction @ plus_direction)
        denominator = max(0.0, 1.0 - dot ** 2)

        if denominator <= _TOLERANCE ** 2:
            if abs(c_plus - dot * c_minus) > _TOLERANCE:
                return _no_solution()
            if abs(abs(c_minus) - 1.0) <= _TOLERANCE:
                solution = _make_solution(observation, visible_minus_cm, visible_plus_cm,
                                          c_minus * minus_direction, tau_energy, tau_momentum, to_lab)
                return _one_or_no_solution(solution)
            return HadronicTauPairSolutionSet(HadronicTauPairSolutionStatus.NON_UNIQUE, ())

        base = ((c_minus - dot * c_plus) * minus_direction +
                (c_plus - dot * c_minus) * plus_direction) / denominator
        discriminant = 1.0 - float(base @ base)
        if discriminant < -_TOLERANCE:
            return _no_solution()
        if discriminant <= _TOLERANCE:
            base_norm = float(np.linalg.norm(base))
            if base_norm < _TOLERANCE:
                return _no_solution()
            solution = _make_solution(observation, visible_minus_cm, visible_plus_cm,
                                      base / base_norm, tau_energy, tau_momentum, to_lab)
            return _one_or_no_solution(solution)

        normal = np.cross(minus_direction, plus_direction) / math.sqrt(denominator)
        scale = math.sqrt(discriminant)
        first = _make_solution(observation, visible_minus_cm, visible_plus_cm,
                               base + scale * normal, tau_energy, tau_momentum, to_lab)
        second = _make_solution(observation, visible_minus_cm, visible_plus_cm,
                                base - scale * normal, tau_energy, tau_momentum, to_lab)
        if first is None or second is None:
            return _no_solution()
        return HadronicTauPairSolutionSet(HadronicTauPairSolutionStatus.TWO_SOLUTIONS, (first, second))


def _bound_unit(value):
    return min(1.0, max(-1.0, value))


def _no_solution():
    return HadronicTauPairSolutionSet(HadronicTauPairSolutionStatus.NO_SOLUTION, ())


def _one_or_no_solution(solution):
    if solution is None:
        return _no_solution()
    return HadronicTauPairSolutionSet(HadronicTauPairSolutionStatus.ONE_SOLUTION, (solution,))


def _make_solution(observation, visible_minus_cm, visible_plus_cm, direction, tau_energy, tau_momentum, to_lab):
    tau_minus_cm = FourMomentum(*(tau_momentum * component for component in direction), tau_energy)
    tau_plus_cm = FourMomentum(-tau_minus_cm.px, -tau_minus_cm.py, -tau_minus_cm.pz, tau_energy)
    neutrino_minus_cm = tau_minus_cm - visible_minus_cm
    neutrino_plus_cm = tau_plus_cm - visible_plus_cm
    if neutrino_minus_cm.energy < -_TOLERANCE or neutrino_plus_cm.energy < -_TOLERANCE or \
       abs(neutrino_minus_cm.mass_squared()) > _TOLERANCE or abs(neutrino_plus_cm.mass_squared()) > _TOLERANCE:
        return None
    tau_minus_lab = tau_minus_cm.boosted(to_lab)
    tau_plus_lab = tau_plus_cm.boosted(to_lab)
    return HadronicTauPairKinematicSolution(
        TauPairKinematicPoint(observation.beams, tau_minus_lab, tau_plus_lab, observation.tau_mass),
        neutrino_minus_cm.boosted(to_lab), neutrino_plus_cm.boosted(to_lab))
