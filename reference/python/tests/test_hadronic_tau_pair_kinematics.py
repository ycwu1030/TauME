from math import cos, sin, sqrt
from pathlib import Path
import unittest

from tauamp_reference.tautau.hadronic_kinematics import (
    HadronicTauPairKinematicSolver,
    HadronicTauPairObservation,
)
from tauamp_reference.tautau.kinematics import BeamState, FourMomentum


ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / 'tests/fixtures/tautau/hadronic_tau_pair_twofold_cases.dat'
TOLERANCE = 5.0e-10


def parse_cases(path):
    cases = []
    current = None
    for raw_line in path.read_text().splitlines():
        parts = raw_line.split()
        if not parts:
            continue
        key, values = parts[0], parts[1:]
        if key == 'case':
            current = {'name': values[0], 'solutions': []}
        elif key in {'electron', 'positron', 'visible_minus', 'visible_plus'}:
            current[key] = FourMomentum(*map(float, values))
        elif key == 'polarizations':
            current[key] = tuple(map(float, values))
        elif key == 'tau_mass':
            current[key] = float(values[0])
        elif key == 'status':
            current[key] = values[0]
        elif key == 'solution':
            values = list(map(float, values))
            current['solutions'].append(tuple(FourMomentum(*values[index:index + 4]) for index in range(0, 16, 4)))
        elif key == 'end':
            cases.append(current)
            current = None
        else:
            raise ValueError(f'unknown fixture field: {key}')
    return cases


def observation(case):
    return HadronicTauPairObservation(
        BeamState(case['electron'], case['positron'], *case['polarizations']),
        case['visible_minus'], case['visible_plus'], case['tau_mass'])


def assert_momentum_close(test_case, actual, expected):
    for field in ('px', 'py', 'pz', 'energy'):
        test_case.assertAlmostEqual(getattr(actual, field), getattr(expected, field), delta=TOLERANCE, msg=field)


def visible_from_tau_rest(direction, tau, tau_mass, visible_mass):
    rest_energy = (tau_mass ** 2 + visible_mass ** 2) / (2.0 * tau_mass)
    rest_momentum = (tau_mass ** 2 - visible_mass ** 2) / (2.0 * tau_mass)
    return FourMomentum(*(rest_momentum * component for component in direction), rest_energy).boosted(tau.spatial() / tau.energy)


class HadronicTauPairKinematicsTest(unittest.TestCase):
    def assert_solution_constraints(self, solution, event):
        total = event.beams.electron + event.beams.positron
        assert_momentum_close(self, solution.point.tau_minus + solution.point.tau_plus, total)
        self.assertAlmostEqual(solution.point.tau_minus.mass_squared(), event.tau_mass ** 2, delta=TOLERANCE)
        self.assertAlmostEqual(solution.point.tau_plus.mass_squared(), event.tau_mass ** 2, delta=TOLERANCE)
        self.assertAlmostEqual(solution.neutrino_minus.mass_squared(), 0.0, delta=TOLERANCE)
        self.assertAlmostEqual(solution.neutrino_plus.mass_squared(), 0.0, delta=TOLERANCE)
        self.assertGreaterEqual(solution.neutrino_minus.energy, -TOLERANCE)
        self.assertGreaterEqual(solution.neutrino_plus.energy, -TOLERANCE)
        assert_momentum_close(self, solution.point.tau_minus - event.visible_minus, solution.neutrino_minus)
        assert_momentum_close(self, solution.point.tau_plus - event.visible_plus, solution.neutrino_plus)

    def test_shared_twofold_fixtures(self):
        self.assertTrue(FIXTURE.is_file())
        for case in parse_cases(FIXTURE):
            event = observation(case)
            result = HadronicTauPairKinematicSolver.solve(event)
            self.assertEqual(result.status.value, case['status'], case['name'])
            self.assertEqual(len(result.solutions), len(case['solutions']), case['name'])
            for actual, expected in zip(result.solutions, case['solutions']):
                for actual_momentum, expected_momentum in zip(
                    (actual.point.tau_minus, actual.point.tau_plus, actual.neutrino_minus, actual.neutrino_plus), expected):
                    assert_momentum_close(self, actual_momentum, expected_momentum)
                self.assert_solution_constraints(actual, event)

    def test_tangent_cones_produce_one_solution(self):
        tau_mass = 1.77686
        energy = 2.1
        tau_momentum = sqrt(energy ** 2 - tau_mass ** 2)
        electron = FourMomentum(0.0, 0.0, energy, energy)
        positron = FourMomentum(0.0, 0.0, -energy, energy)
        tau_minus = FourMomentum(0.0, 0.0, tau_momentum, energy)
        tau_plus = FourMomentum(0.0, 0.0, -tau_momentum, energy)
        event = HadronicTauPairObservation(
            BeamState(electron, positron, 0.0, 0.0),
            visible_from_tau_rest((sin(0.7), 0.0, cos(0.7)), tau_minus, tau_mass, 0.6),
            visible_from_tau_rest((sin(1.0), 0.0, cos(1.0)), tau_plus, tau_mass, 0.8), tau_mass)

        result = HadronicTauPairKinematicSolver.solve(event)

        self.assertEqual(result.status.value, 'one_solution')
        self.assertEqual(len(result.solutions), 1)
        assert_momentum_close(self, result.solutions[0].point.tau_minus, tau_minus)
        self.assert_solution_constraints(result.solutions[0], event)

    def test_compatible_collinear_cones_are_non_unique(self):
        tau_mass = 1.77686
        energy = 2.1
        tau_momentum = sqrt(energy ** 2 - tau_mass ** 2)
        electron = FourMomentum(0.0, 0.0, energy, energy)
        positron = FourMomentum(0.0, 0.0, -energy, energy)
        cosine = 0.5
        event = HadronicTauPairObservation(
            BeamState(electron, positron, 0.0, 0.0),
            FourMomentum(0.0, 0.0, tau_mass ** 2 / (2.0 * (energy - tau_momentum * cosine)),
                         tau_mass ** 2 / (2.0 * (energy - tau_momentum * cosine))),
            FourMomentum(0.0, 0.0, tau_mass ** 2 / (2.0 * (energy + tau_momentum * cosine)),
                         tau_mass ** 2 / (2.0 * (energy + tau_momentum * cosine))), tau_mass)

        result = HadronicTauPairKinematicSolver.solve(event)

        self.assertEqual(result.status.value, 'non_unique')
        self.assertEqual(result.solutions, ())


if __name__ == '__main__':
    unittest.main()
