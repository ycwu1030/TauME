from math import cos, sin, sqrt, tanh
import unittest

from tauamp_reference.tautau.hadronic_kinematics import HadronicTauPairSolutionStatus
from tauamp_reference.tautau.kinematics import BeamState, FourMomentum
from tauamp_reference.tautau.marginalization import average_conditional_observables, ratio_of_average_components
from tauamp_reference.tautau.pion_pair import ElectroweakParameters, pion_pair_components
from tauamp_reference.tautau.visible_pion_pair import (
    KinematicBranchCombination,
    PionPairVisibleEventEvaluation,
    PionPairVisibleEventEvaluator,
)


TOLERANCE = 3.0e-11


def rotate_y(momentum, angle):
    cosine, sine = cos(angle), sin(angle)
    return FourMomentum(cosine * momentum.px + sine * momentum.pz, momentum.py,
                        -sine * momentum.px + cosine * momentum.pz, momentum.energy)


def visible_from_tau_rest(direction, tau, tau_mass, visible_mass):
    rest_energy = (tau_mass ** 2 + visible_mass ** 2) / (2.0 * tau_mass)
    rest_momentum = (tau_mass ** 2 - visible_mass ** 2) / (2.0 * tau_mass)
    return FourMomentum(*(rest_momentum * component for component in direction), rest_energy).boosted(tau.spatial() / tau.energy)


def assert_components_close(test_case, actual, expected):
    for field in ('sm', 'f2_real', 'f2_imaginary', 'f3_real', 'f3_imaginary'):
        test_case.assertAlmostEqual(getattr(actual, field), getattr(expected, field), delta=TOLERANCE, msg=field)


def assert_observables_close(test_case, actual, expected):
    for field in ('f2_real', 'f2_imaginary', 'f3_real', 'f3_imaginary'):
        test_case.assertAlmostEqual(getattr(actual, field), getattr(expected, field), delta=TOLERANCE, msg=field)


class VisiblePionPairTest(unittest.TestCase):
    def setUp(self):
        self.tau_mass = 1.77686
        self.energy = 2.1
        self.pion_mass = 0.13957
        momentum = sqrt(self.energy ** 2 - self.tau_mass ** 2)
        electron = FourMomentum(0.0, 0.0, self.energy, self.energy)
        positron = FourMomentum(0.0, 0.0, -self.energy, self.energy)
        self.beams = BeamState(electron, positron, 0.37, -0.41)
        tau_minus = FourMomentum(momentum * sin(0.91), 0.0, momentum * cos(0.91), self.energy)
        tau_plus = FourMomentum(-tau_minus.px, -tau_minus.py, -tau_minus.pz, self.energy)
        self.pion_minus = visible_from_tau_rest((0.2, -0.3, sqrt(0.87)), tau_minus, self.tau_mass, self.pion_mass)
        self.pion_plus = visible_from_tau_rest((-0.4, 0.1, sqrt(0.83)), tau_plus, self.tau_mass, self.pion_mass)
        self.parameters = ElectroweakParameters(0.313, 0.48, 91.1876, 2.4952)
        self.evaluator = PionPairVisibleEventEvaluator(self.parameters)

    def test_twofold_candidates_carry_equal_weights_and_direct_components(self):
        result = self.evaluator.evaluate(self.beams, self.pion_minus, self.pion_plus, self.tau_mass)

        self.assertEqual(result.status, HadronicTauPairSolutionStatus.TWO_SOLUTIONS)
        self.assertEqual(len(result.entries), 2)
        for entry in result.entries:
            self.assertEqual(entry.origin, 'analytic_branch')
            self.assertAlmostEqual(entry.weight, 0.5, delta=TOLERANCE)
            assert_components_close(self, entry.components,
                                    pion_pair_components(entry.point, self.pion_minus, self.pion_plus, self.parameters))

        weighted = tuple((entry.components, entry.weight) for entry in result.entries)
        expected_average = average_conditional_observables(weighted)
        expected_ratio = ratio_of_average_components(weighted)
        self.assertNotAlmostEqual(expected_average.f3_real, expected_ratio.f3_real, delta=TOLERANCE)
        assert_observables_close(self, result.observable_components(KinematicBranchCombination.AVERAGE_CONDITIONAL_OBSERVABLES),
                                 expected_average)
        assert_observables_close(self, result.observable_components(KinematicBranchCombination.RATIO_OF_AVERAGE_COMPONENTS),
                                 expected_ratio)

    def test_branch_order_does_not_change_either_combination(self):
        result = self.evaluator.evaluate(self.beams, self.pion_minus, self.pion_plus, self.tau_mass)
        reversed_result = PionPairVisibleEventEvaluation(result.status, tuple(reversed(result.entries)))

        for rule in KinematicBranchCombination:
            assert_observables_close(self, result.observable_components(rule), reversed_result.observable_components(rule))

    def test_asymmetric_lab_input_preserves_components(self):
        reference = self.evaluator.evaluate(self.beams, self.pion_minus, self.pion_plus, self.tau_mass)
        beta = (0.0, 0.0, tanh(0.35))
        boosted = self.evaluator.evaluate(
            BeamState(self.beams.electron.boosted(beta), self.beams.positron.boosted(beta), 0.37, -0.41),
            self.pion_minus.boosted(beta), self.pion_plus.boosted(beta), self.tau_mass)

        self.assertEqual(boosted.status, HadronicTauPairSolutionStatus.TWO_SOLUTIONS)
        for actual, expected in zip(boosted.entries, reference.entries):
            assert_components_close(self, actual.components, expected.components)

    def test_one_solution_uses_weight_one(self):
        momentum = sqrt(self.energy ** 2 - self.tau_mass ** 2)
        tau_minus = FourMomentum(0.0, 0.0, momentum, self.energy)
        tau_plus = FourMomentum(0.0, 0.0, -momentum, self.energy)
        production_angle = 0.91
        result = self.evaluator.evaluate(
            BeamState(self.beams.electron, self.beams.positron, 0.0, 0.0),
            rotate_y(visible_from_tau_rest((sin(0.7), 0.0, cos(0.7)), tau_minus, self.tau_mass, self.pion_mass),
                     production_angle),
            rotate_y(visible_from_tau_rest((sin(1.0), 0.0, cos(1.0)), tau_plus, self.tau_mass, self.pion_mass),
                     production_angle), self.tau_mass)

        self.assertEqual(result.status, HadronicTauPairSolutionStatus.ONE_SOLUTION)
        self.assertEqual(len(result.entries), 1)
        self.assertAlmostEqual(result.entries[0].weight, 1.0, delta=TOLERANCE)
        assert_observables_close(self, result.observable_components(KinematicBranchCombination.AVERAGE_CONDITIONAL_OBSERVABLES),
                                 result.observable_components(KinematicBranchCombination.RATIO_OF_AVERAGE_COMPONENTS))

    def test_no_solution_and_non_unique_do_not_produce_observables(self):
        no_solution = self.evaluator.evaluate(self.beams, FourMomentum(0.0, 0.0, 0.0, 1.0), self.pion_plus, self.tau_mass)
        self.assertEqual(no_solution.status, HadronicTauPairSolutionStatus.NO_SOLUTION)
        self.assertEqual(no_solution.entries, ())
        with self.assertRaises(ValueError):
            no_solution.observable_components(KinematicBranchCombination.RATIO_OF_AVERAGE_COMPONENTS)

        momentum = sqrt(self.energy ** 2 - self.tau_mass ** 2)
        cosine = 0.5
        minus_energy = self.tau_mass ** 2 / (2.0 * (self.energy - momentum * cosine))
        plus_energy = self.tau_mass ** 2 / (2.0 * (self.energy + momentum * cosine))
        non_unique = self.evaluator.evaluate(
            BeamState(self.beams.electron, self.beams.positron, 0.0, 0.0),
            FourMomentum(0.0, 0.0, minus_energy, minus_energy), FourMomentum(0.0, 0.0, plus_energy, plus_energy),
            self.tau_mass)
        self.assertEqual(non_unique.status, HadronicTauPairSolutionStatus.NON_UNIQUE)
        self.assertEqual(non_unique.entries, ())
        with self.assertRaises(ValueError):
            non_unique.observable_components(KinematicBranchCombination.AVERAGE_CONDITIONAL_OBSERVABLES)


if __name__ == '__main__':
    unittest.main()
