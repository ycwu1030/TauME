from dataclasses import dataclass
from enum import Enum
from math import isfinite

from .hadronic_kinematics import (
    HadronicTauPairKinematicSolver,
    HadronicTauPairObservation,
    HadronicTauPairSolutionStatus,
)
from .marginalization import average_conditional_observables, ratio_of_average_components
from .pion_pair import pion_pair_components


class KinematicBranchCombination(str, Enum):
    AVERAGE_CONDITIONAL_OBSERVABLES = 'average_conditional_observables'
    RATIO_OF_AVERAGE_COMPONENTS = 'ratio_of_average_components'


@dataclass(frozen=True)
class PionPairHypothesisComponents:
    point: object
    components: object
    weight: float
    origin: str


@dataclass(frozen=True)
class PionPairVisibleEventEvaluation:
    status: HadronicTauPairSolutionStatus
    entries: tuple[PionPairHypothesisComponents, ...]

    def __post_init__(self):
        expected_entries = {
            HadronicTauPairSolutionStatus.ONE_SOLUTION: 1,
            HadronicTauPairSolutionStatus.TWO_SOLUTIONS: 2,
        }.get(self.status, 0)
        if len(self.entries) != expected_entries:
            raise ValueError('visible pion-pair evaluation has incompatible solution entries')
        total_weight = 0.0
        for entry in self.entries:
            if not isfinite(entry.weight) or entry.weight < 0.0:
                raise ValueError('visible pion-pair evaluation weights must be finite and non-negative')
            total_weight += entry.weight
        if self.entries and (not isfinite(total_weight) or total_weight <= 0.0):
            raise ValueError('visible pion-pair evaluation total weight must be finite and positive')

    def observable_components(self, combination):
        weighted = tuple((entry.components, entry.weight) for entry in self.entries)
        if combination == KinematicBranchCombination.AVERAGE_CONDITIONAL_OBSERVABLES:
            return average_conditional_observables(weighted)
        if combination == KinematicBranchCombination.RATIO_OF_AVERAGE_COMPONENTS:
            return ratio_of_average_components(weighted)
        raise ValueError(f'unknown kinematic-branch combination: {combination}')


class PionPairVisibleEventEvaluator:
    def __init__(self, parameters):
        self._parameters = parameters

    def evaluate(self, beams, pion_minus, pion_plus, tau_mass):
        solution_set = HadronicTauPairKinematicSolver.solve(
            HadronicTauPairObservation(beams, pion_minus, pion_plus, tau_mass))
        weight = 1.0 / len(solution_set.solutions) if solution_set.solutions else 0.0
        entries = tuple(
            PionPairHypothesisComponents(
                solution.point,
                pion_pair_components(solution.point, pion_minus, pion_plus, self._parameters),
                weight,
                'analytic_branch')
            for solution in solution_set.solutions)
        return PionPairVisibleEventEvaluation(solution_set.status, entries)
