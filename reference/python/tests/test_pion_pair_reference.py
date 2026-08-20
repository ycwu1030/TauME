from pathlib import Path
import unittest

from tauamp_reference.tautau.pion_pair import BeamState, ElectroweakParameters, FourMomentum, TauPairKinematicPoint, pion_pair_components


ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / 'tests/fixtures/tautau/pion_pair_reference_cases.dat'


def parse_cases(path):
    cases = []
    current = None
    for raw_line in path.read_text().splitlines():
        parts = raw_line.split()
        if not parts:
            continue
        key, values = parts[0], parts[1:]
        if key == 'case':
            current = {'name': values[0]}
        elif key == 'end':
            cases.append(current)
            current = None
        elif key in {'electron', 'positron', 'tau_minus', 'tau_plus', 'pion_minus', 'pion_plus'}:
            current[key] = FourMomentum(*map(float, values))
        elif key == 'parameters':
            current[key] = ElectroweakParameters(*map(float, values))
        elif key == 'polarizations':
            current[key] = tuple(map(float, values))
        elif key == 'tau_mass':
            current[key] = float(values[0])
        elif key == 'expected':
            current[key] = tuple(map(float, values))
        else:
            raise ValueError(f'unknown fixture field: {key}')
    return cases


class PionPairReferenceTest(unittest.TestCase):
    def test_shared_conditional_fixtures(self):
        for case in parse_cases(FIXTURE):
            point = TauPairKinematicPoint(
                BeamState(case['electron'], case['positron'], *case['polarizations']), case['tau_minus'], case['tau_plus'],
                case['tau_mass'])
            actual = pion_pair_components(point, case['pion_minus'], case['pion_plus'], case['parameters'])
            for field, expected in zip(('sm', 'f2_real', 'f2_imaginary', 'f3_real', 'f3_imaginary'), case['expected']):
                self.assertAlmostEqual(getattr(actual, field), expected, delta=3.e-11, msg=f"{case['name']} {field}")


if __name__ == '__main__':
    unittest.main()
