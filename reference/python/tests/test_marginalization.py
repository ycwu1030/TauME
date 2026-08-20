from pathlib import Path
import unittest

from tauamp_reference.tautau.marginalization import LinearComponents, average_conditional_observables, ratio_of_average_components


ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / 'tests/fixtures/tautau/pion_pair_ensemble_cases.dat'


def parse_cases(path):
    cases = []
    current = None
    for raw_line in path.read_text().splitlines():
        parts = raw_line.split()
        if not parts:
            continue
        key, values = parts[0], parts[1:]
        if key == 'case':
            current = {'name': values[0], 'components': []}
        elif key == 'component':
            numbers = list(map(float, values))
            current['components'].append((LinearComponents(*numbers[:5]), numbers[5]))
        elif key == 'expected_average_of_ratios':
            current[key] = tuple(map(float, values))
        elif key == 'expected_ratio_after_average':
            current[key] = tuple(map(float, values))
        elif key == 'end':
            cases.append(current)
            current = None
        else:
            raise ValueError(f'unknown fixture field: {key}')
    return cases


class MarginalizationTest(unittest.TestCase):
    def test_shared_ensemble_fixtures(self):
        for case in parse_cases(FIXTURE):
            for actual, expected in ((average_conditional_observables(case['components']), case['expected_average_of_ratios']),
                                     (ratio_of_average_components(case['components']), case['expected_ratio_after_average'])):
                for field, value in zip(('f2_real', 'f2_imaginary', 'f3_real', 'f3_imaginary'), expected):
                    self.assertAlmostEqual(getattr(actual, field), value, delta=3.e-11, msg=f"{case['name']} {field}")


if __name__ == '__main__':
    unittest.main()
