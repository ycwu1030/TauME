from dataclasses import dataclass
from math import isfinite


@dataclass(frozen=True)
class LinearComponents:
    sm: float
    f2_real: float
    f2_imaginary: float
    f3_real: float
    f3_imaginary: float


@dataclass(frozen=True)
class LinearObservableComponents:
    f2_real: float
    f2_imaginary: float
    f3_real: float
    f3_imaginary: float


def conditional_observables(components):
    if components.sm == 0.0:
        raise ValueError('conditional observable requires a nonzero SM component')
    return LinearObservableComponents(components.f2_real / components.sm, components.f2_imaginary / components.sm,
                                      components.f3_real / components.sm, components.f3_imaginary / components.sm)


def _checked(entries):
    if not entries:
        raise ValueError('linear-component ensemble cannot be empty')
    total = 0.0
    for _, weight in entries:
        if not isfinite(weight) or weight < 0.0:
            raise ValueError('linear-component ensemble weights must be finite and non-negative')
        total += weight
    if not isfinite(total) or total <= 0.0:
        raise ValueError('linear-component ensemble total weight must be finite and positive')
    return total


def average_conditional_observables(entries):
    total = _checked(entries)
    value = [0.0, 0.0, 0.0, 0.0]
    for components, weight in entries:
        conditional = conditional_observables(components)
        for index, component in enumerate((conditional.f2_real, conditional.f2_imaginary, conditional.f3_real,
                                            conditional.f3_imaginary)):
            value[index] += weight * component
    return LinearObservableComponents(*(component / total for component in value))


def ratio_of_average_components(entries):
    _checked(entries)
    total = [0.0] * 5
    for components, weight in entries:
        for index, component in enumerate((components.sm, components.f2_real, components.f2_imaginary,
                                            components.f3_real, components.f3_imaginary)):
            total[index] += weight * component
    return conditional_observables(LinearComponents(*total))
