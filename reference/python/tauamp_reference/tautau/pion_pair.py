from dataclasses import dataclass
import math
import numpy as np

from .marginalization import LinearComponents

_METRIC = np.array([1.0, -1.0, -1.0, -1.0])
_I4 = np.identity(4, dtype=complex)


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


@dataclass(frozen=True)
class BeamState:
    electron: FourMomentum
    positron: FourMomentum
    electron_polarization: float
    positron_polarization: float

    def sqrt_s(self):
        total = self.electron + self.positron
        return math.sqrt(total.energy ** 2 - total.px ** 2 - total.py ** 2 - total.pz ** 2)


@dataclass(frozen=True)
class TauPairKinematicPoint:
    beams: BeamState
    tau_minus: FourMomentum
    tau_plus: FourMomentum
    tau_mass: float


@dataclass(frozen=True)
class ElectroweakParameters:
    electron_charge: float
    sin_theta_w: float
    z_mass: float
    z_width: float


def _normalized(value):
    norm = float(np.linalg.norm(value))
    if norm < 1.e-9:
        raise ValueError('undefined normalized direction')
    return value / norm


def _rotate_to_electron_axis(value, electron):
    z_axis = _normalized(electron.spatial())
    reference = np.array([0.0, 1.0, 0.0])
    if abs(float(reference @ z_axis)) > .9:
        reference = np.array([1.0, 0.0, 0.0])
    x_axis = _normalized(np.cross(reference, z_axis))
    y_axis = np.cross(z_axis, x_axis)
    spatial = value.spatial()
    return FourMomentum(float(spatial @ x_axis), float(spatial @ y_axis), float(spatial @ z_axis), value.energy)


def _pair_cm(point):
    total = point.beams.electron + point.beams.positron
    beta = -total.spatial() / total.energy
    electron_before_rotation = point.beams.electron.boosted(beta)
    return tuple(_rotate_to_electron_axis(value.boosted(beta), electron_before_rotation)
                 for value in (point.beams.electron, point.beams.positron, point.tau_minus, point.tau_plus)), beta, electron_before_rotation


def _spin_axes(electron, tau_minus):
    longitudinal = _normalized(tau_minus.spatial())
    normal = _normalized(-np.cross(electron.spatial(), tau_minus.spatial()))
    return np.cross(normal, longitudinal), normal, longitudinal


def _gamma_matrices():
    gamma0 = np.diag([1.0, 1.0, -1.0, -1.0]).astype(complex)
    pauli = (np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex),
             np.array([[1, 0], [0, -1]], complex))
    result = [gamma0]
    for sigma in pauli:
        matrix = np.zeros((4, 4), complex)
        matrix[:2, 2:] = sigma
        matrix[2:, :2] = -sigma
        result.append(matrix)
    return tuple(result)


_GAMMA = _gamma_matrices()
_GAMMA5 = 1j * _GAMMA[0] @ _GAMMA[1] @ _GAMMA[2] @ _GAMMA[3]


def _slash(vector):
    return sum((_METRIC[index] * vector[index] * _GAMMA[index] for index in range(4)), np.zeros((4, 4), complex))


def _bar(matrix):
    return _GAMMA[0] @ matrix.conjugate().T @ _GAMMA[0]


def _spin_vector(momentum, axis, mass):
    spatial_projection = float(momentum.spatial() @ axis)
    factor = spatial_projection / (mass * (momentum.energy + mass))
    return np.array([spatial_projection / mass, *(axis + momentum.spatial() * factor)])


def _bosons(parameters):
    cosine = math.sqrt(1.0 - parameters.sin_theta_w ** 2)
    scale = parameters.electron_charge / (2.0 * parameters.sin_theta_w * cosine)
    vector = scale * (-.5 + 2.0 * parameters.sin_theta_w ** 2)
    axial = -.5 * scale
    return ((0.0, 0.0, -parameters.electron_charge, 0.0, -parameters.electron_charge, 0.0),
            (parameters.z_mass, parameters.z_width, vector, axial, vector, axial))


def _standard_vertices(vector, axial):
    return tuple(gamma @ (vector * _I4 - axial * _GAMMA5) for gamma in _GAMMA)


def _sigma_q(mu, q):
    return sum((_METRIC[nu] * q[nu] * .5j * (_GAMMA[mu] @ _GAMMA[nu] - _GAMMA[nu] @ _GAMMA[mu]) for nu in range(4)),
               np.zeros((4, 4), complex))


def _tau_vertices(boson, photon, mass, q, sector, form_factor):
    if not photon or sector == 'sm':
        return _standard_vertices(boson[4], boson[5])
    values = []
    for mu in range(4):
        dipole = _sigma_q(mu, q)
        if sector == 'f2':
            values.append(boson[4] * dipole * (1j * form_factor / (2.0 * mass)))
        else:
            values.append(boson[4] * dipole * (form_factor / (2.0 * mass)) @ _GAMMA5)
    return tuple(values)


def _pair_contract(point, cm, left, right, left_photon, right_photon, left_sector, right_sector, left_factor, right_factor):
    electron, positron, tau_minus, tau_plus = cm
    q = tau_minus.vector() + tau_plus.vector()
    electron_density = (_I4 + point.beams.electron_polarization * _GAMMA5) @ _slash(electron.vector()) * .5
    positron_density = (_I4 - point.beams.positron_polarization * _GAMMA5) @ _slash(positron.vector()) * .5
    left_electron = _standard_vertices(left[2], left[3])
    right_electron = tuple(_bar(value) for value in _standard_vertices(right[2], right[3]))
    left_tau = _tau_vertices(left, left_photon, point.tau_mass, q, left_sector, left_factor)
    right_tau = tuple(_bar(value) for value in _tau_vertices(right, right_photon, point.tau_mass, q, right_sector, right_factor))
    lepton = np.empty((4, 4), complex)
    for mu in range(4):
        for nu in range(4):
            lepton[mu, nu] = np.trace(positron_density @ left_electron[mu] @ electron_density @ right_electron[nu])
    s = point.beams.sqrt_s() ** 2
    denominator_left = complex(s - left[0] ** 2, left[0] * left[1])
    denominator_right = complex(s - right[0] ** 2, right[0] * right[1])
    return lepton, tuple(_METRIC[mu] * left_tau[mu] for mu in range(4)), tuple(_METRIC[mu] * right_tau[mu] for mu in range(4)), 1.0 / (denominator_left * denominator_right.conjugate())


def _coefficients(point, cm, contract):
    _, _, tau_minus, tau_plus = cm
    lepton, left_tau, right_tau, weight = contract
    base_minus = _slash(tau_minus.vector()) + point.tau_mass * _I4
    base_plus = _slash(tau_plus.vector()) - point.tau_mass * _I4
    axes = _spin_axes(cm[0], tau_minus)
    values = np.empty((4, 4), complex)
    for minus_index in range(4):
        for plus_index in range(4):
            minus_projector = base_minus if minus_index == 0 else base_minus @ (_I4 + _GAMMA5 @ _slash(_spin_vector(tau_minus, axes[minus_index - 1], point.tau_mass))) * .5
            plus_projector = base_plus if plus_index == 0 else base_plus @ (_I4 + _GAMMA5 @ _slash(_spin_vector(tau_plus, axes[plus_index - 1], point.tau_mass))) * .5
            values[minus_index, plus_index] = weight * sum(lepton[mu, nu] * np.trace(minus_projector @ left_tau[mu] @ plus_projector @ right_tau[nu]) for mu in range(4) for nu in range(4))
    result = np.empty((4, 4), complex)
    result[0, 0] = values[0, 0]
    for i in range(1, 4):
        result[i, 0] = 2.0 * values[i, 0] - result[0, 0]
        result[0, i] = 2.0 * values[0, i] - result[0, 0]
    for i in range(1, 4):
        for j in range(1, 4):
            result[i, j] = 4.0 * values[i, j] - result[0, 0] - result[i, 0] - result[0, j]
    return result


def _sector_matrix(point, cm, model, sector, factor):
    total = np.zeros((4, 4), complex)
    def add(left, right, left_sector, right_sector, left_factor, right_factor):
        nonlocal total
        total += _coefficients(point, cm, _pair_contract(point, cm, model[left], model[right], left == 0, right == 0,
                                                          left_sector, right_sector, left_factor, right_factor))
    if sector == 'sm':
        for left in range(2):
            for right in range(2):
                add(left, right, 'sm', 'sm', 0j, 0j)
    else:
        add(0, 0, sector, 'sm', factor, 0j)
        add(0, 0, 'sm', sector, 0j, factor)
        add(0, 1, sector, 'sm', factor, 0j)
        add(1, 0, 'sm', sector, 0j, factor)
    if np.max(np.abs(total.imag)) > 2.e-9:
        raise RuntimeError('summed production coefficient has an unexpected imaginary residue')
    return total.real


def _pion_analyser(point, cm, beta_to_cm, electron_before_rotation, pion_lab, tau_minus):
    pion_cm = _rotate_to_electron_axis(pion_lab.boosted(beta_to_cm), electron_before_rotation)
    tau = cm[2] if tau_minus else cm[3]
    pion_rest = pion_cm.boosted(-tau.spatial() / tau.energy)
    direction = _normalized(pion_rest.spatial())
    axes = _spin_axes(cm[0], cm[2])
    value = np.array([direction @ axis for axis in axes])
    return value if tau_minus else -value


def _contract(matrix, minus, plus):
    result = matrix[0, 0] + matrix[1:, 0] @ minus + matrix[0, 1:] @ plus
    return float(result + minus @ matrix[1:, 1:] @ plus)


def pion_pair_components(point, pion_minus_lab, pion_plus_lab, parameters):
    cm, beta_to_cm, electron_before_rotation = _pair_cm(point)
    model = _bosons(parameters)
    minus = _pion_analyser(point, cm, beta_to_cm, electron_before_rotation, pion_minus_lab, True)
    plus = _pion_analyser(point, cm, beta_to_cm, electron_before_rotation, pion_plus_lab, False)
    matrices = (_sector_matrix(point, cm, model, 'sm', 0j), _sector_matrix(point, cm, model, 'f2', 1.0 + 0j),
                _sector_matrix(point, cm, model, 'f2', 1j), _sector_matrix(point, cm, model, 'f3', 1.0 + 0j),
                _sector_matrix(point, cm, model, 'f3', 1j))
    return LinearComponents(*(_contract(matrix, minus, plus) for matrix in matrices))
