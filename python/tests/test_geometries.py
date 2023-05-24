from math import pi, sqrt

import numpy as np
import pytest
from magtense.geometries import DemagMatrix, Utility
from magtense.magstatics import Tiles


@pytest.fixture
def point():
    return np.array([[0.0, 0.0, 0.0]])


@pytest.fixture
def tile():
    return Tiles(1)


@pytest.fixture
def valid(tile, point):
    return DemagMatrix.from_mag(tile=tile, point=point)


@pytest.fixture
def invalid():
    with pytest.warns(UserWarning, match="DemagMatrix is not valid."):
        return DemagMatrix(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]))


@pytest.fixture
def rotations():
    R1 = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    R2 = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    R3 = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    return R1, R2, R3


def test_demag_matrix_construction(valid):
    assert valid.values.shape == (3, 3), "DemagMatrix should be 3x3."


def test_demag_matrix_symmetry(valid, invalid):
    assert valid._is_symmetric, "DemagMatrix should be symmetric."
    assert not invalid._is_symmetric, "DemagMatrix should not be symmetric."


def test_demag_matrix_trace(valid, invalid):
    assert valid._is_traceless, "DemagMatrix should be traceless."
    assert not invalid._is_traceless, "DemagMatrix should not be traceless."
    assert not valid._is_unit_trace, "DemagMatrix should not have unit trace."


def test_demag_matrix_validity(valid, invalid):
    assert valid.is_valid, "DemagMatrix should be valid."
    assert not invalid.is_valid, "DemagMatrix should not be valid."


def test_get_euler_angles(rotations):
    R1, R2, R3 = rotations
    assert np.allclose(Utility.get_euler_angles(R1), np.array([pi / 2, 0, 0]))
    assert np.allclose(Utility.get_euler_angles(R2), np.array([0, pi / 2, 0]))
    assert np.allclose(Utility.get_euler_angles(R3), np.array([0, 0, pi / 2]))


def test_invariant():
    N1 = DemagMatrix(np.array([[1, 0, 0], [0, 0, -1], [0, -1, 0]]))
    N2 = DemagMatrix(np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]]))
    N3 = DemagMatrix(np.array([[1, 0, 0], [0, -1, -1], [0, -1, 0]]))
    assert np.allclose(N1.invariant, np.array([1, 1, pi / 2, pi / 4, -pi / 2]))
    assert np.allclose(N2.invariant, np.array([1, 1, -pi / 2, -pi / 4, pi]))
    assert np.allclose(
        N3.invariant, np.array([(sqrt(5) - 1) / 2, 1, -pi / 2, 0.55357436, -pi / 2])
    )
