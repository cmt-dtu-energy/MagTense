from math import asin, atan2
from warnings import warn

import numpy as np

from magtense import magstatics


class Utility:
    """Utility class for geometry-related method."""

    @staticmethod
    def get_euler_angles(R: np.ndarray) -> np.ndarray:
        """
        Returns the Euler angles (in radians) corresponding to the given 3x3 matrix of
        orthogonal eigenvectors. Assumes that the input matrix R is a proper rotation
        matrix (i.e., det(R) = 1). This fucntionality is also implemented in
        scipy.spatial.transform.Rotation.from_matrix(R).as_euler('xyz', degrees=False).
        """
        sy = np.sqrt(R[0, 0] ** 2 + R[1, 0] ** 2)
        if sy < 1e-6:
            # Singular case: cos(y) is close to zero and we can't determine x and z
            x = atan2(-R[1, 2], R[1, 1])
            y = asin(-R[2, 0])
            z = 0
        else:
            x = atan2(R[2, 1], R[2, 2])
            y = atan2(-R[2, 0], sy)
            z = atan2(R[1, 0], R[0, 0])
        return np.array([x, y, z])


class StandardGeometry:
    """
    Standard geometries for magnetostatics, including MNIST problem configuations.
    """

    pass


class DemagMatrix:
    """
    Wrapper for individual demagnetization matrices. These are real, symmetric
    3x3 matrices with trace zero or one, corresponding to a scaling along three
    principal axes, so that when applied to a magnetisation vector they return the
    demagnetization field at an associated point.
    """

    def __init__(self, values: np.ndarray):
        if not values.size == 9:
            raise ValueError("DemagMatrix must be initialized with a 3x3 array.")
        if not values.shape == (3, 3):
            values = values.reshape((3, 3))
        self.values = values
        if not self.is_valid:
            warn(f"DemagMatrix {self} is not valid.")

    def __repr__(self):
        return f"DemagMatrix({self.values})"

    @classmethod
    def from_mag(cls, tile: magstatics.Tiles, point: np.ndarray):
        """Calculate the demagnetization matrix for a given tile and points."""
        return cls(magstatics.get_demag_tensor(tile, point).squeeze())

    @property
    def _is_symmetric(self):
        """Check if the matrix is symmetric."""
        return np.allclose(self.values, self.values.T)

    @property
    def _is_traceless(self):
        """Check if the matrix is traceless."""
        return np.isclose(self.values.trace(), 0.0)

    @property
    def _is_unit_trace(self):
        """Check if the matrix has unit trace."""
        return np.isclose(self.values.trace(), -1.0)

    @property
    def is_valid(self):
        """Check if the matrix is valid."""
        return self._is_symmetric and (self._is_traceless or self._is_unit_trace)

    @property
    def invariant(self):
        """
        Returns the five invariant parameters of the demagnetization tensor N: the two
        largest eigenvalues and the three Euler angles of the orthogonal matrix of
        eigenvectors.
        """
        w, v = np.linalg.eigh(self.values)
        return np.concatenate([w[1:], Utility.get_euler_angles(v)])


class DemagTensor:
    pass
