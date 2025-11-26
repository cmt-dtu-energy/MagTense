from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike, NDArray
from typing import Optional


def rand_spherical_cap(
    cone_angle_degree: float,
    cone_dir: Optional[ArrayLike] = None,
    N: int = 1,
    rng: Optional[np.random.Generator] = None
) -> NDArray[np.float64]:
    """
    Generate random unit vectors within a spherical cap (cone) on the unit sphere.

    Parameters
    ----------
    cone_angle_degree : float
        Half-angle of the cone in degrees.
    cone_dir : array_like of shape (3,), optional
        Direction of the cone axis. Defaults to [0, 0, 1].
    N : int, optional
        Number of random vectors to generate. Defaults to 1.
    rng : numpy.random.Generator, optional
        Random number generator. If None, uses np.random.default_rng().

    Returns
    -------
    r : ndarray of shape (3, N)
        Columns are random unit vectors inside the cone.
    """
    if cone_dir is None:
        cone_dir_arr = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    else:
        cone_dir_arr = np.asarray(cone_dir, dtype=np.float64).reshape(3)

    if rng is None:
        rng = np.random.default_rng()

    # Convert to radians
    cone_angle = np.deg2rad(cone_angle_degree)

    # Generate points on the spherical cap around the north pole
    z: NDArray[np.float64] = rng.random(N) * (1.0 - np.cos(cone_angle)) + np.cos(cone_angle)
    phi: NDArray[np.float64] = rng.random(N) * 2.0 * np.pi
    xy_radius: NDArray[np.float64] = np.sqrt(1.0 - z**2)

    x = xy_radius * np.cos(phi)
    y = xy_radius * np.sin(phi)

    r = np.vstack((x, y, z))  # shape (3, N)

    # If aligned with +z, return directly
    if np.allclose(cone_dir_arr, np.array([0.0, 0.0, 1.0], dtype=np.float64)):
        return r

    # Normalised vectors
    north = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    n_cone = cone_dir_arr / np.linalg.norm(cone_dir_arr)

    # Rotation axis
    u = np.cross(north, n_cone)
    u_norm = np.linalg.norm(u)

    # Handle degenerate case (antiparallel)
    if u_norm == 0.0:
        # 180° rotation around x-axis
        R = np.array([[1.0,  0.0,  0.0],
                      [0.0, -1.0,  0.0],
                      [0.0,  0.0, -1.0]], dtype=np.float64)
        return R @ r

    u = u / u_norm
    rot = np.arccos(np.clip(np.dot(n_cone, north), -1.0, 1.0))

    ux, uy, uz = u

    cross_matrix = np.array([[    0.0, -uz,   uy],
                             [   uz,   0.0, -ux],
                             [  -uy,   ux,  0.0]], dtype=np.float64)

    R = (
        np.cos(rot) * np.eye(3, dtype=np.float64)
        + np.sin(rot) * cross_matrix
        + (1.0 - np.cos(rot)) * np.outer(u, u)
    )

    return R @ r
