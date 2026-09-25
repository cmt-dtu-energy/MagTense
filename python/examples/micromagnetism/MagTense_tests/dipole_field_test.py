"""
Test the far field of a uniformly magnetised cube against the analytical point dipole.

A single cubic prism of side L, magnetised along z, is evaluated at random points on spheres
of increasing radius. Far from the cube the field has to approach that of a point dipole of
moment m = M*L^3, and the way it approaches it is itself a check: a cube has no quadrupole
moment by symmetry, so the leading correction is the octupole and the relative deviation from
the dipole has to fall as R^-4.

That makes this two tests in one. The far field checks the absolute accuracy of the
magnetostatic kernel, and the exponent checks that the near field has the right multipole
structure, which a plain threshold at one distance would not catch.

Running the file executes the test and saves a figure. ``run_test()`` returns the same result
as a list of checks, which is the contract the combined suite in testMagTenseFunctions.py
expects.
"""

# General modules
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Magtense stuff
from magtense.magstatics import Tiles, run_simulation

# Style settings for plots
plt.rcParams['font.size'] = 15
plt.rcParams['text.usetex'] = False

#%% Settings

mu0 = 4 * np.pi * 1e-7

L = 1.2e-6              # Side length of the cube [m]
M_rem = 1.2 / mu0       # Remanent magnetisation, along z [A/m]

# Distances at which to sample, in units of the cube side length. The smallest has to stay
# outside the cube by a wide margin for the multipole expansion to mean anything, and the
# largest is limited by the point at which the deviation reaches the floating point noise
# floor of the kernel itself.
distances = np.array([5, 10, 20, 30, 50])

# Number of random directions per distance. They are drawn from a normal distribution and
# normalised, which is uniform on the sphere, so the test covers all orientations relative to
# the magnetisation rather than only the symmetry axes.
n_points = 40
seed = 42

# The far field is compared at this distance, in units of L
far_field_distance = 30

# How far the field may deviate from the point dipole at far_field_distance [%]
far_field_tol = 1e-3

# The relative deviation has to fall as R^-4. How far the fitted exponent may sit from -4 [-]
exponent_tol = 0.15

output_dir = Path(__file__).resolve().parent


#%% Analytical reference


def dipole_field(r_vec: np.ndarray, m_vec: np.ndarray) -> np.ndarray:
    """H field of a point dipole of moment m_vec at the position r_vec, in A/m."""
    r_mag = np.linalg.norm(r_vec)
    r_hat = r_vec / r_mag
    return (3 * r_hat * np.dot(m_vec, r_hat) - m_vec) / (4 * np.pi * r_mag**3)


def deviation_at(distance: float, directions: np.ndarray) -> np.ndarray:
    """Relative deviation from the point dipole at one distance, one value per point [%].

    The whole field vector is compared rather than only its magnitude, so a field of the right
    strength pointing the wrong way still fails.
    """
    pts = directions * (L * distance)

    tiles = Tiles(n=1, M_rem=M_rem, tile_type=[2], easy_axis=[0, 0, 1])
    tiles.size = ([L, L, L], 0)
    tiles.offset = ([0, 0, 0], 0)
    _, H_pts = run_simulation(tiles, np.asarray(pts))

    # Moment of the cube: the magnetisation times its volume
    m_vec = np.array([0.0, 0.0, M_rem]) * L**3

    deviations = np.empty(len(pts))
    for i, point in enumerate(pts):
        H_dipole = dipole_field(point, m_vec)
        deviations[i] = (np.linalg.norm(H_pts[i, :] - H_dipole)
                         / np.linalg.norm(H_dipole) * 100)
    return deviations


#%% Test


def run_test(plotting: bool = True) -> list[dict]:
    rng = np.random.default_rng(seed)
    directions = rng.normal(size=(n_points, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)

    print(f'Comparing the field of a cube with a point dipole at {n_points} random '
          f'directions per distance')
    print(f"{'R / L':>8} {'max deviation [%]':>20} {'mean deviation [%]':>20}")
    print('-' * 50)

    max_deviations = np.empty(len(distances))
    mean_deviations = np.empty(len(distances))
    for i, distance in enumerate(distances):
        deviations = deviation_at(distance, directions)
        max_deviations[i] = np.max(deviations)
        mean_deviations[i] = np.mean(deviations)
        print(f'{distance:>8} {max_deviations[i]:>20.3e} {mean_deviations[i]:>20.3e}')
    print('-' * 50)

    # The far field itself
    far_field_index = int(np.argmin(np.abs(distances - far_field_distance)))
    far_field_deviation = float(max_deviations[far_field_index])
    print(f'Deviation at R = {far_field_distance}L: {far_field_deviation:.3e} % '
          f'(limit {far_field_tol:.3e} %)')

    # The way it is approached. A cube has no quadrupole moment, so the leading correction is
    # the octupole and the relative deviation goes as R^-4.
    exponent = float(np.polyfit(np.log(distances), np.log(max_deviations), 1)[0])
    exponent_error = abs(exponent + 4)
    print(f'Fitted falloff exponent: {exponent:.3f} (expected -4, '
          f'off by {exponent_error:.3f})')

    if plotting:
        fig, ax = plt.subplots(layout='constrained', figsize=(6, 4))
        ax.loglog(distances, max_deviations, 'o', color='crimson', markersize=8,
                  label='Calculation, largest of the points')
        ax.loglog(distances, mean_deviations, 's', color='steelblue', markersize=6,
                  label='Calculation, mean of the points')
        reference = max_deviations[0] * (distances / distances[0])**-4.0
        ax.loglog(distances, reference, '--', color='forestgreen',
                  label=r'Octupole correction, $R^{-4}$')
        ax.set_xlabel(r'$R \, / \, L$')
        ax.set_ylabel(r'$\text{Deviation from point dipole } [\%]$')
        ax.legend(loc='best', fontsize='small')
        ax.grid(True, which='both', alpha=0.3)

        # Save the figure beside this script and close it so the script does not
        # open an interactive plotting window.
        figure_path = output_dir / 'dipole_field_test.png'
        fig.savefig(figure_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved figure to {figure_path}')

    return [
        {
            'check': f'far field matches a point dipole at R = {far_field_distance}L',
            'value': far_field_deviation,
            'limit': far_field_tol,
            'passed': far_field_deviation < far_field_tol,
        },
        {
            'check': 'deviation falls as R^-4, as the octupole correction',
            'value': exponent_error,
            'limit': exponent_tol,
            'passed': exponent_error < exponent_tol,
        },
    ]


if __name__ == '__main__':
    checks = run_test()
    print('dipole_field_test '
          + ('PASSED' if all(c['passed'] for c in checks) else 'FAILED'))
