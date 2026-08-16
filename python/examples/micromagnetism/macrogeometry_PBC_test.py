"""
Test if periodic boundary conditions as modeled by the macrogeometry method are implemented correctly in the Magtense micromagnetics solver.
SHORT EXPLANATION : If the plotted polar angle is constant at the critical point, the code works.

LONG EXPLANATION :
The system is a pair of uniformly magnetised cubes of sidelength a, placed side by side along x.
There are periodic boundary conditions (PBCs) so that the cubes are repeated at regular intervals of A along the entire x-axis.
Both cubes have a magnetocrystalline easy-axis along z, but the magnetostatic interaction between the tiles and
with their infinite copies produce an effective easy-axis along x as well.
The initial configuration for the two magnetic moments is m1 = (cos π/4, 0, sin π/4) and m2 = (cos π/4, 0, -sin π/4),
i.e. at 45 degrees to the x-axis. By symmetry θ_1 = -θ_2 = θ throughout, so the motion boils down to a single polar angle, θ.

Using the results of the article : "Exact demagnetisation field for periodic one-dimensional array of rectangular prisms" by Durhuus et al.,
one can derive the parameter combinations where the effective anisotropy is zero. In this special case, the magnetic moments will remain at
their starting values, while for all other parameter combinations, they will either align along the x-axis or anti-align along z.
The theoretically predicted critical points is in excellent agreement with the numerical simulation.
The test is to verify the agreement still holds, i.e. θ is essentially constant at the critical point.

Note : The magnetisation and anisotropy constants are chosen so that the anisotropy and pair interaction
mostly cancel, ensuring the θ(t) behaviour is decided by the periodically repeating copies.

Other note : The whole setup is rotated so that the axis-of-periodicity becomes x, y or z in turn, which
validates the macrogeometry along all three axes independently.

The test runs three simulations at the critical spacing, one per axis of periodicity, and requires θ to
stay put. Two further simulations at spacings on either side of the critical point act as controls: they
must move θ towards π/2 and 0 respectively, in the direction predicted by the analytical prefactor below.
Without those controls a solver that simply froze the magnetisation would pass the constant-θ check.
"""

# General modules
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Specialised modules
from scipy.special import polygamma

try:
    from tqdm import tqdm
except ModuleNotFoundError:
    def tqdm(x):
        return x

# Magtense stuff
from magtense.micromag import MicromagProblem

# Style settings for plots
plt.rcParams['font.size'] = 16
plt.rcParams['text.usetex'] = False
plt.rcParams['text.latex.preamble'] = r'\usepackage{physics}'

#%% Settings

# Distance between neighbouring domain copies in macrogeometry
A_crit = 1.18369410e-07 # Critical point where the torque is zero at all theta (constant theta)
A_wide = 2e-7       # Should go towards theta = π/2 (increasing theta)
A_narrow = 0.7e-7   # Should go towards theta = 0 (decreasing theta)

# Which axes to test the PBC along. 0 for x-axis, 1 for y-axis and 2 for z-axis
testAxes = (0, 1, 2)

# The axis used for the two control simulations that check the test is sensitive at all
controlAxis = 0

# Timestepping
t_end = 1e-7                                 # Total time simulated [s]
t_step = 1e-10                               # Requested output spacing [s]
nTimesteps = int(round(t_end / t_step)) + 1  # Include both time endpoints
# The output spacing only sets how densely theta(t) is sampled - RKSuite picks its own
# internal steps. Requesting 1e5 output times made the ODE driver overflow the stack.

# How constant theta has to be at the critical point to count as passed [rad]
angle_tol = 2e-2
# The controls only approach their fixed point asymptotically, so they get a looser bound [rad]
control_tol = 0.1

# Plot settings
results_dir = Path(__file__).resolve().parent / "results"

# Micromagnetic solver settings
cuda = False
cvode = False

#%% Fixed settings (don't change these)

# Constants of nature
mu0 = 4*np.pi * 1e-7    # Vacuum permeability [N/A^2]
gamma = 0    # Gyromagnetic ratio set to zero to turn off precession so motion stays 2D

# Material properties
Ms = 1.0               # Saturation magnetisation in magnetic region [T]
Ms = Ms/mu0            # [A/m]
alpha = 1              # Gilbert damping [-]
eta = alpha/(1 + alpha**2) * 2.21*1e5  # Damping constant [m/(A*s)]
K = 2.7e4              # Uniaxial anisotropy constant [J/m^3]
T = 0                  # Temperature [K]
Aex = 1e-20                # Exchange constant [J/m]

# Applied field
def fct_h_ext(t) -> np.ndarray:
    e_z = np.array([[0, 0, 1]]).T   # Unit vector along z as 3-by-1 column vector
    amplitude = 0    # sine-wave for amplitude
    return np.transpose(amplitude * e_z)    # Has to end with the shape (nt, 3)

# Side lengths of individual micromagnetic cells
a = 1e-8    # [m]

# Number of copies of the simulated domain on each side along the axis of periodicity
n_copies = 1000

#%% Analytical results

# Analytical solution
def PBC_1D_analytical(pts, a, b, c, pos, DeltaX, nExact):
    """Approximate analytical solution for the demagnetisation tensor contribution
    resulting from infinite, distant copies of a prism-cell in a system with periodic
    boundary conditions along the x-axis.

    pts : array of positions to evaluate the field at
    a, b, c : prism shape
    pos : position of the original prism
    DeltaX : distance between neighbouring prism copies along the x-axis
    nExact : the contributions of the original prism and the 2*nExact
             nearest neighbours are excluded, so they can be computed
             exactly elsewhere
    """
    ptsX, posX = pts[:, 0], pos[:, 0]
    x = ptsX[:, np.newaxis] - posX[np.newaxis, :]

    prefactor = b*c / (4*np.pi*DeltaX**2)
    firstTerm = polygamma(1, (a/2 - x)/DeltaX + nExact+1)
    secondTerm = -polygamma(1, (-a/2 - x)/DeltaX + nExact+1)
    thirdTerm = polygamma(1, (a/2 + x)/DeltaX + nExact+1)
    fourthTerm = -polygamma(1, (-a/2 + x)/DeltaX + nExact+1)
    Nxx = prefactor*(firstTerm + secondTerm + thirdTerm + fourthTerm)
    N = np.zeros([*Nxx.shape, 3, 3])
    Nyy = Nzz = -1/2 * Nxx
    N[..., 0, 0] = Nxx
    N[..., 1, 1] = Nyy
    N[..., 2, 2] = Nzz
    return N

def prefactor(q, K, M):
    """The rate of change of potential energy wrt. polar angle is dU/dθ = prefactor * sin(θ)*cos(θ),
    so the sign of the prefactor determines whether the angle goes to π/2 (negative) or 0 (positive) and
    the magnitude decides how quickly.
    q = a/(2*A)"""

    # Pair interaction
    fPair = 2*np.arctan(1/np.sqrt(3)) - np.arctan(3/np.sqrt(11)) - np.arctan(1/(3*np.sqrt(11)))

    # Interaction between a tile and its infinite copies
    fSelf = 3*polygamma(1, 1-q) - 3*polygamma(1, 1+q)

    # Interaction between one tile and the infinite copies of the other tile
    fInter = 1/2*(polygamma(1, 1+q) + polygamma(1, 1-3*q) - polygamma(1, 1+3*q) - polygamma(1, 1-q))

    # Anisotropy contribution
    fAni = -2*np.pi*K/(mu0 * M**2)

    return 2*mu0 * M**2/np.pi * (fPair + q**2 * fSelf + q**2 * fInter + fAni)

# # A_crit found by:
# f = lambda A : np.abs(prefactor(a/(2*A), K, Ms))
# from scipy.optimize import minimize_scalar
# minimize_scalar(f, bracket=[1.18e-7, 1.185e-7, 1.19e-7], tol=1e-20, method='bounded', bounds=[1.18e-7, 1.19e-7])

#%% Magtense computations

def run_case(testAxis: int, A: float):
    """Simulate the two-cube system with periodic copies spaced A apart along testAxis.

    Returns the time array and the polar angles of both moments. The angle is measured
    from the axis of periodicity towards the easy-axis, which is z when the periodicity
    is along x or y and x when the periodicity is along z. Both choices keep the easy-axis
    perpendicular to the axis of periodicity, so the three cases are rotations of one another.
    """
    perpAxis = 2 if testAxis in {0, 1} else 0

    # Initial magnetisations, at +-45 degrees to the axis of periodicity
    theta_p = np.array([np.pi/4, -np.pi/4])
    m0_pv = np.zeros([2, 3])
    m0_pv[:, testAxis] = np.cos(theta_p)
    m0_pv[:, perpAxis] = np.sin(theta_p)

    # Define micromagnetic problem
    grid_type = "uniform"
    grid_L = [a, a, a]                  # Grid size along x, y and z for simulated domain
    res = [1, 1, 1]                     # Grid resolution, i.e. number of cells along x, y and z
    # Make two tiles along the axis of periodicity
    grid_L[testAxis] = 2*a
    res[testAxis] = 2

    # Define macrogeometry.
    # The macrogeometry is a (2nx + 1) X (2ny + 1) X (2nz + 1) grid of domain copies
    n_macro, shiftVec = np.zeros(3), np.zeros(3)
    n_macro[testAxis] = n_copies  # Number of copies on each side along x, y and z
    shiftVec[testAxis] = A  # How far to shift the domain copies along x, y and z

    problem = MicromagProblem(
        res=res,
        A0=Aex,
        Ms=Ms,
        K0=K,
        alpha=eta,
        gamma=gamma,
        m0=m0_pv,
        cuda=cuda,
        cvode=cvode,
        n_macro=n_macro,
        shiftVec=shiftVec,
        grid_L=grid_L,
        grid_type=grid_type,
        usereturnhall=True,
        solver='dynamic',
        # solver='explicit',
        T=T,
    )

    # The FMM demag path segfaults on this two tile, effectively one dimensional geometry when
    # the library is built with USE_FMM3D=1, and an octree is pointless for two tiles anyway.
    # std_problem_3.py disables it the same way.
    problem.use_fmm = 0

    # Easy-axis anisotropy perpendicular to the axis of periodicity
    problem.u_ea = np.zeros([2, 3])
    problem.u_ea[:, perpAxis] = 1

    result = problem.run_simulation(
        t_end=t_end,
        nt=nTimesteps,
        fct_h_ext=fct_h_ext,
        nt_h_ext=2,     # Two samples completely define the constant zero field.
    )

    t_n, M_out = result[:2]
    M_npv = M_out[:, :, 0, :]
    M1_nv, M2_nv = M_npv[:, 0, :], M_npv[:, 1, :]
    # By symmetry theta_2 = -theta_1, so negating the perpendicular component of the second
    # moment should reproduce the first angle.
    theta1_n = np.arctan2(M1_nv[:, perpAxis], M1_nv[:, testAxis])
    theta2_n = np.arctan2(-M2_nv[:, perpAxis], M2_nv[:, testAxis])
    return t_n, theta1_n, theta2_n

#%% Run the simulations

axis_names = ("x", "y", "z")


def run_test(plotting: bool = True) -> list[dict]:
    """Run the macrogeometry PBC test and return the checks it consists of.

    Each check is a dict with the keys 'check', 'value', 'limit' and 'passed', where the test
    passes when value < limit. This is the contract used by testMagTenseFunctions.py.
    """
    # The analytical prefactor decides where theta is heading, and vanishes at the critical spacing
    print('Analytical prefactor (>0 drives theta to 0, <0 drives theta to pi/2, 0 means no torque):')
    for A in (A_crit, A_wide, A_narrow):
        print(f'  A = {A:.4e} m : {prefactor(a/(2*A), K, Ms):+.4e}')

    runs = [(axis, A_crit, 'critical') for axis in testAxes]
    runs += [(controlAxis, A_wide, 'wide'), (controlAxis, A_narrow, 'narrow')]

    print('Run simulations')
    results = {}
    for testAxis, A, label in tqdm(runs):
        results[(testAxis, label)] = run_case(testAxis, A)
    print('Done running simulations')

    checks = []

    # At the critical spacing theta must stay at its initial value of pi/4 along every axis
    for testAxis in testAxes:
        t_n, theta1_n, theta2_n = results[(testAxis, 'critical')]
        drift = float(np.max(np.abs(theta1_n - np.pi/4)))
        asymmetry = float(np.max(np.abs(theta1_n - theta2_n)))
        ok = drift < angle_tol and asymmetry < angle_tol
        verdict = 'works' if ok else 'FAILED'
        print(f'Macrogeometry PBC {verdict} along {axis_names[testAxis]}: '
              f'max |theta - pi/4| = {drift:.2e} rad, '
              f'max |theta_1 + theta_2| = {asymmetry:.2e} rad')
        checks.append({
            'check': f'PBC along {axis_names[testAxis]}: theta constant at critical A',
            'value': drift, 'limit': angle_tol, 'passed': drift < angle_tol,
        })
        checks.append({
            'check': f'PBC along {axis_names[testAxis]}: theta_1 = -theta_2',
            'value': asymmetry, 'limit': angle_tol, 'passed': asymmetry < angle_tol,
        })

    # The controls must actually move, and in the direction the analytical prefactor predicts
    for label, A, target in (('wide', A_wide, np.pi/2), ('narrow', A_narrow, 0.0)):
        t_n, theta1_n, _ = results[(controlAxis, label)]
        moved = float(theta1_n[-1] - theta1_n[0])
        expected = target - np.pi/4
        distance = float(abs(theta1_n[-1] - target))
        right_way = np.sign(moved) == np.sign(expected)
        ok = distance < control_tol and right_way
        verdict = 'works' if ok else 'FAILED'
        print(f'Control at A = {A:.4e} m {verdict}: theta went from {theta1_n[0]:.4f} to '
              f'{theta1_n[-1]:.4f} rad, expected {target:.4f}')
        checks.append({
            'check': f'control, A {label}: theta reaches {target:.3f} rad',
            # A wrong direction is reported as a failure however close the endpoint happens to be
            'value': distance if right_way else float('inf'),
            'limit': control_tol, 'passed': ok,
        })

    if plotting:
        plot_results(results)

    return checks


def plot_results(results) -> None:
    """Two figures: the critical case against its controls, and the drift along each axis."""
    results_dir.mkdir(parents=True, exist_ok=True)

    # ---- The critical spacing against the two controls, periodicity along x ----------------
    fig, ax = plt.subplots(layout='constrained', figsize=(6, 4))
    # Plot control tests with non-critical spacing
    label_wide = r'$A > A_\text{crit}, \text{should go to } \frac{\pi}{2}$'
    label_narrow = r'$A < A_\text{crit}, \text{should go to } 0$'
    labelDict = {'narrow': label_narrow, 'wide': label_wide}
    for label, colour in zip(('narrow', 'wide'), ('darkorange', 'purple'), strict=True):
        t_n, theta1_n, _ = results[(controlAxis, label)]
        ax.plot(t_n * 1e9, theta1_n, color=colour, linestyle='--', label=labelDict[label])
    # Plot test with critical spacing and PBC along x
    t_n, theta1_n, _ = results[(controlAxis, 'critical')]
    ax.plot(t_n * 1e9, theta1_n, color='forestgreen',
            label=r'$A = A_\text{crit}, \text{should be constant}$')
    ax.set_xlabel(r'$\text{Time } [\mathrm{ns}]$')
    ax.set_ylabel(r'$\text{Polar angle}$')
    ax.set_yticks([0, np.pi/4, np.pi/2])
    ax.set_yticklabels(['0', r'$\pi/4$', r'$\pi/2$'])
    ax.legend(loc='best', fontsize='14')

    # Save the figure beside the other micromagnetic example results and close it so the
    # script does not open an interactive plotting window.
    figure_path = results_dir / "macrogeometry_PBC_test_periodic_along_x.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to {figure_path}")

    # ---- The drift at the critical spacing along each axis, against the accepted band ------
    fig, ax = plt.subplots(layout='constrained', figsize=(6, 4))
    critical_colours = ('forestgreen', 'crimson', 'navy')
    for i, testAxis in enumerate(testAxes):
        t_n, theta1_n, _ = results[(testAxis, 'critical')]
        ax.plot(t_n * 1e9, theta1_n - theta1_n[0],
                color=critical_colours[i % len(critical_colours)],
                label=f'$\\text{{PBCs along {axis_names[testAxis]}}}$')
    ax.set_xlabel(r'$\text{Time } [\mathrm{ns}]$')
    ax.set_ylabel(r'$\text{Change in polar angle}$')
    angle_window = 10**(np.floor(np.log10(angle_tol*10)))
    ax.set_yticks([-angle_window, 0, angle_window])
    ax.set_ylim(-angle_window, angle_window)
    ax.set_yticklabels([f'$-10^{{{np.log10(angle_window):.0f}}}$', r'$0$',
                        f'$10^{{{np.log10(angle_window):.0f}}}$'])
    ax.fill_between(np.array([t_n[0]*1e9, t_n[-1]*1e9]),
                    np.array([-angle_tol, -angle_tol]), np.array([angle_tol, angle_tol]),
                    alpha=0.2, color='gray', label='Accepted interval')
    ax.legend(loc='best', fontsize='14')

    figure_path = results_dir / "macrogeometry_PBC_test_periodic_along_x_y_or_z.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to {figure_path}")


if __name__ == '__main__':
    all_checks = run_test()
    print('macrogeometry_PBC_test '
          + ('PASSED' if all(c['passed'] for c in all_checks) else 'FAILED'))
