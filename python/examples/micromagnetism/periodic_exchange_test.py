"""Test the implementation of periodic exchange coupling in a uniform grid

A 6 X 6 X 6 cube with periodic boundary conditions is divided into 8 sections that are internally exchange coupled but
independent of one another. The test is passed if the magnetisation vectors are identical within each connected section but distinct
between sections due to random initial conditions.

The 8 sections are a 4 X 4 X 4 cube of center tiles, the 8 corners, the xFace (4 tiles on both the left and right side),
the yFace (4 tiles on both the front and back face), the zFace (4 tiles on both the top and bottom faces), the
xyEdge (8 tiles on the edges perpendicular to the xy plane), the xyEdge and the yzEdge.
The remaining cells are numerically inert spacers and there is no dipole interaction, so the separate sections are effectively
decoupled, giving 8 tests of the exchange coupling. The center tiles test basic exchange, while the rest check the periodic boundary conditions.

Running the file executes the test and saves a figure. ``run_test()`` returns the same result as a
list of checks, which is the contract the combined suite in testMagTenseFunctions.py expects.
"""

# General modules
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Magtense stuff
from magtense.micromag import MicromagProblem

# Style settings for plots
plt.rcParams['font.size'] = 16
plt.rcParams['text.usetex'] = False
plt.rcParams['text.latex.preamble'] = r'\usepackage{physics}'

#%% Settings

# Timestepping
t_end = 1e-8                      # Total time simulated [s]
t_step = 1e-10                    # Requested output spacing [s]
nTimesteps = int(round(t_end / t_step)) + 1  # Include both time endpoints

# Plot settings
Ndata = 10  # Only the endpoints are needed
dataperiod = int(nTimesteps / Ndata)
results_dir = Path(__file__).resolve().parent / "results"

# Micromagnetic solver settings
cuda=True
cvode=False

#%% Fixed settings (don't change these)

# Constants of nature
mu0 = 4*np.pi * 1e-7    # Vacuum permeability [N/A^2]
gamma = 2.21*1e5    # Gyromagnetic ratio [m/(A*s)]

# Material properties
Ms = 1.0               # Saturation magnetisation in magnetic region [T]
Ms = Ms/mu0            # [A/m]
alpha = 1              # Gilbert damping [-]
eta = alpha/(1 + alpha**2) * gamma  # Damping constant [m/(A*s)]
K = 0                  # Uniaxial anisotropy constant [J/m^3]
T = 0                  # Temperature [K]
Aex = 1e-11            # Exchange constant in the coupled sections [J/m]
Aex_spacer = 1e-25     # Negligible exchange constant in the spacer cells [J/m]

# Applied field (set to zero)
fct_h_ext = lambda t : np.zeros([1, 3])

# Side lengths of individual micromagnetic cells
a = 1e-8

# Define grid
grid_type = "uniform"
grid_L = [6*a, 6*a, 6*a]            # Grid size along x, y and z for simulated domain
res = [6, 6, 6]                     # Grid resolution, i.e. number of cells along x, y and z
Ntiles = np.prod(res)               # Number of tiles

# Determine which tiles are on corners, edge centers and face centers
corner_n, center_n = np.zeros(Ntiles, dtype=bool), np.zeros(Ntiles, dtype=bool)
xFace_n, yFace_n, zFace_n = np.zeros(Ntiles, dtype=bool), np.zeros(Ntiles, dtype=bool), np.zeros(Ntiles, dtype=bool)
xyEdge_n, xzEdge_n, yzEdge_n = np.zeros(Ntiles, dtype=bool), np.zeros(Ntiles, dtype=bool), np.zeros(Ntiles, dtype=bool)
n = 0
for i in range(6):
    for j in range(6):
        for k in range(6):
            if i in {0, 5}:
                if j in {0, 5}:
                    if k in {0, 5}:
                        corner_n[n] = True
                    elif k in {2,3}:
                        xyEdge_n[n] = True
                elif j in {2,3}:
                    if k in {0, 5}:
                        xzEdge_n[n] = True
                    elif k in {2,3}:
                        xFace_n[n] = True
            elif i in {2,3}:
                if j in {0, 5}:
                    if k in {0, 5}:
                        yzEdge_n[n] = True
                    elif k in {2,3}:
                        yFace_n[n] = True
                elif j in {2,3}:
                    if k in {0, 5}:
                        zFace_n[n] = True
                    elif k in {2,3}:
                        center_n[n] = True
            n += 1

# Identify the coupled magnetic sections and the numerically inert spacers.
nonZeroMagnetisation_n = center_n + corner_n + xFace_n + yFace_n + zFace_n + xyEdge_n + xzEdge_n + yzEdge_n
zeroMagnetisation_n = np.logical_not(nonZeroMagnetisation_n)

# Ms must remain nonzero because MagTense forms exchange and anisotropy factors
# by dividing by Ms.  The tiny spacer exchange decouples the eight sections
# without introducing 0/0 coefficients in the uniform exchange stencil.
Ms_n = Ms * np.ones(Ntiles)
Aex_n = Aex * np.ones((Ntiles, 1))
Aex_n[zeroMagnetisation_n, 0] = Aex_spacer

# Initial magnetisation directions
rng = np.random.default_rng(7)
m0_pv = rng.uniform(-1, 1, [Ntiles, 3])    # Reproducible random directions
m0_pv = m0_pv / np.sqrt(np.sum(m0_pv**2, axis=1, keepdims=True))           # Normalise

#%% Magtense computation

# A section counts as exchange coupled when all its moments agree to within this [-]
error_threshold = 1e-4

sections_s = [center_n, corner_n, xFace_n, yFace_n, zFace_n, xyEdge_n, xzEdge_n, yzEdge_n]
sectionNames_s = ['center', 'corner', 'xFace', 'yFace', 'zFace', 'xyEdge', 'xzEdge', 'yzEdge']


def run_test(plotting: bool = True) -> list[dict]:
    """Run the periodic exchange test and return one check per connected section.

    Each check is a dict with the keys 'check', 'value', 'limit' and 'passed', where the test
    passes when value < limit. This is the contract used by testMagTenseFunctions.py.
    """
    # Setup micromagnetics problem
    problem = MicromagProblem(
        res=res,
        A0=Aex_n,
        Ms=Ms_n[:, np.newaxis],
        K0=K,
        alpha=eta,
        gamma=gamma,
        m0=m0_pv,
        cuda=cuda,
        cvode=cvode,
        grid_L=grid_L,
        grid_type=grid_type,
        usereturnhall=True,
        solver='dynamic',
        T=T,
        usedemag=0,
    )

    # Run simulation
    print('Run simulation')
    result = problem.run_simulation(
        t_end=t_end,
        nt=nTimesteps,
        fct_h_ext=fct_h_ext,
        nt_h_ext=2,  # Two samples completely define the constant zero field.
    )
    print('Done running simulation')

    # Get important data
    M_out = result[1]
    M_pv = M_out[-1, :, 0, :]

    # Test if exchange coupled magnetic moments are indeed identical
    section_errors = []
    for s in range(len(sections_s)):
        section_n = sections_s[s]
        sectionName = sectionNames_s[s]
        Mavg_dv = np.mean(M_pv[section_n, :], axis=0, keepdims=True)
        error_p = np.sqrt(np.sum((M_pv[section_n, :] - Mavg_dv)**2, axis=1))/Ms
        max_error = float(np.max(error_p))
        section_errors.append(max_error)
        if max_error < error_threshold:
            print(f'Exchange coupling works for {sectionName} tiles')
        else:
            print(f'Exchange coupling test failed for {sectionName} tiles')

    if plotting:
        # A logarithmic bar plot shows directly whether every connected section lies
        # below the numerical acceptance threshold.
        bar_colours = [
            'forestgreen' if error < error_threshold else 'firebrick'
            for error in section_errors
        ]
        fig, ax = plt.subplots(layout='constrained', figsize=(9, 4))
        ax.bar(sectionNames_s, section_errors, color=bar_colours)
        ax.axhline(error_threshold, color='black', linestyle='--', label='Pass threshold')
        ax.set_yscale('log')
        ax.set_xlabel('Periodically connected section')
        ax.set_ylabel(r'Maximum $|\mathbf{M}_i - \langle\mathbf{M}\rangle|/M_s$')
        plt.setp(ax.get_xticklabels(), rotation=25, ha='right')
        ax.legend()

        # Save the validation figure beside the other micromagnetic example results
        # and close it so the script does not open an interactive plotting window.
        results_dir.mkdir(parents=True, exist_ok=True)
        figure_path = results_dir / 'periodic_exchange_test.png'
        fig.savefig(figure_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved figure to {figure_path}')

    return [
        {
            'check': f'{name} tiles exchange coupled',
            'value': error,
            'limit': error_threshold,
            'passed': error < error_threshold,
        }
        for name, error in zip(sectionNames_s, section_errors, strict=True)
    ]


if __name__ == '__main__':
    checks = run_test()
    print('periodic_exchange_test '
          + ('PASSED' if all(c['passed'] for c in checks) else 'FAILED'))
