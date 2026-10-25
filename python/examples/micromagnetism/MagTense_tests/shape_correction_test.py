"""
Test if the shape correction field is implemented correctly in the Magtense micromagnetics solver.

The system is a uniform grid of micromagnetic tiles forming a rectangular prism.

A second check uses an unstructured prism mesh running from 0 to L whose first layer of cells is
refined, so the cells at the two ends differ in size. The macrogeometry prism of the shape
correction has to be centred on the outer faces of that mesh. For a uniform magnetisation the
field of the cells and that of the macrogeometry prism then cancel exactly, and what is left is
the field of the sample. Centring it on the midpoint of the end-cell centres instead shifts it by
a/8 and leaves a field error of about 20 %.

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
plt.rcParams['font.size'] = 15
plt.rcParams['text.usetex'] = False
plt.rcParams['text.latex.preamble'] = r'\usepackage{physics}'

#%% Settings

# Define grid
a = 1e-9    # Side lengths of individual micromagnetic cell [m]
grid_type = "uniform"
res = [2, 3, 12]                    # Grid resolution, i.e. number of cells along x, y and z
grid_L = [res[0]*a, res[1]*a, res[2]*a]            # Grid size along x, y and z for simulated domain
# res = [1, 1, 1]                    # Grid resolution, i.e. number of cells along x, y and z
# grid_L = [res[0]*a, res[1]*a, res[2]*a]            # Grid size along x, y and z for simulated domain
Ntiles = np.prod(res)   # Number of tiles

# Sample- and macrogeometry shapes (assumes rectangular prisms)
# The sample is intentionally macroscopic so the demagnetisation field becomes uniform across the simulated domain.
macroShape = np.array(grid_L)
aSample = 3  # Length of sample along x
# aSample = 6  # Length of sample along x
bSample = 1  # Transverse dimension of sample
sampleShape = np.array([aSample, bSample, bSample])

# Applied field (set to zero)
def fct_h_ext(t) -> np.ndarray:
    e_z = np.array([[0, 0, 1]]).T   # Unit vector along z as 3-by-1 column vector
    amplitude = 0    # sine-wave for amplitude
    return np.transpose(amplitude * e_z)    # Has to end with the shape (nt, 3)

# Timestepping
# t_end = 1e-7                        # Total time simulated [s]
t_end = 1e-8                        # Total time simulated [s]
t_step = 5e-11                      # Requested output spacing [s]
nTimesteps = int(round(t_end / t_step)) + 1  # Include both time endpoints

# Plot settings
Ndata = 200
dataperiod = int(nTimesteps / Ndata)
output_dir = Path(__file__).resolve().parent

# How far the magnetisation may deviate from uniform to count as passed [-]
uniformity_tol = 1e-4
# How far the mean polar angle may drift over the simulation to count as passed [rad].
# With the shape correction in place the demagnetisation field of the sample is cancelled
# exactly by the uniaxial anisotropy, so theta should not move at all.
drift_tol = 1e-5

# Micromagnetic solver settings
cuda=False
cvode=False

#%% Fixed settings (don't change these)

# Constants of nature
mu0 = 4*np.pi * 1e-7    # Vacuum permeability [N/A^2]
# gamma = 2.21*1e5    # Gyromagnetic ratio [m/(A*s)]
gamma = 0    # Gyromagnetic ratio [m/(A*s)]

# Material properties
Ms = 0.1               # Saturation magnetisation in magnetic region [T]
Ms = Ms/mu0            # [A/m]
alpha = 1              # Gilbert damping [-]
# eta = alpha/(1 + alpha**2) * gamma  # Damping constant [m/(A*s)]
eta = alpha/(1 + alpha**2) * 2.21*1e5  # Damping constant [m/(A*s)]
T = 0                  # Temperature [K]
# A moderate exchange retains a uniform state without making RKSuite stiff on
# the 1 nm mesh.  It also avoids the zero-exchange 0/0 normalisation case.
Aex = 1e-14            # Exchange constant [J/m]

# Initial magnetisation directions
m0_v = np.array([1/np.sqrt(2), 0, 1/np.sqrt(2)])    # xz-plane at 45 degrees to horizontal
# m0_v = np.array([2/np.sqrt(2), 0, 1/np.sqrt(2)])    # xz-plane at 45 degrees to horizontal
# m0_v = np.random.uniform(-1, 1, [Ntiles, 3])   # Uniformly random direction
m0_v = m0_v / np.sqrt(np.sum(m0_v**2))           # Normalise
m0_pv = np.tile(m0_v, [Ntiles, 1])          # Give alle tiles the same initial magnetisation

#%% Magnetocrystalline anisotropy

def Kshape(a, b):
    """The shape dependent, uniaxial anisotropy constant from the whole samples demagnetisation field.
    Assumes equal length along y and z, i.e. b=c"""
    fx = b**2 / (a * np.sqrt(a**2 + 2*b**2))
    Nxx = 2/np.pi * np.arctan(fx)
    return mu0*Ms**2/4 * (1 - 3*Nxx)
Ksh = Kshape(aSample, bSample)

# Uniaxial anisotropy constant (equated to shape anisotropy) [J/m^3]
K = -Ksh

# Magnetocrystalline anisotropy axis
# u_pv = np.tile(np.array([0, 0, 1]), [Ntiles, 1])
u_pv = np.tile(np.array([1, 0, 0]), [Ntiles, 1])

#%% Unstructured mesh with unequal end cells

# The shape correction places the macrogeometry prism on the centre of the simulated domain. On an
# unstructured mesh that centre has to be taken from the outer faces of the cells: the centres of
# the end cells only give it when the cells at the two ends have the same size. This mesh is a
# 4 x 2 x 2 block of 1 nm cubes running from 0 to L, whose first layer along x is refined into
# cells of half the size, so the centres of the end cells are a/4 and a/2 from the two faces.
mesh_base = (4, 2, 2)           # Base cells along x, y and z
mesh_grid_L = np.array(mesh_base) * a
# For a uniform magnetisation the field of the cells and the macrogeometry prism cancel exactly,
# by superposition, wherever the prism sits on the domain. What is left is the field of the sample,
# so it has to be uniform to the single precision the tensor is stored in.
mesh_field_tol = 1e-5


def refined_end_mesh() -> tuple[np.ndarray, np.ndarray]:
    """Cell centres and side lengths of the mesh with a refined first layer along x."""
    pts, abc = [], []
    for k in range(mesh_base[2]):
        for j in range(mesh_base[1]):
            for i in range(mesh_base[0]):
                corner = np.array([i, j, k]) * a
                if i == 0:
                    for sub in np.ndindex(2, 2, 2):
                        pts.append(corner + (np.array(sub) + 0.5) * a / 2)
                        abc.append(np.full(3, a / 2))
                else:
                    pts.append(corner + 0.5 * a)
                    abc.append(np.full(3, a))
    return np.array(pts), np.array(abc)


def sample_demag_tensor() -> np.ndarray:
    """Demagnetisation tensor at the centre of the sample prism, which has b = c."""
    Nxx_s = 2/np.pi * np.arctan(bSample**2 / (aSample * np.sqrt(aSample**2 + 2*bSample**2)))
    Nyy_s = (1 - Nxx_s) / 2
    return np.diag([Nxx_s, Nyy_s, Nyy_s])


def run_mesh_check() -> float:
    """Relative deviation of the demagnetisation field from that of the sample, on the mesh."""
    pts, abc = refined_end_mesh()
    ntot = len(pts)
    problem = MicromagProblem(
        res=(ntot, 1, 1),
        grid_type='unstructuredPrisms',
        grid_pts=pts,
        grid_abc=abc,
        grid_L=list(mesh_grid_L),
        A0=Aex * np.ones((ntot, 1)),
        Ms=Ms * np.ones((ntot, 1)),
        K0=0,
        alpha=eta,
        gamma=0,
        m0=np.tile(m0_v, [ntot, 1]),
        macroShape=mesh_grid_L,
        sampleShape=sampleShape,
        cuda=cuda,
        cvode=cvode,
        # The shape correction is evaluated at the cell centres, so the cells have to be as well
        # for the field of the cells and of the macrogeometry prism to cancel
        useavgn=False,
        usereturnhall=True,
        solver='dynamic',
        exch_presize=64,
    )
    problem.use_fmm = 0
    # Only the field of the initial state is needed
    result = problem.run_simulation(
        t_end=1e-15, nt=2, fct_h_ext=fct_h_ext, nt_h_ext=2,
    )
    H_dem = np.asarray(result[5][0, :, 0, :])
    H_sample = -Ms * sample_demag_tensor() @ m0_v
    return float(np.abs(H_dem - H_sample).max() / np.abs(H_sample).max())


#%% Magtense computation

def run_test(plotting: bool = True) -> list[dict]:
    """Run the shape correction test and return the checks it consists of.

    Each check is a dict with the keys 'check', 'value', 'limit' and 'passed', where the test
    passes when value < limit. This is the contract used by testMagTenseFunctions.py.
    """
    problem = MicromagProblem(
        res=res,
        A0=Aex,
        Ms=Ms,
        K0=K,
        alpha=eta,
        gamma=gamma,
        m0=m0_pv,
        macroShape=macroShape,
        sampleShape=sampleShape,
        cuda=cuda,
        cvode=cvode,
        # The uniaxial anisotropy below is set from the analytical demagnetisation factor of the
        # sample prism, which is what the tensor evaluated at the cell centres reproduces. With
        # the averaged tensor (the library default) the cancellation is only good to 1.1e-4 rad
        # instead of 2.6e-8, so this test opts out explicitly to keep the tight drift limit.
        useavgn=False,
        grid_L=grid_L,
        grid_type=grid_type,
        usereturnhall=True,
        solver='dynamic',
        T=T,
    )
    # The FMM demag path segfaults on this small grid when the library is built with
    # USE_FMM3D=1, and an octree is pointless for 72 tiles anyway. std_problem_3.py
    # disables it the same way.
    problem.use_fmm = 0

    problem.u_ea = u_pv

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
    t_out, M_out = result[:2]
    t_n = t_out
    M_npv = M_out[:, :, 0, :]

    # Check if the magnetisation is still uniform at the end of the simulation
    Mavg_dv = np.mean(M_npv[-1, :, :], axis=0, keepdims=True)
    error_p = np.sqrt(np.sum((M_npv[-1, :, :] - Mavg_dv)**2, axis=1))/Ms
    uniformity = float(np.max(error_p))
    if uniformity < uniformity_tol:
        print('Magnetisation is uniform')
    else:
        print('Magnetisation is NOT uniform')

    # The shape correction cancels the demagnetisation field of the sample against the uniaxial
    # anisotropy, so the mean polar angle should not move over the simulation at all
    tPlot_n = t_n[::dataperiod] * 1e9   # Switch from s to ns
    Mavg_nv = np.mean(M_npv, axis=1)
    theta_n = np.arctan2(Mavg_nv[:, 2], Mavg_nv[:, 0])[::dataperiod]
    drift = float(np.max(np.abs(theta_n - theta_n[0])))
    verdict = 'stays put' if drift < drift_tol else 'DRIFTS'
    print(f'Mean polar angle {verdict}: max change = {drift:.3e} rad')

    if plotting:
        # Make plot
        fig, ax = plt.subplots(layout='constrained', figsize=(6, 4))
        ax.plot(tPlot_n, theta_n-theta_n[0], color='crimson',
                label=r'Calculation (should be constant)')
        ax.set_xlabel(r'$\text{Time } [\mathrm{ns}]$')
        ax.set_ylabel(r'$\text{Change in polar angle}$')
        angle_window = 10**(np.floor(np.log10(drift_tol*10)))
        ax.set_yticks([-angle_window, 0, angle_window])
        ax.set_ylim(-angle_window, angle_window)
        ax.set_yticklabels([f'$-10^{{{np.log10(angle_window):.0f}}}$', r'$0$',
                            f'$10^{{{np.log10(angle_window):.0f}}}$'])
        ax.fill_between(np.array([tPlot_n[0], tPlot_n[-1]]),
                        np.array([-drift_tol, -drift_tol]), np.array([drift_tol, drift_tol]),
                        alpha=0.2, color='forestgreen', label='Accepted interval')
        ax.legend(loc='best')

        # Save the validation figure beside this script and close it so the script
        # does not open an interactive plotting window.
        figure_path = output_dir / "shape_correction_test.png"
        fig.savefig(figure_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved figure to {figure_path}")

    mesh_deviation = run_mesh_check()
    verdict = 'is the sample field' if mesh_deviation < mesh_field_tol else 'is NOT the sample field'
    print(f'Unstructured mesh with unequal end cells: the field {verdict}, '
          f'max deviation = {mesh_deviation:.3e}')

    return [
        {
            'check': 'magnetisation stays uniform',
            'value': uniformity,
            'limit': uniformity_tol,
            'passed': uniformity < uniformity_tol,
        },
        {
            'check': 'mean polar angle does not drift',
            'value': drift,
            'limit': drift_tol,
            'passed': drift < drift_tol,
        },
        {
            'check': 'unstructured mesh, unequal end cells: field is the sample field',
            'value': mesh_deviation,
            'limit': mesh_field_tol,
            'passed': mesh_deviation < mesh_field_tol,
        },
    ]


if __name__ == '__main__':
    checks = run_test()
    print('shape_correction_test '
          + ('PASSED' if all(c['passed'] for c in checks) else 'FAILED'))
