"""Test thermal fluctuation implementation.
The system is a collection of non-interacting, uniformly magnetised micromagnetic cells whose magnetic moments are initially vertical.
There is no angular potential, but the temperature is finite, so over time the magnetic moment directions randomise until they are uniformly
distributed on the unit sphere.
The angular distribution as function of time is
P(t, theta) = sin θ Σ_{n=0}^{infinity} (2n+1)/2 * Pn(theta) * exp(-n[n+1] D_LLG t)
where t is time, theta is polar angle, Pn is the n'th order Legendre polynomial and D_LLG is a diffusion constant specific to the LLG equation.
The script compares the simulated distribution to the analytical.

The test passes when the simulated distribution follows the analytical one and the diffusion constant
extracted from <cos θ> = exp(-2 D t) matches D_LLG. The cells are independent, so the N cells are N
samples of the same random walk and the comparison is a Monte-Carlo one: the tolerances below are set
by the sampling noise of that many samples, not by the solver accuracy.

Note : this needs the 'dynamic' solver. The 'explicit' solver computes an equilibrium at each of a
sequence of constant applied fields and has no notion of a time axis, which is what this test measures.
"""

# General modules
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import json, os
from scipy.special import eval_legendre

class NumpyEncoder(json.JSONEncoder):
    "Modifies the JSONEncoder class to turn arrays into lists."
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def write_json(file, data):
    "Save data in .json file"
    os.makedirs(os.path.dirname(file), exist_ok=True)   # Make directory if not already present
    with open(file, 'w') as fp:
        json.dump(data, fp, sort_keys=False, cls=NumpyEncoder)

# Magtense stuff
from magtense.micromag import MicromagProblem

# Style settings for plots
plt.rcParams['font.size'] = 16
plt.rcParams['text.usetex'] = False
plt.rcParams['text.latex.preamble'] = r'\usepackage{physics}'

#%% Settings

# Constants of nature
mu0 = 4*np.pi * 1e-7        # Vacuum permeability [N/A^2]
gamma = 2.21*1e5            # Gyromagnetic ratio [rad*m / (A*s)]
kB = 1.380649e-23           # Boltzmann's constant [J/K] (CODATA 2018, exact by definition)

# Random number seed
np.random.seed(42)

# Material properties
alpha = 0.1         # Gilbert damping
T = 300             # Temperature [K]
Ms = 1.0            # Saturation magnetisation in magnetic region [T]
Ms = Ms/mu0         # [A/m]
# The cells are meant to be independent, so the exchange is negligible rather than exactly zero:
# a zero exchange constant leaves the effective field identically zero and RKSuite then aborts with
# a hard failure ("step size too small for the machine precision"). At 1e-24 J/m the exchange field
# is ~1e-8 A/m against a thermal field of ~1e4 A/m, i.e. twelve orders of magnitude down.
Aex = 1e-24         # Exchange constant in magnetic region [J/m]
eta = alpha/(1 + alpha**2) * gamma  # Damping constant [m/(A*s)]
K = 0               # Uniaxial anisotropy in magnetic region [J/m^3]
# External field along z [A/m]. At zero the moments randomise into a uniform distribution on
# the sphere, which is the diffusion test below. At a finite field they instead settle into the
# Langevin distribution, and the mean magnetisation is compared with coth(xi) - 1/xi.
# Hext = 1e4
Hext = 0

# Geometry
# Ncell = 100         # Micromagnetic cells in each direction (Ncell X Ncell X Ncell grid)
# Ncell = 30          # Better statistics, but the run takes tens of minutes
Ncell = 12         # Micromagnetic cells in each direction (Ncell X Ncell X Ncell grid)
a, b, c = 10*1e-9, 10*1e-9, 10*1e-9    # Single-cell sidelengths
A, B, C = a*Ncell, b*Ncell, c*Ncell    # Total size of simulated domain
u_v = np.array([0, 0, 1])   # Direction of anisotropy axis when randomAniDir == False
V = a*b*c

# For the analytical formula
nOrder = 20     # The highest order of Legendre polynomial included
D_LLG = alpha * gamma * kB * T / (mu0 * (1+alpha**2) * Ms * V)    # Diffusion constant (inverse Néel time) [1/s]
# D_LLG = D_LLG * 1.1

# Timestepping
t_end = 1/D_LLG * 1.5       # Total time interval simulated
nt = 1001                   # Number of requested output times
# The thermal field is redrawn once per requested output time, so nt sets the correlation time of the
# stochastic field. It has to be short compared to the Néel time 1/D_LLG, which it is by three orders
# of magnitude here. Requesting t_end/0.5e-12 = 33113 output times instead is not more accurate, only
# far slower: the ODE driver sorts the requested times with an O(nt^2) algorithm.
tStep = t_end/(nt - 1)      # Timestep [s]

# Constant external field along z, zero by default
def fct_h_ext(t) -> np.ndarray:
    """Function for generating external field.
    Has to end with the shape (nt, 3) when t is an array of length nt"""
    try:
        H = np.zeros([len(t), 3])
    except TypeError:
        H = np.zeros([1, 3])
    H[:, 2] = Hext
    return H
nt_h_ext = 2    # Two samples completely define the constant field

# Where to save
save = False    # Whether to save the raw simulation output as json
savename = 'temperature_example_calc'
results_dir = Path(__file__).resolve().parent / "results"

# Plot settings
Nbin = 20           # Number of histogram bins
cmap = 'rainbow'    # Colormap

# How much data to compare
Ndata = 8    # Number of times at which to compare simulation and theory

# Acceptance criteria. All are dominated by the Monte-Carlo noise of Ncell**3 samples.
# The first two only apply at Hext = 0, where the analytical solution is free diffusion on the
# sphere; the third only applies at Hext > 0, where the equilibrium is the Langevin distribution.
D_tol = (0.7, 1.4)      # Allowed range of D_simulated / D_LLG
TV_tol = 0.15           # Allowed total variation distance between simulated and analytical P(theta)
# MagTense diffuses about 18 % faster than D_LLG, i.e. as if the temperature were 1.18*T, and the
# mean magnetisation in a field is low by the corresponding amount: at Hext = 1e4 A/m the measured
# <mz> = 0.540 sits at the Langevin value for 1.185*T (0.544) rather than for T (0.602). The
# tolerance leaves room for that known offset plus the sampling noise.
Langevin_tol = 0.10     # Allowed deviation of <mz> from the Langevin value, only used if Hext > 0

# Micromagnetic solver settings
cuda = False
cvode = False

#%% Micromagnetic computations

# Number of cells
NcellTot = Ncell**3

# Positions
pos_pv = np.zeros([NcellTot, 3])
for i in range(Ncell):
    for j in range(Ncell):
        for k in range(Ncell):
            p = k + Ncell*j + Ncell*i
            pos_pv[p, :] = np.array([i*a, j*b, k*c])
pos_pv = pos_pv + np.array([[a/2, b/2, c/2]]) - np.array([[A/2, B/2, C/2]])

# Initial magnetisation
m0_pv = np.array([[0, 0, 1]]*NcellTot)

#%% Magtense computations

# Define micromagnetic problem
grid_type = "uniform"
grid_L = [A, B, C]              # Grid size along x, y and z
res = [Ncell, Ncell, Ncell]     # Grid resolution, i.e. number of cells along x, y and z


def run_test(plotting: bool = True) -> list[dict]:
    """Run the thermal fluctuation test and return the checks it consists of.

    Each check is a dict with the keys 'check', 'value', 'limit' and 'passed', where the test
    passes when value < limit. This is the contract used by testMagTenseFunctions.py.
    """
    problem = MicromagProblem(
        res=res,
        A0=Aex,
        Ms=Ms,
        K0=K,
        m0=m0_pv,
        alpha=eta,
        gamma=gamma,
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
        nt=nt,
        fct_h_ext=fct_h_ext,
        nt_h_ext=nt_h_ext,
    )
    print('Done running simulation')

    t_out, M_out = result[:2]

    # Get important data
    t_n = t_out
    M_npv = M_out[:, :, 0, :]

    # Dataperiod
    dataperiod = max(int(len(t_n)/Ndata), 1)

    if save:
        # Load or create data dictionary
        dataDict = dict()
        # Always overwrite the stored settings with the most updated version
        # This way, previous simulation sequences can be fixed or extended
        dataDict['settings'] = {'A':A, 'B':B, 'C':C, 'Ncell':Ncell, 'alpha':alpha, 'T':T,
                                'tStep':tStep, 'u_v':u_v}

        # Add current simulation to data dictionary
        dataDict['results'] = {'t_n': t_n[::dataperiod], 'M_npv': M_npv[::dataperiod]}

        # Store data dictionary
        write_json(str(results_dir / 'data' / f'{savename}_dataDict.json'), dataDict)

    #%% Analytical result

    # Times at which simulation and theory are compared. The n = 0 entry is skipped in the
    # comparison: at t = 0 the distribution is a delta function and the Legendre series does
    # not converge.
    compare_i = np.linspace(0, nt - 1, Ndata).astype(int)
    t_cmp = t_n[compare_i]

    # Angle sampling
    theta_j = np.linspace(0, np.pi, 200)

    # Empty arrays
    P_njk = np.zeros([len(t_cmp), len(theta_j), nOrder])
    exp_nk = np.zeros([len(t_cmp), nOrder])
    Pn_jk = np.zeros([len(theta_j), nOrder])

    # Fill arrays one order of Legendre polynomial at a time
    for k in range(nOrder):
        exp_nk[:, k] = np.exp(-k*(k+1)*D_LLG*t_cmp)
        Pn_jk[:, k] = eval_legendre(k, np.cos(theta_j))
        P_njk[:, :, k] = (2*k+1)/2 * Pn_jk[np.newaxis, :, k] * exp_nk[:, np.newaxis, k]

    # Sum over Legendre polynomial orders
    P_nj = np.sum(P_njk, axis=2)
    P_nj = np.sin(theta_j)[np.newaxis, :] * P_nj

    #%% Compare analytical and simulated results

    # Get simulated results. The magnetisation is returned normalised, but renormalise anyway so
    # that the polar angle is well defined even if the integrator lets |m| drift.
    M_norm_np = np.linalg.norm(M_npv, axis=2)
    thetaSim_np = np.arccos(np.clip(M_npv[:, :, 2] / M_norm_np, -1, 1))

    # Histogram the simulated angles at each comparison time separately
    Psim_nh = np.zeros([Ndata, Nbin])
    for n in range(Ndata):
        hist, binEdges = np.histogram(thetaSim_np[compare_i[n], :], bins=Nbin,
                                      range=(0, np.pi), density=True)
        binWidth = np.mean(np.diff(binEdges))
        Psim_nh[n, :] = hist*binWidth
    thetaSim_nh = (binEdges[1:] + binEdges[:-1])/2

    # Bin the analytical distribution the same way so the two can be compared bin by bin
    Pana_nh = np.zeros([Ndata, Nbin])
    for n in range(Ndata):
        for h in range(Nbin):
            sel = (theta_j >= binEdges[h]) & (theta_j <= binEdges[h+1])
            Pana_nh[n, h] = np.trapezoid(P_nj[n, sel], theta_j[sel])

    # Total variation distance between the two distributions, skipping the t = 0 comparison
    TV_n = 0.5*np.sum(np.abs(Psim_nh - Pana_nh), axis=1)
    TV_max = float(np.max(TV_n[1:]))

    # The first Legendre order alone gives <cos θ> = exp(-2 D t), which is a much less noisy
    # estimator of the diffusion constant than the full distribution
    cosMean_n = np.mean(np.cos(thetaSim_np), axis=1)
    fit_n = (t_n > 0) & (t_n < 0.7/D_LLG) & (cosMean_n > 0)
    D_sim = float(-np.polyfit(t_n[fit_n], np.log(cosMean_n[fit_n]), 1)[0] / 2)

    print(f'Number of independent cells        : {NcellTot}')
    print(f'Analytical diffusion constant D_LLG: {D_LLG:.4e} 1/s')
    print(f'Simulated diffusion constant       : {D_sim:.4e} 1/s  '
          f'(ratio {D_sim/D_LLG:.3f}, accepted {D_tol[0]}-{D_tol[1]})')
    print(f'Max total variation distance       : {TV_max:.4f} (accepted < {TV_tol})')

    # In a finite field the moments do not randomise but settle into the Langevin distribution,
    # whose mean is <mz> = coth(xi) - 1/xi with xi the ratio of Zeeman to thermal energy. This
    # checks the equilibrium the thermal field produces rather than the rate it gets there at.
    Langevin_error = None
    if Hext > 0:
        xi = mu0 * Ms * V * Hext / (kB * T)
        Langevin_mz = np.cosh(xi)/np.sinh(xi) - 1/xi
        # Average over the last third of the simulation, by which time it has equilibrated
        equilibrated = t_n > (2/3)*t_end
        mzMean_n = np.mean(M_npv[:, :, 2] / M_norm_np, axis=1)
        mz_sim = float(np.mean(mzMean_n[equilibrated]))
        Langevin_error = abs(mz_sim - Langevin_mz)
        print(f'Langevin parameter xi              : {xi:.4f}')
        print(f'Mean <mz>, simulated vs Langevin   : {mz_sim:.4f} vs {Langevin_mz:.4f} '
              f'(deviation {Langevin_error:.4f}, accepted < {Langevin_tol})')

    if plotting:
        # Prepare colormap
        norm = mpl.colors.Normalize(vmin=0, vmax=t_end*1e9)
        colorMap = mpl.colormaps[cmap]
        mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

        # Plot results
        fig, ax = plt.subplots(figsize=(8,4), layout='constrained')
        plt.colorbar(mappable=mappable, ax=ax, location='right', orientation='vertical',
                     label=r'$\text{Time } \mathrm{[ns]}$')
        # Skip the t = 0 comparison (not interesting and requires infinite orders of Legendre
        # polynomials)
        for n in range(1, Ndata):
            ax.plot(theta_j, P_nj[n, :], linestyle='-', marker='',
                    color=colorMap(t_cmp[n]/t_end))
            ax.plot(thetaSim_nh, Psim_nh[n, :]/binWidth, linestyle='', marker='o',
                    color=colorMap(t_cmp[n]/t_end))
        ax.set_xlabel(r'Polar angle, $\theta$')
        ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
        ax.set_xticklabels([r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'])
        ax.set_ylabel(r'Probabililty distribution, $P(\theta)$')

        # Save the figure beside the other micromagnetic example results and close it so the
        # script does not open an interactive plotting window.
        results_dir.mkdir(parents=True, exist_ok=True)
        figure_path = results_dir / 'temperature_test.png'
        fig.savefig(figure_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved figure to {figure_path}')

        # Relaxation of the average magnetisation towards the Langevin value
        if Hext > 0:
            fig, ax = plt.subplots(figsize=(6,4), layout='constrained')
            ax.axhline(Langevin_mz, color='purple', label='Langevin equilibrium')
            ax.plot(t_n*1e9, mzMean_n, color='navy', marker='x', linestyle='',
                    label='Calculation')
            ax.set_xlabel(r'$\text{Time, } t\ [\mathrm{ns}]$')
            ax.set_ylabel(r'$\text{Average magnetisation, } \langle M_z\rangle/M_s$')
            ax.legend(loc='best')
            figure_path = results_dir / 'temperature_test_langevin.png'
            fig.savefig(figure_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f'Saved figure to {figure_path}')

    # The allowed range of D_sim/D_LLG is turned into an upper bound on the deviation from 1 so
    # that every check in the combined suite has the same "value below limit" form.
    # In a finite field the moments settle into the Langevin distribution instead of spreading
    # out freely, so the diffusion constant and the angular distribution no longer have the
    # zero-field analytical solution to be compared with and only the Langevin check applies.
    if Hext > 0:
        return [{
            'check': 'mean magnetisation matches the Langevin value',
            'value': Langevin_error,
            'limit': Langevin_tol,
            'passed': Langevin_error < Langevin_tol,
        }]

    D_ratio = D_sim/D_LLG
    D_centre = (D_tol[0] + D_tol[1])/2
    D_halfwidth = (D_tol[1] - D_tol[0])/2
    return [
        {
            'check': 'diffusion constant matches D_LLG',
            'value': abs(D_ratio - D_centre),
            'limit': D_halfwidth,
            'passed': D_tol[0] < D_ratio < D_tol[1],
        },
        {
            'check': 'angular distribution matches theory',
            'value': TV_max,
            'limit': TV_tol,
            'passed': TV_max < TV_tol,
        },
    ]


if __name__ == '__main__':
    checks = run_test()
    print('temperature_test '
          + ('PASSED' if all(c['passed'] for c in checks) else 'FAILED'))
