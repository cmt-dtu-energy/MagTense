"""Test thermal fluctuation implementation.
The system is a collection of non-interacting, uniformly magnetised micromagnetic cells whose magnetic moments are initially vertical.
There is no angular potential, but the temperature is finite, so over time the magnetic moment directions randomise until they are uniformly
distributed on the unit sphere.
The angular distribution as function of time is
P(t, theta) = sin θ Σ_{n=0}^{infinity} (2n+1)/2 * Pn(theta) * exp(-n[n+1] D_LLG t)
where t is time, theta is polar angle, Pn is the n'th order Legendre polynomial and D_LLG is a diffusion constant specific to the LLG equation.
The script compares the simulated distribution to the analytical."""

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
kB = 1.3806488 * 10**(-23)  # Boltzmann's constant [J/K]

# Random number seed
np.random.seed(42)

# Material properties
alpha = 0.1         # Gilbert damping
T = 300             # Temperature [K]
Ms = 1.0            # Saturation magnetisation in magnetic region [T]
Ms = Ms/mu0         # [A/m]
Aex = 1e-24         # Exchange constant in magnetic region [J/m]
eta = alpha/(1 + alpha**2) * gamma  # Damping constant [m/(A*s)]
K = 0               # Uniaxial anisotropy in magnetic region [J/m^3]
# Hext = 1e4          # External field along z [A/m]
Hext = 0          # External field along z [A/m]

# Geometry
Ncell = 30         # Micromagnetic cells in each direction (Ncell X Ncell X Ncell grid)
# Ncell = 20        # Micromagnetic cells in each direction (Ncell X Ncell X Ncell grid)
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
nt = 1001                   # Number of timesteps
tStep = t_end/(nt - 1)      # Timestep [s]
# tStep = 1*1e-10           # Timestep [s]
# nt = round(t_end/tStep)     # Number of timesteps

# Zero external field
def fct_h_ext(t) -> np.ndarray:
    """Function for generating external field.
    Has to end with the shape (nt, 3) when t is an array of length nt"""
    try:
        H = np.zeros([len(t), 3])
    except TypeError:
        H = np.zeros([1, 3])
    H[:, 2] = Hext
    return H
nt_h_ext = 2    # Number of time coordinates at which to evaluate fct_h_ext

# Where to save
save = False    # Whether to save the raw simulation output as json
savename = 'temperature_example_calc'
results_dir = Path(__file__).resolve().parent / "results"

# Plot settings
show = True         # Whether to show the plots
Nbins = 20          # Number of histogram bins
cmap = 'rainbow'    # Colormap
colors_n = ['navy', 'forestgreen', 'crimson', 'magenta', 'darkorange', 'purple', 'darkcyan', 'darkgoldenrod']

# How much data to save
Ndata = 6    # Number of datapoints to save

# Micromagnetic solver settings
cuda=True
cvode=False

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

problem = MicromagProblem(
    res=res,
    A0=Aex,
    Ms=Ms,
    K0=K,
    m0=m0_pv,
    alpha=eta,
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
pts = result[2]
H_exc, H_ext, H_dem, H_ani = result[3:7]

#%% Store results

# Get important data
tSim_n = t_out
M_npv = M_out[:, :, 0, :]
Bext_nv = H_ext[:, 0, 0, :] * mu0

# Data selector
dataSelector = np.linspace(0, nt-1, Ndata).astype(int)

if save:
    # Load or create data dictionary
    dataDict = dict()
    # Always overwrite the stored settings with the most updated version
    # This way, previous simulation sequences can be fixed or extended
    dataDict['settings'] = {'A':A, 'B':B, 'C':C, 'Ncell':Ncell, 'alpha':alpha, 'T':T, 'tStep':tStep, 'u_v':u_v}

    # Add current simulation to data dictionary
    dataDict[f'results'] = {'t_n': tSim_n[dataSelector], 'M_npv': M_npv[dataSelector]}

    # Store data dictionary
    write_json(f'data/{savename}_dataDict.json', dataDict)


#%% Analytical result

# Time and angle sampling
t_n = np.linspace(0, t_end, Ndata)
theta_j = np.linspace(0, np.pi, 200)

# Empty arrays
P_njk = np.zeros([len(t_n), len(theta_j), nOrder])
exp_nk = np.zeros([len(t_n), nOrder])
Pn_jk = np.zeros([len(theta_j), nOrder])

# Fill arrays one order of Legendre polynomial at a time
for k in range(nOrder):
    exp_nk[:, k] = np.exp(-k*(k+1)*D_LLG*t_n)
    Pn_jk[:, k] = eval_legendre(k, np.cos(theta_j))
    P_njk[:, :, k] = (2*k+1)/2 * Pn_jk[np.newaxis, :, k] * exp_nk[:, np.newaxis, k]

# Sum over Legendre polynomial orders
P_nj = np.sum(P_njk, axis=2)
P_nj = np.sin(theta_j)[np.newaxis, :] * P_nj

#%% Compare analytical and simulated results

# Get simulated results
thetaSim_np = np.arccos(M_npv[dataSelector, :, 2])
Psim_nh = np.zeros([Ndata, Nbins])
for n in range(Ndata):
    hist, binEdges = np.histogram(thetaSim_np[n, :], bins=Nbins, range=(0, np.pi), density=True)
    binWidth = np.mean(np.diff(binEdges))
    # Psim_nh[n, :] = hist*binWidth
    Psim_nh[n, :] = hist
thetaSim_nh = (binEdges[1:] + binEdges[:-1])/2

# Prepare colormap
norm = mpl.colors.Normalize(vmin=0, vmax=t_end*1e9)
colorMap = mpl.colormaps[cmap]
mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

# Plot results
fig, ax = plt.subplots(figsize=(8,4), layout='constrained')
plt.colorbar(mappable=mappable, ax=ax, location = 'right', orientation = 'vertical', label=r'$\text{Time } \mathrm{[ns]}$')
for n in range(1, Ndata):   # Skip the t = 0 comparison (not interesting and requires infinite orders of Legendre polynomials)
    ax.plot(theta_j, P_nj[n, :], linestyle='-', marker='', color=colorMap(t_n[n]/t_end))
    ax.plot(thetaSim_nh, Psim_nh[n, :], linestyle='', marker='o', color=colorMap(t_n[n]/t_end))
ax.set_xlabel(r'Polar angle, $\theta$')
ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
ax.set_xticklabels([r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'])
ax.set_ylabel(r'Probabililty distribution, $P(\theta)$')

# Save the figure beside the other micromagnetic example results and close it so the script
# does not open an interactive plotting window.
results_dir.mkdir(parents=True, exist_ok=True)
figure_path = results_dir / 'temperature_test.png'
fig.savefig(figure_path, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'Saved figure to {figure_path}')


# Plot relaxation of average magnetisation
if Hext > 0:
    Mavg_n = np.mean(M_npv[:, :, 2], axis=1)
    xi = mu0 * Ms * V * Hext / (kB * T)
    fig, ax = plt.subplots(figsize=(6,4), layout='constrained')
    ax.axhline(np.cosh(xi)/np.sinh(xi) - 1/xi, color='purple')
    ax.plot(tSim_n*1e9, Mavg_n, color='navy', marker='x')
    ax.set_xlabel(r'$\text{Time, } t [ns]$')
    ax.set_ylabel(r'$\text{Average magnetisation, } <M_z>/M_s$')
    fig.show()