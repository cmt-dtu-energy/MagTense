"""
Test if the shape correction field is implemented correctly in the Magtense micromagnetics solver.

The system is a uniform grid of micromagnetic tiles forming a rectangular prism.
"""

# General modules
import numpy as np
import matplotlib.pyplot as plt

# Magtense stuff
from magtense.micromag import MicromagProblem

# Style settings for plots
plt.rcParams['font.size'] = 16
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

# # Applied field
# Bamp = 1        # Field amplitude [T]
# Bdir_v = np.array([1/np.sqrt(2), 0, 1/np.sqrt(2)])  # Magnetic field direction
# Bdir_v = Bdir_v / np.sqrt(np.sum(Bdir_v**2))        # Normalise
# fct_h_ext = lambda t : Bamp/mu0 * Bdir_v[:, np.newaxis]    # Constant applied field as a 3-by-1 column vector [A/m]
# Applied field (set to zero)
def fct_h_ext(t) -> np.ndarray:
    e_z = np.array([[0, 0, 1]]).T   # Unit vector along z as 3-by-1 column vector
    amplitude = 0    # sine-wave for amplitude
    return np.transpose(amplitude * e_z)    # Has to end with the shape (nt, 3)

# Timestepping
# t_end = 1e-7                        # Total time simulated [s]
t_end = 1e-8                        # Total time simulated [s]
t_step = 1e-12                      # Timestep [s]
nTimesteps = int(t_end / t_step)    # Number of timesteps

# Plot settings
Ndata = 200
dataperiod = int(nTimesteps / Ndata)
show = True

# Micromagnetic solver settings
cuda=False
cvode=True

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
Aex = 1e-10                # Exchange constant [J/m]

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

#%% Magtense computation

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
    grid_L=grid_L,
    grid_type=grid_type,
    usereturnhall=True,
    solver='dynamic',
    T=T,
)
problem.u_ea = u_pv

# Run simulation
print('Run simulation')
result = problem.run_simulation(
    t_end=t_end,
    nt=nTimesteps,
    fct_h_ext=fct_h_ext,
    nt_h_ext=100,
)
print('Done running simulation')


#%% Plot results

# Get important data
t_out, M_out = result[:2]
t_n = t_out
M_npv = M_out[:, :, 0, :]

# Check if the magnetisation is still uniform at the end of the simulation
Mavg_dv = np.mean(M_npv[-1, :, :], axis=0, keepdims=True)
error_p = np.sqrt(np.sum((M_npv[-1, :, :] - Mavg_dv)**2, axis=1))/Ms
if np.max(error_p) < 1e-4:
    print('Magnetisation is uniform')

# Select data to plot
tPlot_n = t_n[::dataperiod] * 1e9   # Switch from s to ns
Mavg_nv = np.mean(M_npv, axis=1)
theta_n = np.arctan2(Mavg_nv[:, 2], Mavg_nv[:, 0])[::dataperiod]

# Make plot
fig, ax = plt.subplots(layout='constrained', figsize=(6, 4))
ax.plot(tPlot_n, theta_n, color='forestgreen')
ax.set_xlabel(r'$\text{Time } [\mathrm{ns}]$')
ax.set_ylabel(r'$\text{Polar angle}$')
ax.set_yticks([0, np.pi/4, np.pi/2])
ax.set_yticklabels(['0', r'$\pi/4$', r'$\pi/2$'])
if show:
    fig.show()