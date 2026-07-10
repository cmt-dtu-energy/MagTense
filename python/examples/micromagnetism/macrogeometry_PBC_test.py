"""
Test if periodic boundary conditions as modeled by the macrogeometry method are implemented correctly in the Magtense micromagnetics solver.
SHORT EXPLANATION : If the plotted polar angle is constant, the code works.

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

Other note : By changing testAxis, the entire setup is rotated so the axis-of-periodicity becomes y or z, allowing
independent validation along x, y and z.
"""

# General modules
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
A = 1.1837e-7   # Critical point where the torque is zero at all theta (constant theta)
# A = 2e-7      # Should go towards theta = π/2 (increasing theta)
# A = 0.7e-7    # Should go towards theta = 0 (decreasing theta)

# Which axis to test PBC
testAxis = 0    # 0 for x-axis, 1 for y-axis and 2 for z-axis

# Timestepping
t_end = 1e-7                        # Total time simulated [s]
# t_end = 1e-8                      # Total time simulated [s]
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
gamma = 0    # Gyromagnetic ratio set to zero to turn off precession so motion stays 2D

# Material properties
Ms = 1.0               # Saturation magnetisation in magnetic region [T]
Ms = Ms/mu0            # [A/m]
alpha = 1              # Gilbert damping [-]
eta = alpha/(1 + alpha**2) * 2.21*1e5  # Damping constant [m/(A*s)]
K = 2.7e4              # Uniaxial anisotropy constant [J/m^3]
T = 0                  # Temperature [K]
Aex = 0                # Exchange constant [J/m]

# Applied field
def fct_h_ext(t) -> np.ndarray:
    e_z = np.array([[0, 0, 1]]).T   # Unit vector along z as 3-by-1 column vector
    amplitude = 0    # sine-wave for amplitude
    return np.transpose(amplitude * e_z)    # Has to end with the shape (nt, 3)

# Side lengths of individual micromagnetic cells
a = 1e-8    # [m]

# Initial magnetisations
theta_p = np.array([np.pi/4, -np.pi/4])
if testAxis == 0:
    mx_p = np.cos(theta_p)
    my_p = np.array([0, 0])
    mz_p = np.sin(theta_p)
elif testAxis == 1:
    mx_p = np.array([0, 0])
    my_p = np.cos(theta_p)
    mz_p = np.sin(theta_p)
elif testAxis == 2:
    mx_p = np.sin(theta_p)
    my_p = np.array([0, 0])
    mz_p = np.cos(theta_p)
else:
    print("testAxis must be 0, 1 or 2 for the x, y or z axis respectively")
m0_pv = np.vstack([mx_p, my_p, mz_p]).T     # Initial magnetisation directions

#%% Magtense computations

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
n_macro[testAxis] = 1000 # Number of copies on each side along x, y and z
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

problem.u_ea = np.zeros([2,3])
if testAxis in {0, 1}:
    problem.u_ea=np.array([[0, 0, 1], [0, 0, 1]])   # Easy-axis anisotropy along z
elif testAxis == 2:
    problem.u_ea=np.array([[1, 0, 0], [1, 0, 0]])   # Easy-axis anisotropy along x

# Run simulation
print('Run simulation')
result = problem.run_simulation(
    t_end=t_end,
    nt=nTimesteps,
    fct_h_ext=fct_h_ext,
    nt_h_ext=100,     # Number of datapoints used to fit the zero-function (error when set to 1)
)
print('Done running simulation')

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

#%% Plot results

# Get important data
t_out, M_out = result[:2]
t_n = t_out
M_npv = M_out[:, :, 0, :]

# Select data to plot
tPlot_n = t_n[::dataperiod] * 1e9   # Switch from s to ns
M1_nv = M_npv[:, 0, :]
M2_nv = M_npv[:, 1, :]
theta1_n = np.arctan2(M1_nv[:, 2], M1_nv[:, 0])[::dataperiod]
theta2_n = np.arctan2(-M2_nv[:, 2], M2_nv[:, 0])[::dataperiod]

# Make plot
fig, ax = plt.subplots(layout='constrained', figsize=(6, 4))
ax.plot(tPlot_n, theta1_n, color='forestgreen')
ax.set_xlabel(r'$\text{Time } [\mathrm{ns}]$')
ax.set_ylabel(r'$\text{Polar angle}$')
ax.set_yticks([0, np.pi/4, np.pi/2])
ax.set_yticklabels(['0', r'$\pi/4$', r'$\pi/2$'])
if show:
    fig.show()
