"""Test micromagnetics solver with macrogeometry, finite temperature and the shape correction field"""

# General modules
import numpy as np
import matplotlib.pyplot as plt

# Specialised modules
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

#%% Settings

# Constants of nature
mu0 = 4*np.pi * 1e-7    # Vacuum permeability [N/A^2]
gamma = 2.21*1e5        # Gyromagnetic ratio [rad*m / (A*s)]

# Random number seed
np.random.seed(42)

# Parameters to sweep over
Nparams = 4    # Number of distinct parameter sets
tStep_s = [0.5*1e-12]*Nparams
Bamp_s = [1]*Nparams
Bfreq_s = [1e8]*Nparams
wGap_s = [0, 1e-15]*2
# Ncell_s = [9]*Nparams   # Micromagnetic cells
Ncell_s = [6]*Nparams     # Micromagnetic cells
nRep_s = [2]*Nparams      # Repetitions of the simulated domain: (2n+1) X (2n+1) X (2n+1) grid of domain copies
alpha_s = [1]*Nparams
a, b, c = 7*1e-9, 7*1e-9, 7*1e-9    # Single-cell sidelengths
A_s, B_s, C_s = [a * Ncell_s[s] for s in range(Nparams)], [b * Ncell_s[s] for s in range(Nparams)], [c * Ncell_s[s] for s in range(Nparams)]  # [m]
include_shape_correction_s = [True]*Nparams
# include_shape_correction_s = [False]*Nparams
aSample_s, bSample_s, cSample_s = [1]*Nparams, [1]*Nparams, [3]*Nparams
randomAniDir_s = [False]*Nparams    # Whether to randomise the direction of the anisotropy axes
u_v = np.array([1, 0, 0])   # Direction of anisotropy axis when randomAniDir_s == False
Ngrain_s = [1]*Nparams   # Number of grains along x, y and z with different anisotropy axes. Only relevant if randomAniDir_s == True
T_s = [0]*2 + [300]*2  # Temperature [K]
s0 = 0  # Which parameter set to start with
# sMax = 5
sMax = Nparams

# Where to save
save = True    # Whether to save
base_savename = 'example_calc'

# Plot settings
show = True                      # Whether to show the plots
dataperiodPlot = 2               # Sampling rate for plotted datapoints

# Timestepping
Ncycles = 1.3               # Number of field cycles to simulate
nTimestepsField = 2000      # Number of time coordinates to evaluate field at for the purpose of interpolation

# How much data to save
Ndata = 1000 * Ncycles    # Number of datapoints to save

# Micromagnetic solver settings
cuda=False      # Set to False because I only implemented the shape correction for this case
cvode=False     # When cvode=True, loss of moment normalisation is a computational bottleneck

print(f'Run the simulation sequence : {base_savename}')
# print(f'Of Nparams = {Nparams} parameter combinations, run from index s0 = {s0} to sMax = {sMax}')
print(f'Nparams = {Nparams}, s0 = {s0}, sMax = {sMax}')
for s in tqdm(range(s0, sMax)):
    # Material properties
    Ms = 1.5               # Saturation magnetisation in magnetic region [T]
    Ms = Ms/mu0            # [A/m]
    Aex = 1e-11            # Exchange constant in magnetic region [J/m]
    alpha = alpha_s[s]     # Gilbert damping [-]
    eta = alpha/(1 + alpha**2) * gamma  # Damping constant [m/(A*s)]
    K = 8*1e3              # Uniaxial anisotropy in magnetic region [J/m^3]
    T = T_s[s]             # Temperature [K]

    # Applied field
    Bamp = Bamp_s[s]    # Field amplitude [T]
    Bfreq = Bfreq_s[s]  # Field frequency [Hz]
    def fct_h_ext(t) -> np.ndarray:
        e_z = np.array([[0, 0, 1]]).T   # Unit vector along z as 3-by-1 column vector
        # NOTE : sign change to counter weird sign change in MagTense
        amplitude = np.sin(2 * np.pi * Bfreq * t) * Bamp/mu0    # sine-wave for amplitude
        return np.transpose(amplitude * e_z)    # Has to end with the shape (nt, 3)

    # Timestepping
    t_end = 1/Bfreq * Ncycles           # Total time simulated [s]
    t_step = tStep_s[s]                 # Timestep [s]
    nTimesteps = int(t_end / t_step)    # Number of timesteps

    # Initial magnetisation
    direction = 'random'        # Directions of initial magnetisation. 'random' => uniformly random polar angles
    stdDir = np.pi/3     # Standard deviation for the gaussian distributed shifts in polar angle (2*stdDir for azimuthal)

    # Geometry
    A, B, C = A_s[s], B_s[s], C_s[s]  # Side lengths of simulation domain
    Ncell = Ncell_s[s]   # Number of micromagnetic cells along each axis
    NcellX, NcellY, NcellZ = Ncell, Ncell, Ncell  # Number of micromagnetic cells along x, y and z
    wGap = wGap_s[s]                 # Thickness of gap layer
    Ngrain = Ngrain_s[s]      # Number of grains along each axis
    NgrainX, NgrainY, NgrainZ = Ngrain, Ngrain, Ngrain

    # Macrogeometry
    nRep = nRep_s[s]
    Nrep = (2*nRep + 1)
    NrepX, NrepY, NrepZ = Nrep, Nrep, Nrep   # Repetitions in each direction
    nRepX, nRepY, nRepZ = nRep, nRep, nRep   # Repetitions in each direction on both sides

    # Sample shape (sidelengths of equivalent rectangular prism)
    include_shape_correction = include_shape_correction_s[s]
    # The absolute values are unimportant as long as they greatly exceed the macrogeometry size,
    # but relative values are important in determining sample shape.
    aSample, bSample, cSample = aSample_s[s], bSample_s[s], cSample_s[s]

    # Where to save
    savename = f'{base_savename}_{s}'


    #%% Micromagnetic computations

    # Side lengths of individual micromagnetic cells
    a, b, c = A / NcellX, B / NcellY, C / NcellZ

    # Number of cells
    NcellTot = NcellX * NcellY * NcellZ

    # Positions
    pos_pv = np.zeros([NcellTot, 3])
    for i in range(NcellX):
        for j in range(NcellY):
            for k in range(NcellZ):
                p = k + NcellZ*j + NcellY*i
                pos_pv[p, :] = np.array([i*a, j*b, k*c])
    pos_pv = pos_pv + np.array([[a/2, b/2, c/2]]) - np.array([[A/2, B/2, C/2]])

    # Initial magnetisations
    theta0, phi0 = np.random.uniform(0, np.pi), np.random.uniform(0, 2*np.pi)
    theta_p = theta0 + np.random.normal(loc = 0, scale = stdDir, size=[NcellTot])
    phi_p = phi0 + np.random.normal(loc = 0, scale = 2*stdDir, size=[NcellTot])
    Mx_p = Ms * np.sin(theta_p) * np.cos(phi_p)
    My_p = Ms * np.sin(theta_p) * np.sin(phi_p)
    Mz_p = Ms * np.cos(theta_p)
    M_pv = np.vstack([Mx_p, My_p, Mz_p]).T
    # Reformat arrays in tiles object to be consistent with M
    m0_pv = M_pv / Ms    # Initial magnetisation directions

    # Information about macrogeometry and sample shape
    Xshift, Yshift, Zshift = A + wGap, B + wGap, C + wGap   # How far to translate domain copies in each direction
    n_macro = np.array([nRepX, nRepY, nRepZ])
    shiftVec = np.array([Xshift, Yshift, Zshift])
    if include_shape_correction:
        macroShape = np.array([Xshift * NrepX, Yshift * NrepY, Zshift * NrepZ])
        sampleShape = np.array([aSample, bSample, cSample])
    else:
        macroShape, sampleShape = None, None

    # Uniform anisotropy
    u_pv = np.tile(u_v, [NcellTot, 1])

    # Boundary conditions on exchange
    if wGap == 0:
        exchPBC = np.array([True, True, True])
    else:
        exchPBC = np.array([False, False, False])

    #%% Magtense computations

    # Define micromagnetic problem
    # grid_type = "unstructuredPrisms"
    grid_type = "uniform"
    grid_L = [A, B, C]              # Grid size along x, y and z
    res = [NcellX, NcellY, NcellZ]  # Grid resolution, i.e. number of cells along x, y and z

    problem = MicromagProblem(
        res=res,
        A0=Aex,
        Ms=Ms,
        K0=K,
        alpha=eta,
        m0=m0_pv,
        cuda=cuda,
        cvode=cvode,
        n_macro=n_macro,
        shiftVec=shiftVec,
        macroShape=macroShape,
        sampleShape=sampleShape,
        grid_L=grid_L,
        grid_type=grid_type,
        usereturnhall=True,
        solver='dynamic',
        T=T,
        exchPBC=exchPBC,
    )

    problem.u_ea=u_pv

    # Run simulation
    print('Run simulation')
    result = problem.run_simulation(
        t_end=t_end,
        nt=nTimesteps,
        fct_h_ext=fct_h_ext,
        nt_h_ext=nTimestepsField,
    )
    print('Done running simulation')

    t_out, M_out = result[:2]
    pts = result[2]
    H_exc, H_ext, H_dem, H_ani = result[3:7]


    #%% Plot results

    # Get important data
    t_n = t_out
    M_npv = M_out[:, :, 0, :]
    Mavg_nv = np.mean(M_npv, axis=1)
    Bext_nv = H_ext[:, 0, 0, :] * mu0

    tPlot_n = t_n[::dataperiodPlot] * 1e9
    Mplot_nv = Mavg_nv[::dataperiodPlot]
    Bplot_nv = Bext_nv[::dataperiodPlot]

    fig, ax = plt.subplots(layout='constrained', figsize=(6, 4))
    ax.plot(Bplot_nv[:, 2], Mplot_nv[:, 2], color='forestgreen')
    ax.set_xlabel(r'$B_z$ [T]')
    ax.set_ylabel(r'$m_z$')
    if show:
        fig.show()
    if save:
        fig.savefig(f'figures/{savename}_Hysteresis.pdf', format='pdf')

    fig, (ax1, ax2) = plt.subplots(layout='constrained', ncols=2, figsize=(8, 4))
    ax1.plot(tPlot_n, Bplot_nv[:, 2], color='crimson')
    ax2.plot(tPlot_n, Mplot_nv[:, 2], color='navy')
    ax1.set_xlabel(r'$t$ [ns]')
    ax1.set_ylabel(r'$B_z$ [T]')
    ax2.set_xlabel(r'$t$ [ns]')
    ax2.set_ylabel(r'$m_z$')
    if show:
        fig.show()
    if save:
        fig.savefig(f'figures/{savename}_Time_evolution.pdf', format='pdf')