"""muMag standard problem 3: the flower/vortex energy crossover of a cube.

The cube is relaxed in the flower and in the vortex state at a sequence of edge lengths, in units
of the exchange length. The energy of the two states crosses at the single domain limit. This is
the Python counterpart of
matlab/examples/Micromagnetism/mumag_micromag_Std_problem_3/Standard_problem_3.m and takes the
same options.

The problem can be run on three kinds of mesh, selected with mesh_type:

  'uniform'            a regular grid of res[0] x res[1] x res[2] cells. This is the default and
                       the only mesh type that uses the res argument.
  'unstructuredPrisms' a grid of unstructured Cartesian prisms read from the text file mesh_file
                       (one row per prism: centre [x,y,z] and side lengths [a,b,c] in metres).
                       The mesh is scaled to each edge length in L_loop.
  'tetrahedron'        a tetrahedral mesh generated at each edge length with create_tetra_mesh,
                       with mesh_res_param cells per edge length. The nodes and the connectivity
                       are handed to MagTense, which analyses the mesh and builds the exchange
                       operator itself.

For the two unstructured meshes the mesh is the only thing given to MagTense: it computes the
cell volumes, the exchange operator and the energies itself. The energies of all mesh types come
back in problem.E_out, evaluated in Fortran with the cell volumes of the mesh.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from magtense.micromag import MicromagProblem
from magtense.utils import create_tetra_mesh, tetra_volumes

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
MESH_DIR = (
    REPOSITORY_ROOT / "documentation" / "examples_mumag_validation" / "Validation_standard_problem_3"
)
DEFAULT_MESH_FILE = MESH_DIR / "Std_prob_3_unstructured_cartesian_grains_9_mesh_4_ref_2.txt"
MESH_TYPES = ("uniform", "unstructuredPrisms", "tetrahedron")


def std_prob_3(
    res: tuple[int, int, int] = (10, 10, 10),
    L_loop: np.ndarray | None = None,
    mesh_type: str = "uniform",
    mesh_file: Path | str = DEFAULT_MESH_FILE,
    mesh_res_param: float = 5,
    cuda: bool = True,
    cvode: bool = False,
    use_minimizer: bool = True,
    plotting: bool = True,
    plot_details: bool = False,
    figpath: Path | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Run the muMag standard problem 3 for a range of cube sizes.

    Args:
        res: cells along x, y and z of the uniform grid. Ignored for the other mesh types.
        L_loop: the cube edge lengths to simulate, in units of the exchange length. The default
            is the ten values of the Matlab example, linspace(8, 9, 10).
        mesh_type: 'uniform', 'unstructuredPrisms' or 'tetrahedron', see the module docstring.
        mesh_file: the text file with the unstructured Cartesian mesh, used with
            mesh_type='unstructuredPrisms'.
        mesh_res_param: tetrahedra per edge length, used with mesh_type='tetrahedron'.
        cuda: use CUDA for the calculations.
        cvode: use CVODE for the numerical time evolution.
        use_minimizer: relax each state with the energy minimizer (solver 'explicit') instead of
            integrating the Landau-Lifshitz equation in time at zero field over t_end (the
            'dynamic' solver). The minimizer ignores t_end and stops when the largest torque is
            below problem.min_tol. Either way the number of effective-field evaluations spent is
            printed, which is the cost to compare.
        plotting: show (or save, with figpath) the energy of the two states against L.
        plot_details: also show the magnetization and the energies at each L.
        figpath: directory to save the figures in. None shows them interactively.

    Returns:
        L_loop: the cube edge lengths in units of the exchange length that were simulated.
        E_arr: the equilibrium energies with shape (4, len(L_loop), 2), reduced by Km*V with
            Km = mu0 Ms^2 / 2. The first axis is (demagnetisation, exchange, anisotropy,
            external), the order the Matlab example has always used, and the last axis selects
            the flower state (0) or the vortex state (1). The single domain limit is the L where
            the two total energies cross.
    """
    if mesh_type not in MESH_TYPES:
        raise ValueError(f"mesh_type must be one of {MESH_TYPES}, got {mesh_type!r}")

    mu0 = 4 * np.pi * 1e-7
    if L_loop is None:
        L_loop = np.linspace(8, 9, 10)
    L_loop = np.atleast_1d(np.asarray(L_loop, dtype=np.float64))

    # The material parameters. These are common to all mesh types
    alpha = 1e3
    gamma = 0.0
    Ms = 1000e3
    K0 = 0.1 * 0.5 * mu0 * Ms**2
    A0 = 1.74532925199e-10

    # The exchange length and the energy scale of the mumag problem
    lex = np.sqrt(A0 / (0.5 * mu0 * Ms**2))
    Km = 0.5 * mu0 * Ms**2

    def h_ext_fct(t) -> np.ndarray:
        return np.atleast_2d(t).T * np.array([0, 0, 0])

    # The unstructured Cartesian mesh is stored for a single edge length and is scaled to each
    # value in L_loop. The edge length it was made at is recovered from the total volume of the
    # prisms, which is the cube
    if mesh_type == "unstructuredPrisms":
        data = np.loadtxt(mesh_file, comments="%")
        pos_mesh = data[:, 0:3]
        dims_mesh = data[:, 3:6]
        L_mesh = np.sum(np.prod(dims_mesh, axis=1)) ** (1 / 3)
        print(f"Loaded {len(pos_mesh)} prisms meshed at L = {L_mesh / lex:.4g} l_ex")

    E_arr = np.zeros(shape=(4, len(L_loop), 2))
    n_feval = np.zeros(shape=(len(L_loop), 2), dtype=int)

    for i, L in enumerate(L_loop):
        print(f"ITERATION: {i + 1} / {len(L_loop)}")
        grid_L = np.array([lex, lex, lex]) * L

        # Setup the mesh. Each branch gives the mesh arguments of the problem, the cell centres,
        # pts, which are used below to set up the vortex state, and the cell volumes, which the
        # mean magnetization has to be weighted by
        if mesh_type == "uniform":
            ntot = int(np.prod(res))
            mesh_kwargs = {}
            x, y, z = np.meshgrid(
                np.linspace(-1, 1, res[0]), np.linspace(-1, 1, res[1]), np.linspace(-1, 1, res[2]),
                indexing="ij",
            )
            # The tiles are ordered with x running fastest, so the arrays are flattened that way
            pts = np.stack([a.swapaxes(0, 2).reshape(-1) for a in (x, y, z)], axis=1)
            volumes = np.full(ntot, np.prod(grid_L) / ntot)

        elif mesh_type == "unstructuredPrisms":
            ntot = len(pos_mesh)
            print(f"Prisms N_grid = {ntot}")
            # Scale the stored mesh to the current edge length
            scale = L / (L_mesh / lex)
            grid_pts = pos_mesh * scale
            grid_abc = dims_mesh * scale
            mesh_kwargs = {"grid_type": mesh_type, "grid_pts": grid_pts, "grid_abc": grid_abc}
            pts = grid_pts
            volumes = np.prod(grid_abc, axis=1)

        else:
            # Create the tetrahedral mesh, with a cell size that scales with the cube. The mesh is
            # all MagTense needs: it runs the mesh analysis and builds the exchange operator
            nodes, elements, pts = create_tetra_mesh(grid_L, lex * L / mesh_res_param)
            ntot = len(pts)
            print(f"Tetra N_grid = {ntot}")
            mesh_kwargs = {
                "grid_type": mesh_type, "grid_nod": nodes, "grid_ele": elements.T, "grid_pts": pts,
            }
            volumes = tetra_volumes(nodes, elements)

        # Setup the problem. The equilibrium is found either by integrating the LL equation in
        # time at zero field over t_end (the 'dynamic' solver) or by the energy minimizer, which
        # ignores the time window and stops when the largest torque is below problem.min_tol.
        # The energies come back in problem.E_out in both cases.
        problem = MicromagProblem(
            res=res if mesh_type == "uniform" else (ntot, 1, 1),
            grid_L=grid_L,
            solver="explicit" if use_minimizer else "dynamic",
            A0=A0,
            Ms=Ms,
            K0=K0,
            alpha=alpha,
            gamma=gamma,
            cuda=cuda,
            cvode=cvode,
            usereturnhall=True,
            **mesh_kwargs,
        )
        problem.u_ea[:, 2] = 1

        for j in range(2):
            if j == 0:
                print("Flower state")
                m0 = np.zeros((ntot, 3))
                m0[:, 2] = 1
                t_end = 10e-9

            else:
                print("Vortex state")
                # A vortex in the xz-plane about the centre of the cube. The centre is taken from
                # the extent of the cell centres, so the mesh can sit anywhere
                center = (pts.min(axis=0) + pts.max(axis=0)) / 2
                angle = np.arctan2(pts[:, 2] - center[2], pts[:, 0] - center[0])
                m0 = np.zeros((ntot, 3))
                m0[:, 0] = np.sin(angle)
                m0[:, 2] = -np.cos(angle)
                t_end = 200e-9

            problem.m0 = m0

            # For the LL relaxation a convergence check at every output time lets the integration
            # stop as soon as the magnetization is stationary instead of running to t_end.
            nt = 50
            problem.t_conv = np.linspace(0, t_end, nt)
            problem.nt_conv = nt
            problem.conv_tol = np.repeat(1e-6, nt)

            # The minimizer treats every row of the field table as a separate field to relax at,
            # so it gets a single row; the dynamic solver interpolates the table in time.
            t, M_out = problem.run_simulation(
                t_end=t_end, nt=nt, fct_h_ext=h_ext_fct, nt_h_ext=1 if use_minimizer else 2
            )[:2]

            # The energy terms come back from Fortran in J as (time, field, term) with the terms
            # in the order exchange, external, demag, anisotropy, evaluated from the fields the
            # solver ran on and the cell volumes of the mesh. Divided by Km*V they are the
            # reduced energies of the mumag problem.
            E_red = problem.E_out[:, 0, :] / (Km * np.prod(grid_L))
            E_exc, E_ext, E_dem, E_ani = E_red.T
            # dem, exc, ani, ext - the order the Matlab example has always used
            E_arr[:, i, j] = np.array([E_dem[-1], E_exc[-1], E_ani[-1], E_ext[-1]])

            n_feval[i, j] = int(problem.n_feval.sum())
            if use_minimizer:
                print(
                    f"   Minimizer: {n_feval[i, j]} field evaluations, "
                    f"{int(problem.min_iter.sum())} iterations, "
                    f"status {int(problem.min_status[-1])}, E/(Km V) = {E_arr[:, i, j].sum():.6f}"
                )
            else:
                print(
                    f"   LL relaxation: {n_feval[i, j]} field evaluations, "
                    f"E/(Km V) = {E_arr[:, i, j].sum():.6f}"
                )

            if plot_details:
                _plot_details(t, M_out[:, :, 0, :], pts, volumes, E_red, i, j, figpath)

    print(f"Total field evaluations: {n_feval.sum()}")

    if plotting:
        plt.figure()
        plt.plot(L_loop, np.sum(E_arr[:, :, 0], axis=0), ".", markersize=12, label="Flower state")
        plt.plot(L_loop, np.sum(E_arr[:, :, 1], axis=0), ".", markersize=12, label="Vortex state")
        plt.xlabel("L [$l_{ex}$]")
        plt.ylabel("E [-]")
        plt.legend()
        plt.grid(True)
        plt.title(f"Standard problem 3, mesh: {mesh_type}")
        if figpath is None:
            plt.show()
        else:
            figpath.mkdir(parents=True, exist_ok=True)
            plt.savefig(figpath / "3_solution.png")

    return L_loop, E_arr


def _plot_details(
    t: np.ndarray,
    M: np.ndarray,
    pts: np.ndarray,
    volumes: np.ndarray,
    E_red: np.ndarray,
    i: int,
    j: int,
    figpath: Path | None,
) -> None:
    """The starting and the ending magnetization, and the mean magnetization and the energies
    against time, for one relaxation. Works for every mesh type, since it only uses the cell
    centres and the cell volumes."""
    state = "flower" if j == 0 else "vortex"
    weights = volumes / volumes.sum()
    m_avg = np.einsum("tnc,n->tc", M, weights)

    fig = plt.figure(figsize=(12, 8))
    for k, (label, m) in enumerate((("Starting magnetization", M[0]), ("Ending magnetization", M[-1]))):
        ax = fig.add_subplot(2, 2, 2 * k + 1, projection="3d")
        ax.quiver(pts[:, 0], pts[:, 1], pts[:, 2], m[:, 0], m[:, 1], m[:, 2], length=0.5 * np.ptp(pts[:, 0]) / 10)
        ax.set_title(label)
        ax.set_box_aspect((np.ptp(pts[:, 0]), np.ptp(pts[:, 1]), np.ptp(pts[:, 2])))

    ax = fig.add_subplot(2, 2, 2)
    for c, (label, color) in enumerate((("$<m_x>$", "r"), ("$<m_y>$", "g"), ("$<m_z>$", "b"))):
        ax.plot(t, m_avg[:, c], color + "d", label=label)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Reduced magnetization, $m_i$ [-]")
    ax.legend()
    ax.grid(True)

    ax = fig.add_subplot(2, 2, 4)
    for c, label in enumerate((r"$E_{exc}$", r"$E_{ext}$", r"$E_{dem}$", r"$E_{ani}$")):
        ax.plot(t, E_red[:, c] - E_red[0, c], ".", label=label)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Energy change [-]")
    ax.legend()
    ax.grid(True)

    fig.suptitle(f"Standard problem 3, L index {i + 1}, {state} state")
    fig.tight_layout()
    if figpath is None:
        plt.show()
    else:
        figpath.mkdir(parents=True, exist_ok=True)
        fig.savefig(figpath / f"3_details_L{i + 1}_{state}.png")
    plt.close(fig)


if __name__ == "__main__":
    std_prob_3(
        mesh_type="uniform",   # 'uniform', 'unstructuredPrisms' or 'tetrahedron'
        cuda=True,
        cvode=False,
        use_minimizer=True,    # set False to relax by integrating the LL equation in time instead
        plotting=True,
        plot_details=False,
        figpath=None,
    )
