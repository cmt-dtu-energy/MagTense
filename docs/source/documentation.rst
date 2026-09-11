Documentation
========================================

This page lists the public Python API. The authoritative source is always the
code itself, which carries full docstrings:

* `magstatics.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/src/magtense/magstatics.py>`_
  - the magnetostatic framework
* `micromag.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/src/magtense/micromag.py>`_
  - the micromagnetic framework
* `utils.py <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/src/magtense/utils.py>`_
  - plotting and helper routines

For the Matlab interface, the corresponding entry points are the MEX-files
listed under :ref:`Matlab` and the helper functions in
`matlab/util <https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/util>`_.

========================================
Magnetostatics
========================================

----------------------------------------
Tiles class
----------------------------------------

``magtense.magstatics.Tiles`` holds a collection of ``n`` magnetic tiles. Every
per-tile property is exposed as an array attribute, and ``tiles[i]`` returns a
view of a single tile.

::

    Tiles(
        n: int,
        center_pos: list[float] | None = None,
        dev_center: list[float] | None = None,
        size: list[float] | None = None,
        vertices: list[float] | None = None,
        M_rem: float | list[float] | None = None,
        easy_axis: list[float] | None = None,
        mag_angle: list[float] | None = None,
        mu_r_ea: float | list[float] | None = None,
        mu_r_oa: float | list[float] | None = None,
        tile_type: int | list[int] | None = None,
        offset: list[float] | None = None,
        rot: list[float] | None = None,
        color: list[float] | None = None,
        magnet_type: list[int] | None = None,
    )

Attributes: ``n``, ``center_pos``, ``dev_center``, ``size``, ``vertices``,
``tile_type``, ``offset``, ``rot``, ``M``, ``u_ea``, ``u_oa1``, ``u_oa2``,
``mu_r_ea``, ``mu_r_oa``, ``M_rem``, ``color``, ``magnet_type``,
``stfcn_index``, ``incl_it``, ``M_rel``, ``use_sym``, ``sym_op``. Their meaning
is described under :ref:`Geometry`, :ref:`Magnetization` and
:ref:`Other parameters`.

Methods:

::

    set_easy_axis(
        val: list | None = None,
        idx: int | None = None,
        seed: int = 42
    ) -> None

Sets the easy axis from a polar angle in :math:`[0,\pi]` and an azimuth in
:math:`[0, 2\pi]`, drawing them at random if ``val`` is omitted.

::

    refine_prism(
        idx: int | float | list,
        mat: list
    ) -> None

Subdivides one or more prism tiles into ``mat = [nx, ny, nz]`` sub-prisms.

----------------------------------------
Magnetostatic functions
----------------------------------------

::

    grid_config(
        spots: list | np.ndarray,
        area: list | np.ndarray,
        filled_pos: list | np.ndarray | None = None,
        n_pts: tuple = (20, 20, 1),
        mode: str = "uniform",
        n_tiles: int | None = None,
        mag_angles: list | None = None,
        B_rem: float = 1.2,
        seed: int = 42
    ) -> tuple[Tiles, np.ndarray]

Builds a grid of prism tiles together with a set of evaluation points.

::

    run_simulation(
        tiles: Tiles,
        pts: np.ndarray,
        obs_size: np.ndarray = None,
        max_error: float = 1e-5,
        max_it: int = 500,
        T: float = 300.0,
        mu_r: float = 20,
    ) -> tuple[Tiles, np.ndarray]

Iterates the magnetization of the tiles and returns them together with the
H-field in the evaluation points. ``obs_size`` is an optional argument which
gives the size of the observation cell, when using the volume-averaged prism 
(tile_type = 8) tensor.

::

    iterate_magnetization(
        tiles: Tiles,
        max_error: float = 1e-5,
        max_it: int = 500,
        T: float = 300.0,
        mu_r: float = 20
    ) -> Tiles

Iterates through the tiles to determine their influence on each other.

::

    get_demag_tensor(
        tiles: Tiles,
        pts: np.ndarray,
        obs_size: np.ndarray = None
    ) -> np.ndarray

Returns the demagnetization tensor of the tiles at the evaluation points.

::

    get_H_field(
        tiles: Tiles,
        pts: np.ndarray,
        demag_tensor: np.ndarray | None = None
    ) -> np.ndarray

Calculates the demagnetizing field of a magnetic setup. Passing a
precalculated ``demag_tensor`` avoids recomputing it when the geometry does not
change.

::

    get_H_field_fmm(
        tiles, pts, eps=1e-6, nterms_in=10, cells_per_node=10,
        nlmin=0, nlmax=2, ifunif=1, do_target=0, do_FI=1, n_pts=-1
    ) -> np.ndarray

FMM-backed H-field. This is a test version that uses a single dipole per tile
and requires a build with ``USE_FMM3D=1``.

::

    get_rotmat(rot: np.ndarray) -> np.ndarray

Rotation matrix corresponding to the three rotation angles of a tile.

========================================
Micromagnetic problem class
========================================

``magtense.micromag.MicromagProblem`` defines and runs a micromagnetic problem.
Every constructor argument is documented in
:ref:`Micromagnetic parameter reference`; the parameters that are not
constructor arguments are plain attributes and are listed there as well.

::

    MicromagProblem(
        res: list[int],
        grid_L: list[float] = (500e-9, 125e-9, 3e-9),
        grid_nnod: int = 0,
        grid_type: str | None = "uniform",
        prob_mode: str | None = "new",
        solver: str | None = "dynamic",
        m0=None, A0=None, Ms=None, K0=None, K1=None, K2=None, T=None,
        K0_arr=None, CrysAxis=None,
        alpha: float = 4.42e3,
        gamma: float = 2.21e5,
        max_T0: float = 2.0,
        nt_conv: int = 1,
        conv_tol: float = 1e-4,
        tol: float = 1e-4,
        thres: float = 1e-6,
        setTimeDis: int = 10,
        dem_thres: float = 0.0,
        demag_approx: str | None = None,
        cv: float = 0.0,
        grid_pts=None, grid_abc=None,
        exch_val=None, exch_rows=None, exch_cols=None,
        exch_nval: int = 1, exch_nrow: int = 1, exch_ncols: int = 1,
        exch_intpn: str | None = "extended",
        exch_meth: str | None = "directlaplacianneumann",
        exch_weigh: float = 8,
        exch_presize: int = 12,
        demigstp: int = 0,
        passexch: int = 0,
        filename: str = "t",
        cuda: bool = False,
        cvode: bool = False,
        usedemag: bool = True,
        useavgn: bool = True,
        usereturnhall: bool = False,
        precision: bool = False,
        n_threads: int = 1,
        N_ave: tuple[int] = (1, 1, 1),
        t_alpha: np.ndarray = np.zeros(1),
        alpha_fct=lambda t: np.atleast_2d(t).T * 0,
        sampleShape=None, exchPBC=None, macroShape=None,
        shiftVec=np.zeros(3), n_macro=np.zeros(3),
        hysteresis_solver: str = "static",
        rng_seed: int = 0,
    ) -> None

----------------------------------------
Micromagnetic run methods
----------------------------------------

::

    run_simulation(
        t_end: float,
        nt: int,
        fct_h_ext: Callable,
        nt_h_ext: int
    ) -> list

A single solve, used for both the ``dynamic`` and the ``explicit`` solver.
``fct_h_ext`` is evaluated on ``nt_h_ext`` uniformly spaced times between 0 and
``t_end`` and must return an ``(nt_h_ext, 3)`` array. The returned list is
described in :ref:`Micromagnetic output`.

::

    run_hysteresis(H_ext: np.ndarray) -> list

A predefined sequence of applied fields, given as an ``(n,4)`` array of
``[t, Hx, Hy, Hz]`` rows. Requires ``hysteresis_solver="static"``.

::

    run_hysteresis_adaptive(
        H_start: np.ndarray,
        H_end: np.ndarray,
        dH_initial: float,
        dH_min: float,
        dH_max: float,
        max_steps: int,
        dM_min: float = 1e-3,
        dM_target: float = 1e-2,
        dM_reject: float = 5e-2,
        dH_grow: float = 1.5,
        dH_shrink: float = 0.75,
        switch_refine_dH: float | None = None,
    ) -> list

A field sweep in which the solver picks the field steps itself. Requires
``hysteresis_solver="adaptive"`` and ``solver="explicit"``. See
:ref:`Adaptive hysteresis`.

========================================
Plotting utilities
========================================

::

    create_plot(
        tiles: Tiles | None = None,
        eval_pts: list | None = None,
        field: np.ndarray | None = None,
        spots: list | None = None,
        area: list | None = None,
        show: bool = True,
    ) -> None

Plots the tiles together with the calculated H-field at the evaluation points
as a quiver plot, and optionally a grid.

**Tile types**: 1 = cylinder, 2 = prism, 3 = circ_piece, 4 = circ_piece_inv,
5 = tetrahedron, 6 = sphere, 7 = spheroid, 8 = avg_prism, 10 = ellipsoid

::

    plot_M_thin_film(
        m: np.ndarray,
        res: list[int],
        title: str | None = None,
        scale: float | None = None,
        width: float = 0.002,
        headwidth: float = 3,
        headlength: float = 5,
        figpath: Path | None = None,
    ) -> None

Plots the in-plane magnetization of a thin-film micromagnetic problem.

::

    plot_magfield(
        field: np.ndarray,
        magnet: np.ndarray | None = None,
        vmax: float = 1
    ) -> None

Plots the components of a field on a regular grid as images.

Further helpers in ``magtense.utils`` include ``plot_M_avg_seq``,
``plot_field``, ``plot_grid``, ``plot_Halbach``, ``load_COMSOL`` and
``validation``, as well as the individual tile plotting routines
``plot_cube``, ``plot_sphere``, ``plot_spheroid``, ``plot_tetrahedron``,
``plot_cylindrical`` and ``plot_circpiece``.
