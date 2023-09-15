Documentation
========================================

========================================
Python
========================================

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Tiles [`source <https://github.com/cmt-dtu-energy/MagTense/blob/00179ccaa29a5c452de1aa1f6991df2bdc9ed9e1/python/src/magtense/magstatics.py#L9>`_]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    self.init(
        n: int,
        center_pos: Optional[List] = None,
        dev_center: Optional[List] = None,
        size: Optional[List] = None,
        vertices: Optional[List] = None,
        tile_type: Union[int, List, None] = None,
        offset: Optional[List] = None,
        rot: Optional[List] = None,
        M_rem: Union[int, List, None] = None,
        easy_axis: Optional[List] = None,
        color: Optional[List] = None,
        magnet_type: Optional[List] = None,
        mag_angle: Optional[List] = None,
    ) -> None

    self.set_easy_axis(
        val: Optional[List] = None,
        idx: Optional[int] = None,
        seed: int = 42
    ) -> None

    self.refine_prism(
        idx: Union[int, List],
        mat: List
    ) -> None


----------------------------------------
Functions
----------------------------------------

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
grid_config [`source <https://github.com/cmt-dtu-energy/MagTense/blob/00179ccaa29a5c452de1aa1f6991df2bdc9ed9e1/python/src/magtense/magstatics.py#L528>`_]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    grid_config(
        spots: Union[List, np.ndarray],
        area: Union[List, np.ndarray],
        filled_pos: Optional[List] = None,
        n_pts: List = [20, 20, 1],
        mode: str = "uniform",
        n_tiles: Optional[int] = None,
        mag_angles: Optional[List] = None,
        B_rem: float = 1.2,
        seed: int = 42
    ) -> tuple[Tiles, np.ndarray]


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
iterate_magnetization [`source <https://github.com/cmt-dtu-energy/MagTense/blob/00179ccaa29a5c452de1aa1f6991df2bdc9ed9e1/python/src/magtense/magstatics.py#L666>`_]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Iterate through tiles to determine their influence on each other.

::

    iterate_magnetization(
        tiles: Tiles,
        max_error: float = 1e-5,
        max_it: int = 500,
        T: float = 300.,
        mu_r: float = 20
    ) -> Tiles

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
get_demag_tensor [`source <https://github.com/cmt-dtu-energy/MagTense/blob/00179ccaa29a5c452de1aa1f6991df2bdc9ed9e1/python/src/magtense/magstatics.py#L724>`_]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Get demagnetization tensor of tiles and the specified evaluation points.

::

    get_demag_tensor(
        tiles: Tiles,
        pts: np.ndarray
    ) -> np.ndarray

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
get_H_field [`source <https://github.com/cmt-dtu-energy/MagTense/blob/00179ccaa29a5c452de1aa1f6991df2bdc9ed9e1/python/src/magtense/magstatics.py#L763>`_]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Calculate the demagnetizing field strength of a magnetic setup.

::
    
    get_H_field(
        tiles: Tiles,
        pts: np.ndarray,
        demag_tensor: Optional[np.ndarray] = None
    ) -> np.ndarray

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
run_simulation [`source <https://github.com/cmt-dtu-energy/MagTense/blob/00179ccaa29a5c452de1aa1f6991df2bdc9ed9e1/python/src/magtense/magstatics.py#L601>`_]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Run MagTense with the Fortran source code as Python module.

::

    run_simulation(
        tiles: Tiles,
        pts: np.ndarray,
        max_error: float = 1e-5,
        max_it: int = 500,
        T: float = 300.,
        console: bool = True
    ) -> tuple[Tiles, np.ndarray]:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
create_plot [`source <https://github.com/cmt-dtu-energy/MagTense/blob/00179ccaa29a5c452de1aa1f6991df2bdc9ed9e1/python/src/magtense/utils.py#L452>`_]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Creates a plot with the iterated tiles and the calculated magnetic field H at the
evaluation points as quiver plot. Additionally, an optional grid can be displayed.
**Tile types**: 1 = cylinder, 2 = prism, 3 = circ_piece, 4 = circ_piece_inv, 5 = tetrahedron, 6 = sphere, 7 = spheroid, 10 = ellipsoid

::

    create_plot(
        tiles: Optional[Tiles] = None,
        eval_pts: Optional[np.ndarray] = None,
        field: Optional[np.ndarray] = None,
        spots: Optional[List] = None,
        area: Optional[List] = None,
    ) -> None
