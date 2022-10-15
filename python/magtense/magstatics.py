import numpy as np

from typing import Optional, List, Tuple, Union
from pkg_resources import resource_filename

from magtense.lib import magtensesource
from magtense.utils.geo import euler_to_rot_axis, get_rotmat


class Tiles:
    '''
    Input to Fortran derived type MagTile

    Args:
        center_pos: r0, theta0, z0
        dev_center: dr, dtheta, dz
        size: a, b, c
        vertices: v1, v2, v3, v4 as column vectors.
        M: Mx, My, Mz
        u_ea: Easy axis.
        u_oa: Other axis.
        mu_r_ea: Relative permeability in easy axis.
        mu_r_oa: Relative permeability in other axis.
        M_rem: Remanent magnetization.
        tile_type: 1 = cylinder, 2 = prism, 3 = circ_piece, 4 = circ_piece_inv,
                   5 = tetrahedron, 6 = sphere, 7 = spheroid, 10 = ellipsoid
        offset: Offset of global coordinates.
        rot: Rotation in local coordinate system.
        color: Color in visualization.
        magnet_type: 1 = hard magnet, 2 = soft magnet, 3 = soft + constant mu_r
        stfcn_index: default index into the state function.
        incl_it: If equal to zero the tile is not included in the iteration.
        use_sym: Whether to exploit symmetry.
        sym_op: 1 for symmetry and -1 for anti-symmetry respectively to the planes.
        M_rel: Change in magnetization during last iteration of iterate_magnetization().
        grid_pos: Position in grid if a 3-D grid is used to place tiles.
        n: Number of tiles in the simulation.
    '''
    def __init__(self, 
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
        grid_pos: Optional[List] = None,
        mag_angle: Optional[List] = None,
    ) -> None:
        self._center_pos = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self._dev_center = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self._size = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self._vertices = np.zeros(shape=(n,3,4), dtype=np.float64, order='F')
        self._M = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self._u_ea = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self._u_oa1 = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self._u_oa2 = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self._mu_r_ea = np.ones(shape=(n), dtype=np.float64, order='F')
        self._mu_r_oa = np.ones(shape=(n), dtype=np.float64, order='F')
        self._M_rem = np.zeros(shape=(n), dtype=np.float64, order='F')
        self._tile_type = np.ones(shape=(n), dtype=np.int32, order='F')
        self._offset = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self._rot = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self._color = np.zeros(shape=(n,3), dtype=np.float64, order='F')
        self._magnet_type = np.ones(shape=(n), dtype=np.int32, order='F')
        self._stfcn_index = np.ones(shape=(n), dtype=np.int32, order='F')
        self._incl_it = np.ones(shape=(n), dtype=np.int32, order='F')
        self._use_sym = np.zeros(shape=(n), dtype=np.int32, order='F')
        self._sym_op = np.ones(shape=(n,3), dtype=np.float64, order='F')
        self._M_rel = np.zeros(shape=(n), dtype=np.float64, order='F')
        
        self._grid_pos = np.zeros(shape=(n,3))
        self._n = n

        if center_pos is not None: self.center_pos = center_pos
        if dev_center is not None: self.dev_center = dev_center
        if size is not None: self.size = size
        if vertices is not None: self.vertices = vertices
        if tile_type is not None: self.tile_type = tile_type
        if offset is not None: self.offset = offset
        if rot is not None: self.rot = rot
        if M_rem is not None: self.M_rem = M_rem
        if easy_axis is not None: self.u_ea = easy_axis
        if color is not None: self.color = color
        if magnet_type is not None: self.magnet_type = magnet_type
        if grid_pos is not None: self.grid_pos = grid_pos
        if mag_angle is not None: self.set_mag_angle(mag_angle)

    def __str__(self):
        res = f'{self.n} tiles are present in the setup:\n'
        for i in range(self.n):
            res += f'Tile {i} at grid position ({self.grid_pos[i][0]},' \
                + f'{self.grid_pos[i][1]},{self.grid_pos[i][2]}) with coordinates' \
                + f' x={self.offset[i][0]}, y={self.offset[i][1]}, z={self.offset[i][2]}.\n'
        return res

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, val):
        self._n = val

    @property
    def grid_pos(self):
        return self._grid_pos

    @grid_pos.setter
    def grid_pos(self, val):
        if isinstance(val, Tuple):
            self._grid_pos[val[1]] = np.asarray(val[0])
        else:
            assert len(val) == self.n
            self._grid_pos = np.asarray(val)

    @property
    def center_pos(self):
        return self._center_pos

    @center_pos.setter
    def center_pos(self, val):
        if isinstance(val, Tuple):
            self._center_pos[val[1]] = np.asarray(val[0])
        else:
            if isinstance(val[0], (int, float)):
                self._center_pos[:] = np.asarray(val)
            elif len(val) == self.n:
                self._center_pos = np.asarray(val)

    @property
    def dev_center(self):
        return self._dev_center

    @dev_center.setter
    def dev_center(self, val):
        if isinstance(val, Tuple):
            self._dev_center[val[1]] = np.asarray(val[0])
        else:
            if isinstance(val[0], (int, float)):
                self._dev_center[:] = np.asarray(val)
            elif len(val) == self.n:
                self._dev_center = np.asarray(val)

    @property
    def size(self):
        return self._size

    @size.setter
    def size(self, val):
        if isinstance(val, Tuple):
            self._size[val[1]] = np.asarray(val[0])
        else:
            if isinstance(val[0], (int, float)):
                self._size[:] = np.asarray(val)
            elif len(val) == self.n:
                self._size = np.asarray(val)

    @property
    def vertices(self):
        return self._vertices

    @vertices.setter
    def vertices(self, val):
        if isinstance(val, Tuple):
            vert = np.asarray(val[0])
            if vert.shape == (3,4):
                self._vertices[val[1]] = vert
            elif vert.shape == (4,3):
                self._vertices[val[1]] = vert.T
            else:
                raise ValueError("Four 3-D vertices have to be defined!")
        else:
            vert = np.asarray(val)
            if vert.shape == (3,4):
                self._vertices[:] = vert
            elif vert.shape == (4,3):
                self._vertices[:] = vert.T
            else:
                assert vert[0].shape == (3,4)
                self._vertices = vert

    @property
    def tile_type(self):
        return self._tile_type

    @tile_type.setter
    def tile_type(self, val):
        if isinstance(val, Tuple):
            self._tile_type[val[1]] = val[0]
        else:
            if isinstance(val, (int, float)):
                self._tile_type = np.asarray([val for _ in range(self.n)])
            elif len(val) == self.n:
                self._tile_type = np.asarray(val)

    @property
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, val):
        if isinstance(val, Tuple):
            self._offset[val[1]] = np.asarray(val[0])
        else:
            if isinstance(val[0], (int, float)):
                self._offset[:] = np.asarray(val)
            elif len(val) == self.n:
                self._offset = np.asarray(val)

    @property
    def rot(self):
        return self._rot

    @rot.setter
    def rot(self, val):
        if isinstance(val, Tuple):
            self._rot[val[1]] = euler_to_rot_axis(val[0]) \
                if self.tile_type[val[1]] == 7 else np.asarray(val[0])
        else:
            if isinstance(val[0], (List, np.ndarray)):
                for i in range(self.n):
                    self._rot[i] = euler_to_rot_axis(val[i]) \
                        if self.tile_type[i] == 7 else np.asarray(val[i])
            else:
                self._rot[0] = euler_to_rot_axis(val) \
                    if self.tile_type[0] == 7 else np.asarray(val)

    @property
    def M(self):
        return self._M

    @M.setter
    def M(self, val):
        self._M = val

    @property
    def u_ea(self):
        return self._u_ea

    @property
    def u_oa1(self):
        return self._u_oa1

    @property
    def u_oa2(self):
        return self._u_oa2

    @u_ea.setter
    def u_ea(self, val):
        if isinstance(val, Tuple):
            self._u_ea[val[1]] = np.around(val[0] / np.linalg.norm(val[0]), decimals=9)
            self._M[val[1]] = self.M_rem[val[1]] * self.u_ea[val[1]]
        else:
            if isinstance(val[0], (int, float)): val = [val for _ in range(self.n)]
            for i, ea in enumerate(val):
                self._u_ea[i] = np.around(ea / np.linalg.norm(ea), decimals=9)
                self._M[i] = self.M_rem[i] * self.u_ea[i]
                oa_1 = np.array([val[i][1], -val[i][0], 0])
                self._u_oa1[i] = np.around(oa_1 / np.linalg.norm(oa_1), decimals=9)
                self._u_oa2[i] = np.around(np.cross(self.u_ea[i], self.u_oa1[i]), decimals=9)

    @property
    def mu_r_ea(self):
        return self._mu_r_ea

    @mu_r_ea.setter
    def mu_r_ea(self, val):
        if isinstance(val, Tuple):
            self._mu_r_ea[val[1]] = val[0]
        else:
            if isinstance(val, (int, float)):
                self._mu_r_ea = np.asarray([val for _ in range(self.n)])
            elif len(val) == self.n:
                self._mu_r_ea = np.asarray(val)

    @property
    def mu_r_oa(self):
        return self._mu_r_oa

    @mu_r_oa.setter
    def mu_r_oa(self, val):
        if isinstance(val, Tuple):
            self._mu_r_oa[val[1]] = val[0]
        else:
            if isinstance(val, (int, float)):
                self._mu_r_oa = np.asarray([val for _ in range(self.n)])
            elif len(val) == self.n:
                self._mu_r_oa = np.asarray(val)

    @property
    def M_rem(self):
        return self._M_rem

    @M_rem.setter
    def M_rem(self, val):
        if isinstance(val, Tuple):
            self._M_rem[val[1]] = val[0]
        else:
            if isinstance(val, (int, float)):
                self._M_rem = np.asarray([val for _ in range(self.n)])
            elif len(val) == self.n:
                self._M_rem = np.asarray(val)

    @property
    def color(self):
        return self._color

    @color.setter
    def color(self, val):
        if isinstance(val, Tuple):
            self._color[val[1]] = np.asarray(val[0])
        else:
            if isinstance(val[0], (int, float)):
                self._color[:] = np.asarray(val)
            elif len(val) == self.n:
                self._color = np.asarray(val)

    @property
    def magnet_type(self):
        return self._magnet_type

    @magnet_type.setter
    def magnet_type(self, val):
        if isinstance(val, Tuple):
            self._magnet_type[val[1]] = val[0]
        else:
            if isinstance(val, (int, float)):
                self._magnet_type = np.asarray([val for _ in range(self.n)])
            elif len(val) == self.n:
                self._magnet_type = np.asarray(val)

    @property
    def stfcn_index(self):
        return self._stfcn_index

    @stfcn_index.setter
    def stfcn_index(self, val):
        if isinstance(val, Tuple):
            self._stfcn_index[val[1]] = val[0]
        else:
            if isinstance(val, (int, float)):
                self._stfcn_index = np.asarray([val for _ in range(self.n)])
            elif len(val) == self.n:
                self._stfcn_index = np.asarray(val)

    @property
    def incl_it(self):
        return self._incl_it

    @incl_it.setter
    def incl_it(self, val):
        if isinstance(val, Tuple):
            self._incl_it[val[1]] = val[0]
        else:
            if isinstance(val, (int, float)):
                self._incl_it = np.asarray([val for _ in range(self.n)])
            elif len(val) == self.n:
                self._incl_it = np.asarray(val)

    @property
    def M_rel(self):
        return self._M_rel

    @M_rel.setter
    def M_rel(self, val):
        if isinstance(val, Tuple):
            self._M_rel[val[1]] = val[0]
        else:
            if isinstance(val, (int, float)):
                self._M_rel = [val for _ in range(self.n)]
            elif len(val) == self.n:
                self._M_rel = val

    @property
    def use_sym(self):
        return self._use_sym

    @property
    def sym_op(self):
        return self._sym_op

    def set_mag_angle(self, val=None, seed=42):
        if val is None:
            rng = np.random.default_rng(seed)
            rand = rng.random(size=(self.n,2))
        
        # polar angle [0, pi], azimuth [0, 2*pi]
        for i in range(self.n):
            if val is None:
                self.set_mag_angle_i([np.pi, 2*np.pi] * rand[i], i)
            elif isinstance(val[0], (int, float)):
                self.set_mag_angle_i(val, i)
            else:
                self.set_mag_angle_i(val[i], i)

    def set_mag_angle_i(self, val, i):
        if isinstance(val, (int, float)):
            raise ValueError('Both spherical angles have to be set!')
        else:
            polar_angle, azimuth = val
            self.u_ea = ([np.sin(polar_angle) * np.cos(azimuth),
                          np.sin(polar_angle) * np.sin(azimuth),
                          np.cos(polar_angle)], i)
            self._u_oa1[i] = [np.sin(polar_angle) * np.sin(azimuth),
                              np.sin(polar_angle) * (-np.cos(azimuth)), 0]
            self._u_oa2[i] = [0.5*np.sin(2*polar_angle) * np.cos(azimuth),
                              0.5*np.sin(2*polar_angle) * np.sin(azimuth),
                              -1*np.sin(polar_angle)**2]
    
    def add_tiles(self, n):
        self._center_pos = np.append(self._center_pos, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._dev_center = np.append(self._dev_center, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._size = np.append(self._size, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._vertices = np.append(self._vertices, np.zeros(shape=(n,3,4), dtype=np.float64, order='F'), axis = 0)
        self._M = np.append(self._M, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._u_ea = np.append(self._u_ea, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._u_oa1 = np.append(self._u_oa1, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._u_oa2 = np.append(self._u_oa2, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._mu_r_ea = np.append(self._mu_r_ea, np.ones(shape=(n), dtype=np.float64, order='F'), axis = 0)
        self._mu_r_oa = np.append(self._mu_r_oa, np.ones(shape=(n), dtype=np.float64, order='F'), axis = 0)
        self._M_rem = np.append(self._M_rem, np.zeros(shape=(n), dtype=np.float64, order='F'), axis = 0)
        self._tile_type = np.append(self._tile_type, np.ones(n, dtype=np.int32, order='F'), axis = 0)
        self._offset = np.append(self._offset, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._rot = np.append(self._rot, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._color = np.append(self._color, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._magnettype = np.append(self.magnet_type, np.ones(n, dtype=np.int32, order='F'), axis = 0)
        self._stfcn_index = np.append(self.stfcn_index, np.ones(shape=(n), dtype=np.int32, order='F'), axis = 0)
        self._incl_it = np.append(self.incl_it, np.ones(shape=(n), dtype=np.int32, order='F'), axis = 0)
        self._use_sym = np.append(self.use_sym, np.zeros(shape=(n), dtype=np.int32, order='F'), axis = 0)
        self._sym_op = np.append(self.sym_op, np.ones(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._M_rel = np.append(self.M_rel, np.zeros(shape=(n), dtype=np.float64, order='F'), axis = 0)
        self._grid_pos = np.append(self._grid_pos, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
        self._n += n

    def refinement_prism(self, idx, mat):
        if isinstance(idx, (int, float)):
            self.refinement_prism_i(idx, mat)
        else:
            for i in idx:
                self.refinement_prism_i(i, mat)

    def refinement_prism_i(self, i, mat):
        old_n = self.n
        self.add_tiles(mat[0] * mat[1] * mat[2] - 1)
        R = get_rotmat(self.rot[i])

        x_off, y_off, z_off = self.size[i] / mat
        center_0  = -self.size[i] / 2 + ( self.size[i] / mat ) / 2
        center_pts = [center_0 + [x_off * j, y_off * k, z_off * l] \
            for j in range(mat[0]) for k in range(mat[1]) for l in range(mat[2])]
        ver_cube = (np.dot(R, np.array(center_pts).T)).T + self.offset[i]
        
        for j in range(mat[0]):
            for k in range(mat[1]):
                for l in range(mat[2]):
                    cube_idx = j * mat[1] * mat[2] + k * mat[2] + l
                    idx = cube_idx + old_n
                    if idx == self.n:
                        self._offset[i] = ver_cube[cube_idx]
                        self._size[i] = self.size[i] / mat
                    else:
                        self._size[idx] = self.size[i] / mat
                        self._M[idx] = self.M[i]
                        self._M_rel[idx] = self.M_rel[i]
                        self._color[idx] = self.color[i]
                        self._magnettype[idx] = self.magnet_type[i]
                        self._mu_r_ea[idx] = self.mu_r_ea[i]
                        self._mu_r_oa[idx] = self.mu_r_oa[i]
                        self._rot[idx] = self.rot[i]
                        self._tile_type[idx] = self.tile_type[i]
                        self._u_ea[idx] = self.u_ea[i]
                        self._u_oa1[idx] = self.u_oa1[i]
                        self._u_oa2[idx] = self.u_oa2[i]
                        self._offset[idx] = ver_cube[cube_idx]
                        self._incl_it[idx] = self.incl_it[i]
    

def grid_config(
    spots: Union[List, np.ndarray],
    area: Union[List, np.ndarray],
    filled_pos: Optional[List] = None,
    n_pts: List = [20, 20, 1],
    mode: str = "uniform",
    n_tiles: Optional[int] = None,
    mag_angles: Optional[List] = None,
    B_rem: float = 1.2,
    seed: int = 42
) -> tuple[Tiles, np.ndarray]:
    assert len(spots) == 3
    assert len(area) == 3
    assert len(n_pts) == 3
    
    # Evaluation points in grid structure
    area = np.asarray(area)
    if mode == "uniform":
        seg = area / np.asarray(n_pts)
        pts = [[(i + 0.5) * seg[0], (j + 0.5) * seg[1], (k + 0.5) * seg[2]]
               for i in range(0, n_pts[0])
               for j in range(0, n_pts[1])
               for k in range(0, n_pts[2])]

    elif mode == "center":
        seg = (area[0] / 10) / n_pts[0]
        seg_angle = np.pi / n_pts[1]
        seg_layer = area[2] / n_pts[2]
        pts = [[i * seg, j * seg_angle, k * seg_layer]
               for i in range(0, n_pts[0])
               for j in range(0, n_pts[1])
               for k in range(0, n_pts[2])]
    
    else:
        raise NotImplementedError()

    eval_pts = np.asarray(pts, dtype=np.float64, order='F')

    # Tiles and Grid
    spots = np.asarray(spots)
    tile_size = area / spots
    rng = np.random.default_rng(seed)
    
    # Fill grid randomly
    if filled_pos is None:
        filled_pos = []
        if n_tiles is None:
            n_tiles = rng.integers(np.prod(spots))
        
        cnt_tiles = 0
        while cnt_tiles < n_tiles:
            new_pos = [rng.integers(spots[0]),
                       rng.integers(spots[1]),
                       rng.integers(spots[2])]
            
            if new_pos not in filled_pos:
                filled_pos.append(new_pos)
                cnt_tiles += 1
    
    if len(filled_pos) == 0:
        tiles = None        
    else:
        tiles = Tiles(
            n=len(filled_pos),
            size=tile_size,
            tile_type=2,
            M_rem=B_rem/(4*np.pi*1e-7),
            color=[1, 0, 0],
            grid_pos=filled_pos,
        )

        for i, pos in enumerate(filled_pos):
            if np.greater_equal(np.asarray(pos), spots).any():
                raise ValueError(f'Desired position {pos} is not in the grid!')
            tiles.offset = (np.around((pos * tile_size) + tile_size/2, decimals=9), i) 

        tiles.set_mag_angle(mag_angles)

    return tiles, eval_pts


def run_simulation(
    tiles: Tiles,
    pts: np.array,
    max_error: float = 1e-5,
    max_it: int = 500,
    T: float = 300.,
    console: bool = True
) -> tuple[Tiles, np.ndarray]:
    '''
    Function for running MagTense with the Fortran source code as Python module

    Args:
        tiles: Magnetic tiles to producde magnetic field.
        pts: Evaluation points.
        max_error: Iteration stops if magnetization change below this value.
        max_it: Maximum number of performed iterations.
        T: Temperature for the state function if required.
        console: Boolean if output in console.
    '''
    DATA_PATH = resource_filename('magtense', 'utils/data/data_stateFcn.csv')
    data_stateFcn = np.genfromtxt(DATA_PATH, delimiter=';', dtype=np.float64)

    H_out, M_out, Mrel_out = magtensesource.fortrantopythonio.runsimulation( 
        centerpos=tiles.center_pos,
        dev_center=tiles.dev_center,
        tile_size=tiles.size,
        vertices=tiles.vertices,
        mag=tiles.M,
        u_ea=tiles.u_ea,
        u_oa1=tiles.u_oa1,
        u_oa2=tiles.u_oa2,
        mu_r_ea=tiles.mu_r_ea,
        mu_r_oa=tiles.mu_r_oa,
        mrem=tiles.M_rem,
        tiletype=tiles.tile_type,
        offset=tiles.offset,
        rotangles=tiles.rot,
        color=tiles.color,
        magnettype=tiles.magnet_type,
        statefunctionindex=tiles.stfcn_index,
        includeiniteration=tiles.incl_it,
        exploitsymmetry=tiles.use_sym,
        symmetryops=tiles.sym_op,
        mrel=tiles.M_rel,
        pts=pts,
        data_statefcn=data_stateFcn,
        n_statefcn=1,
        t=T,
        maxerr=max_error,
        nitemax=max_it,
        iteratesolution=True,
        returnsolution=True,
        console=console
    )

    updated_tiles = tiles
    updated_tiles.M = M_out
    updated_tiles.M_rel = Mrel_out
    
    return updated_tiles, H_out


def iterate_magnetization(tiles, max_error=1e-5, max_it=500, T=300., mu_r=20):
    '''
    Iterates through the given tiles to determine their influence on each other.  
    Updated tiles are returned.

    Args:
        max_error: Iteration stops if magnetization change is below this value.
        max_it: Maximum number of performed iterations.
        T: Temperature for the state function if required.
        mu_r: Relative permeability of soft tiles to select the non-linear M-H-curve.
    '''    
    DATA_PATH = resource_filename('magtense', f'utils/data/Fe_mur_{mu_r}_Ms_2_1.csv')
    data_stateFcn = np.genfromtxt(DATA_PATH, delimiter=';', dtype=np.float64)

    M_out, Mrel_out = magtensesource.fortrantopythonio.iteratetiles(
        centerpos=tiles.center_pos,
        dev_center=tiles.dev_center,
        tile_size=tiles.size,
        vertices=tiles.vertices,
        mag=tiles.M,
        u_ea=tiles.u_ea,
        u_oa1=tiles.u_oa1,
        u_oa2=tiles.u_oa2,
        mu_r_ea=tiles.mu_r_ea,
        mu_r_oa=tiles.mu_r_oa,
        mrem=tiles.M_rem,
        tiletype=tiles.tile_type,
        offset=tiles.offset,
        rotangles=tiles.rot,
        color=tiles.color,
        magnettype=tiles.magnet_type,
        statefunctionindex=tiles.stfcn_index,
        includeiniteration=tiles.incl_it,
        exploitsymmetry=tiles.use_sym,
        symmetryops=tiles.sym_op,
        mrel=tiles.M_rel,
        data_statefcn=data_stateFcn,
        n_statefcn=1,
        t=T,
        maxerr=max_error,
        nitemax=max_it
    )

    updated_tiles = tiles
    updated_tiles.M = M_out
    updated_tiles.M_rel = Mrel_out

    return updated_tiles


def get_demag_tensor(tiles, pts):
    '''
    Returns the demagnetization tensor N of the given tiles and the specified evaluation points.
    '''
    demag_tensor = magtensesource.fortrantopythonio.getnfromtiles(
        centerpos=tiles.center_pos,
        dev_center=tiles.dev_center,
        tile_size=tiles.size,
        vertices=tiles.vertices,
        mag=tiles.M,
        u_ea=tiles.u_ea,
        u_oa1=tiles.u_oa1,
        u_oa2=tiles.u_oa2,
        mu_r_ea=tiles.mu_r_ea,
        mu_r_oa=tiles.mu_r_oa,
        mrem=tiles.M_rem,
        tiletype=tiles.tile_type,
        offset=tiles.offset,
        rotangles=tiles.rot,
        color=tiles.color,
        magnettype=tiles.magnet_type,
        statefunctionindex=tiles.stfcn_index,
        includeiniteration=tiles.incl_it,
        exploitsymmetry=tiles.use_sym,
        symmetryops=tiles.sym_op,
        mrel=tiles.M_rel,
        pts=pts
    )

    return demag_tensor


def get_H_field(tiles, pts, demag_tensor=None):
    '''
    Returns the H field at the specified evaluation points of the given tiles.  
    Optionally, a precalculated demagnetization tensor N can be handed over.
    This prevents unnecessary and expensive recalculation of the demag tensor
    if the geometry of the setup does not change.
    '''
    useN = False if demag_tensor is None else True
    if demag_tensor is None:
        demag_tensor = np.zeros(shape=(tiles.n,len(pts),3,3), dtype=np.float64, order='F')
        
    H_out = magtensesource.fortrantopythonio.gethfromtiles(
        centerpos=tiles.center_pos,
        dev_center=tiles.dev_center,
        tile_size=tiles.size,
        vertices=tiles.vertices,
        mag=tiles.M, u_ea=tiles.u_ea,
        u_oa1=tiles.u_oa1,
        u_oa2=tiles.u_oa2,
        mu_r_ea=tiles.mu_r_ea,
        mu_r_oa=tiles.mu_r_oa,
        mrem=tiles.M_rem,
        tiletype=tiles.tile_type,
        offset=tiles.offset,
        rotangles=tiles.rot,
        color=tiles.color,
        magnettype=tiles.magnet_type,
        statefunctionindex=tiles.stfcn_index,
        includeiniteration=tiles.incl_it,
        exploitsymmetry=tiles.use_sym,
        symmetryops=tiles.sym_op,
        mrel=tiles.M_rel,
        pts=pts,
        n=demag_tensor,
        usestoredn=useN
    )                 

    return H_out
