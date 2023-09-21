from typing import List, Optional, Tuple, Union

import numpy as np
from pkg_resources import resource_filename

from magtense.lib import magtensesource


class Tile:
    """
    This is a proxy class to enable Tiles to be accessed as if they were a list of Tile
    objects. So, for example,

    t = Tiles(5)
    print(t[0].center_pos)  # [0, 0, 0]
    t[0].center_pos = [1, 1, 1]
    print(t[0].center_pos)  # [1, 1, 1]
    """

    def __init__(self, parent, index):
        self._parent = parent
        self._index = index

    def __getattr__(self, name):
        # Get the attribute from the parent
        attr = getattr(self._parent, name)
        # If the attribute is a np.ndarray, return the element at our index
        if isinstance(attr, np.ndarray):
            return attr[self._index]
        # Otherwise, just return the attribute itself
        return attr

    def __setattr__(self, name, value):
        if name in ["_parent", "_index"]:
            # For these attributes, set them normally
            super().__setattr__(name, value)
        else:
            # For all other attributes, set the value in the parent's attribute array
            attr = getattr(self._parent, name)
            if isinstance(attr, np.ndarray):
                attr[self._index] = value


class Tiles:
    """
    Input to Fortran derived type MagTile.

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
        n: Number of tiles in the simulation.
    """

    def __init__(
        self,
        n: int,
        center_pos: Optional[List] = None,
        dev_center: Optional[List] = None,
        size: Optional[List] = None,
        vertices: Optional[List] = None,
        mu_r_ea: Union[int, List, None] = None,
        mu_r_oa: Union[int, List, None] = None,
        tile_type: Union[int, List, None] = None,
        offset: Optional[List] = None,
        rot: Optional[List] = None,
        M_rem: Union[int, List, None] = None,
        easy_axis: Optional[List] = None,
        color: Optional[List] = None,
        magnet_type: Optional[List] = None,
        mag_angle: Optional[List] = None,
    ) -> None:
        self._center_pos = np.zeros(shape=(n, 3), dtype=np.float64, order="F")
        self._dev_center = np.zeros(shape=(n, 3), dtype=np.float64, order="F")
        self._size = np.zeros(shape=(n, 3), dtype=np.float64, order="F")
        self._vertices = np.zeros(shape=(n, 3, 4), dtype=np.float64, order="F")
        self._M = np.zeros(shape=(n, 3), dtype=np.float64, order="F")
        self._u_ea = np.zeros(shape=(n, 3), dtype=np.float64, order="F")
        self._u_oa1 = np.zeros(shape=(n, 3), dtype=np.float64, order="F")
        self._u_oa2 = np.zeros(shape=(n, 3), dtype=np.float64, order="F")
        self._mu_r_ea = np.ones(shape=(n), dtype=np.float64, order="F")
        self._mu_r_oa = np.ones(shape=(n), dtype=np.float64, order="F")
        self._M_rem = np.zeros(shape=(n), dtype=np.float64, order="F")
        self._tile_type = np.ones(shape=(n), dtype=np.int32, order="F")
        self._offset = np.zeros(shape=(n, 3), dtype=np.float64, order="F")
        self._rot = np.zeros(shape=(n, 3), dtype=np.float64, order="F")
        self._color = np.zeros(shape=(n, 3), dtype=np.float64, order="F")
        self._magnet_type = np.ones(shape=(n), dtype=np.int32, order="F")
        self._stfcn_index = np.ones(shape=(n), dtype=np.int32, order="F")
        self._incl_it = np.ones(shape=(n), dtype=np.int32, order="F")
        self._use_sym = np.zeros(shape=(n), dtype=np.int32, order="F")
        self._sym_op = np.ones(shape=(n, 3), dtype=np.float64, order="F")
        self._M_rel = np.zeros(shape=(n), dtype=np.float64, order="F")
        self._n = n

        if center_pos is not None:
            self.center_pos = center_pos
        if dev_center is not None:
            self.dev_center = dev_center
        if size is not None:
            self.size = size
        if vertices is not None:
            self.vertices = vertices
        if mu_r_ea is not None:
            self.mu_r_ea = mu_r_ea
        if mu_r_oa is not None:
            self.mu_r_oa = mu_r_oa
        if tile_type is not None:
            self.tile_type = tile_type
        if offset is not None:
            self.offset = offset
        if rot is not None:
            self.rot = rot
        if M_rem is not None:
            self.M_rem = M_rem
        if easy_axis is not None:
            self.u_ea = easy_axis
        else:
            self.set_easy_axis(mag_angle)
        if color is not None:
            self.color = color
        if magnet_type is not None:
            self.magnet_type = magnet_type

    def __str__(self):
        res = ""
        for i in range(self.n):
            res += f"Tile_{i} with coordinates {self.offset[i]}.\n"
        return res

    def __getitem__(self, index):
        if index < 0 or index >= self.n:
            raise IndexError("Index out of bounds")
        return Tile(self, index)

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, val):
        self._n = val

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
            if vert.shape == (3, 4):
                self._vertices[val[1]] = vert
            elif vert.shape == (4, 3):
                self._vertices[val[1]] = vert.T
            else:
                raise ValueError("Four 3-D vertices have to be defined!")
        else:
            vert = np.asarray(val)
            if vert.shape == (3, 4):
                self._vertices[:] = vert
            elif vert.shape == (4, 3):
                self._vertices[:] = vert.T
            else:
                assert vert[0].shape == (3, 4)
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
                self._tile_type = np.asarray(
                    [val for _ in range(self.n)], dtype=np.int32
                )
            elif len(val) == self.n:
                self._tile_type = np.asarray(val, dtype=np.int32)

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
            self._rot[val[1]] = (
                _euler_to_rot_axis(val[0])
                if self.tile_type[val[1]] == 7
                else np.asarray(val[0])
            )
        else:
            if isinstance(val[0], (List, np.ndarray)):
                for i in range(self.n):
                    self._rot[i] = (
                        _euler_to_rot_axis(val[i])
                        if self.tile_type[i] == 7
                        else np.asarray(val[i])
                    )
            else:
                self._rot[0] = (
                    _euler_to_rot_axis(val)
                    if self.tile_type[0] == 7
                    else np.asarray(val)
                )

    @property
    def M(self):
        return self._M

    @M.setter
    def M(self, val):
        if isinstance(val, Tuple):
            self._M[val[1]] = val[0]
        else:
            if isinstance(val[0], (int, float)):
                assert len(val) == 3
                self._M = np.asarray([val for _ in range(self.n)])
            elif len(val) == self.n:
                assert len(val[0]) == 3
                self._M = np.asarray(val)

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
            self.M = (self.M_rem[val[1]] * self.u_ea[val[1]], val[1])
        else:
            if isinstance(val[0], (int, float)):
                val = [val for _ in range(self.n)]
            for i, ea in enumerate(val):
                self._u_ea[i] = np.around(ea / np.linalg.norm(ea), decimals=9)
                self.M = (self.M_rem[i] * self.u_ea[i], i)
                if ea[1] != 0 or ea[2] != 0:
                    w = np.array([1, 0, 0])
                else:
                    w = np.array([0, 1, 0])
                self._u_oa1[i] = np.around(np.cross(self.u_ea[i], w), decimals=9)
                self._u_oa2[i] = np.around(
                    np.cross(self.u_ea[i], self.u_oa1[i]), decimals=9
                )

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

    def set_easy_axis(
        self, val: Optional[List] = None, idx: Optional[int] = None, seed: int = 42
    ) -> None:
        """
        polar angle [0, pi], azimuth [0, 2*pi]
        """
        if val is None:
            rng = np.random.default_rng(seed)
            rand = rng.random(size=(self.n, 2))

        if idx is None:
            for i in range(self.n):
                if val is None:
                    self._set_ea_i([np.pi, 2 * np.pi] * rand[i], i)
                elif isinstance(val[0], (int, float)):
                    self._set_ea_i(val, i)
                else:
                    self._set_ea_i(val[i], i)
        else:
            i_val = [np.pi, 2 * np.pi] * rand[idx] if val is None else val
            self._set_ea_i(i_val, idx)

    def _set_ea_i(self, val, i):
        if isinstance(val, (int, float)):
            raise ValueError("Both spherical angles have to be set!")
        else:
            polar_angle, azimuth = val
            self.u_ea = (
                [
                    np.sin(polar_angle) * np.cos(azimuth),
                    np.sin(polar_angle) * np.sin(azimuth),
                    np.cos(polar_angle),
                ],
                i,
            )
            self._u_oa1[i] = [
                np.sin(polar_angle) * np.sin(azimuth),
                np.sin(polar_angle) * (-np.cos(azimuth)),
                0,
            ]
            self._u_oa2[i] = [
                0.5 * np.sin(2 * polar_angle) * np.cos(azimuth),
                0.5 * np.sin(2 * polar_angle) * np.sin(azimuth),
                -1 * np.sin(polar_angle) ** 2,
            ]

    def _add_tiles(self, n: int) -> None:
        self._center_pos = np.append(
            self._center_pos,
            np.zeros(shape=(n, 3), dtype=np.float64, order="F"),
            axis=0,
        )
        self._dev_center = np.append(
            self._dev_center,
            np.zeros(shape=(n, 3), dtype=np.float64, order="F"),
            axis=0,
        )
        self._size = np.append(
            self._size, np.zeros(shape=(n, 3), dtype=np.float64, order="F"), axis=0
        )
        self._vertices = np.append(
            self._vertices,
            np.zeros(shape=(n, 3, 4), dtype=np.float64, order="F"),
            axis=0,
        )
        self._M = np.append(
            self._M, np.zeros(shape=(n, 3), dtype=np.float64, order="F"), axis=0
        )
        self._u_ea = np.append(
            self._u_ea, np.zeros(shape=(n, 3), dtype=np.float64, order="F"), axis=0
        )
        self._u_oa1 = np.append(
            self._u_oa1, np.zeros(shape=(n, 3), dtype=np.float64, order="F"), axis=0
        )
        self._u_oa2 = np.append(
            self._u_oa2, np.zeros(shape=(n, 3), dtype=np.float64, order="F"), axis=0
        )
        self._mu_r_ea = np.append(
            self._mu_r_ea, np.ones(shape=(n), dtype=np.float64, order="F"), axis=0
        )
        self._mu_r_oa = np.append(
            self._mu_r_oa, np.ones(shape=(n), dtype=np.float64, order="F"), axis=0
        )
        self._M_rem = np.append(
            self._M_rem, np.zeros(shape=(n), dtype=np.float64, order="F"), axis=0
        )
        self._tile_type = np.append(
            self._tile_type, np.ones(n, dtype=np.int32, order="F"), axis=0
        )
        self._offset = np.append(
            self._offset, np.zeros(shape=(n, 3), dtype=np.float64, order="F"), axis=0
        )
        self._rot = np.append(
            self._rot, np.zeros(shape=(n, 3), dtype=np.float64, order="F"), axis=0
        )
        self._color = np.append(
            self._color, np.zeros(shape=(n, 3), dtype=np.float64, order="F"), axis=0
        )
        self._magnet_type = np.append(
            self._magnet_type, np.ones(n, dtype=np.int32, order="F"), axis=0
        )
        self._stfcn_index = np.append(
            self._stfcn_index, np.ones(shape=(n), dtype=np.int32, order="F"), axis=0
        )
        self._incl_it = np.append(
            self._incl_it, np.ones(shape=(n), dtype=np.int32, order="F"), axis=0
        )
        self._use_sym = np.append(
            self._use_sym, np.zeros(shape=(n), dtype=np.int32, order="F"), axis=0
        )
        self._sym_op = np.append(
            self._sym_op, np.ones(shape=(n, 3), dtype=np.float64, order="F"), axis=0
        )
        self._M_rel = np.append(
            self._M_rel, np.zeros(shape=(n), dtype=np.float64, order="F"), axis=0
        )
        self._n += n

    def refine_prism(self, idx: Union[int, List], mat: List):
        if isinstance(idx, (int, float)):
            self._refine_prism_i(idx, mat)
        else:
            for i in idx:
                self._refine_prism_i(i, mat)

    def _refine_prism_i(self, i, mat):
        assert self.tile_type[i] == 2
        old_n = self.n
        self._add_tiles(np.prod(mat) - 1)
        R_mat = get_rotmat(self.rot[i])

        x_off, y_off, z_off = self.size[i] / mat
        c_0 = -self.size[i] / 2 + (self.size[i] / mat) / 2
        c_pts = [
            c_0 + [x_off * j, y_off * k, z_off * m]
            for j in range(mat[0])
            for k in range(mat[1])
            for m in range(mat[2])
        ]
        ver_cube = (np.dot(R_mat, np.array(c_pts).T)).T + self.offset[i]

        for j in range(mat[0]):
            for k in range(mat[1]):
                for m in range(mat[2]):
                    cube_idx = j * mat[1] * mat[2] + k * mat[2] + m
                    idx = cube_idx + old_n
                    if idx == self.n:
                        self._offset[i] = ver_cube[cube_idx]
                        self._size[i] = self.size[i] / mat
                    else:
                        self._size[idx] = self.size[i] / mat
                        self._M[idx] = self.M[i]
                        self._M_rel[idx] = self.M_rel[i]
                        self._color[idx] = self.color[i]
                        self._magnet_type[idx] = self.magnet_type[i]
                        self._mu_r_ea[idx] = self.mu_r_ea[i]
                        self._mu_r_oa[idx] = self.mu_r_oa[i]
                        self._rot[idx] = self.rot[i]
                        self._tile_type[idx] = self.tile_type[i]
                        self._u_ea[idx] = self.u_ea[i]
                        self._u_oa1[idx] = self.u_oa1[i]
                        self._u_oa2[idx] = self.u_oa2[i]
                        self._offset[idx] = ver_cube[cube_idx]
                        self._incl_it[idx] = self.incl_it[i]


def _euler_to_rot_axis(euler):
    """
    Converting Euler angles to rotation axis.
    For spheroids, the rotation axis has to be set rather than Euler angles.
    Rotation in MagTense performed in local coordinate system:
    Euler - (1) Rot_X_L, (2) Rot_Y_L', (3) Rot_Z_L''
    """
    # symm-axis of geometry points in the direction of z-axis of L
    # Rotates given rotation axis with pi/2 around y_L''
    # Moves x-axis of L'' to z-axis
    # ax[0] = ax[2], ax[1] = ax[1], ax[2] = -ax[0]
    R_y = np.array(
        [
            [np.cos(np.pi / 2), 0, np.sin(np.pi / 2)],
            [0, 1, 0],
            [-np.sin(np.pi / 2), 0, np.cos(np.pi / 2)],
        ]
    )
    ax = np.dot(R_y, np.asarray(euler).T)

    # Calculate the spherical coordinates: yaw and pitch
    # x_L'' has to match ax
    # Perform negative yaw around x_L and pitch around y_L'
    # The azimuthal angle is offset by pi/2 (zero position of x_L'')
    rot_x = -np.arctan2(ax[1], ax[0])
    rot_y = np.arccos(ax[2] / np.sqrt(ax[0] ** 2 + ax[1] ** 2 + ax[2] ** 2)) - np.pi / 2

    return np.array([rot_x, rot_y, 0])


def get_rotmat(rot):
    """
    G to L in local coordinate system
    TODO: Check rotation from local to global: (1) Rot_X, (2) Rot_Y, (3) Rot_Z
    """
    rot_x = (
        [1, 0, 0],
        [0, np.cos(rot[0]), -np.sin(rot[0])],
        [0, np.sin(rot[0]), np.cos(rot[0])],
    )
    rot_y = (
        [np.cos(rot[1]), 0, np.sin(rot[1])],
        [0, 1, 0],
        [-np.sin(rot[1]), 0, np.cos(rot[1])],
    )
    rot_z = (
        [np.cos(rot[2]), -np.sin(rot[2]), 0],
        [np.sin(rot[2]), np.cos(rot[2]), 0],
        [0, 0, 1],
    )
    return np.asarray(rot_x) @ np.asarray(rot_y) @ np.asarray(rot_z)


def grid_config(
    spots: Union[List, np.ndarray],
    area: Union[List, np.ndarray],
    filled_pos: Optional[List] = None,
    n_pts: List = [20, 20, 1],
    mode: str = "uniform",
    n_tiles: Optional[int] = None,
    mag_angles: Optional[List] = None,
    B_rem: float = 1.2,
    seed: int = 42,
) -> tuple[Tiles, np.ndarray]:
    spots = np.asarray(spots)
    area = np.asarray(area)
    rng = np.random.default_rng(seed)

    if mode == "uniform":
        seg = area / np.asarray(n_pts)
        pts = [
            [(i + 0.5) * seg[0], (j + 0.5) * seg[1], (k + 0.5) * seg[2]]
            for i in range(0, n_pts[0])
            for j in range(0, n_pts[1])
            for k in range(0, n_pts[2])
        ]

    elif mode == "center":
        seg = (area[0] / 10) / n_pts[0]
        seg_angle = np.pi / n_pts[1]
        seg_layer = area[2] / n_pts[2]
        pts = [
            [i * seg, j * seg_angle, k * seg_layer]
            for i in range(0, n_pts[0])
            for j in range(0, n_pts[1])
            for k in range(0, n_pts[2])
        ]

    else:
        raise NotImplementedError()

    if filled_pos is None:
        if n_tiles is None:
            n_tiles = 1 + rng.integers(np.prod(spots))
        elif n_tiles < 1 or n_tiles > np.prod(spots):
            raise ValueError("n_tiles out of range!")

        # Generate unique linear indices
        linear_indices = rng.choice(np.prod(spots), size=n_tiles, replace=False)

        # Convert the linear indices to 3D coordinates
        filled_pos = np.empty((n_tiles, 3), dtype=int)
        for idx, linear_index in enumerate(linear_indices):
            filled_pos[idx, 2] = linear_index // (spots[0] * spots[1])
            filled_pos[idx, 1] = (linear_index % (spots[0] * spots[1])) // spots[0]
            filled_pos[idx, 0] = linear_index % spots[0]
        filled_pos = filled_pos.tolist()
    elif len(filled_pos) == 0:
        raise ValueError("filled_pos is empty!")

    tiles = Tiles(
        n=len(filled_pos),
        size=area / spots,
        tile_type=2,
        M_rem=B_rem / (4 * np.pi * 1e-7),
        color=[1, 0, 0],
        mag_angle=mag_angles,
        offset=[(np.asarray(pos) + 0.5) * (area / spots) for pos in filled_pos],
    )

    eval_pts = np.asarray(pts, dtype=np.float64, order="F")

    return tiles, eval_pts


def run_simulation(
    tiles: Tiles,
    pts: np.ndarray,
    max_error: float = 1e-5,
    max_it: int = 500,
    T: float = 300.0,
    mu_r: float = 20,
    console: bool = True,
) -> tuple[Tiles, np.ndarray]:
    """
    Run magnetostatic simulation to calculate the demagnetizing field strength.

    Args:
        tiles: Magnetic tiles to produce magnetic field.
        pts: Evaluation points.
        max_error: Iteration stops if magnetization change below this value.
        max_it: Maximum number of performed iterations.
        T: Temperature for the state function if required.
        console: Boolean if output in console.

    Returns:
        Updated tiles.
        Demagnetizing field strength in evaluation points.
    """
    DATA_PATH = resource_filename("magtense", f"mat/Fe_mur_{mu_r}_Ms_2_1.csv")
    data_stateFcn = np.genfromtxt(DATA_PATH, delimiter=";", dtype=np.float64)

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
        console=console,
    )

    tiles.M = M_out
    tiles.M_rel = Mrel_out

    return tiles, H_out


def iterate_magnetization(
    tiles: Tiles,
    max_error: float = 1e-5,
    max_it: int = 500,
    T: float = 300.0,
    mu_r: float = 20,
) -> Tiles:
    """
    Iterate through tiles to determine their influence on each other.

    Args:
        tiles: Magnetic tiles to be iterated over.
        max_error: Iteration stops if magnetization change is below this value.
        max_it: Maximum number of performed iterations.
        T: Temperature for the state function if required.
        mu_r: Relative permeability of soft tiles to select the non-linear M-H-curve.

    Returns:
        Updated tiles.
    """
    DATA_PATH = resource_filename("magtense", f"mat/Fe_mur_{mu_r}_Ms_2_1.csv")
    data_stateFcn = np.genfromtxt(DATA_PATH, delimiter=";", dtype=np.float64)

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
        nitemax=max_it,
    )

    tiles.M = M_out
    tiles.M_rel = Mrel_out

    return tiles


def get_demag_tensor(tiles: Tiles, pts: np.ndarray) -> np.ndarray:
    """
    Get demagnetization tensor of tiles and the specified evaluation points.

    Args:
        tiles: Magnetic tiles to produce magnetic field.
        pts: Evaluation points.

    Returns:
        Demagnetization tensor.
    """
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
        pts=pts,
    )

    return demag_tensor


def get_H_field(
    tiles: Tiles, pts: np.ndarray, demag_tensor: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Calculate the demagnetizing field strength of a magnetic setup.

    Args:
        tiles: Magnetic tiles to produce magnetic field.
        pts: Evaluation points.
        demag_tensor: Optional precalculated demagnetization tensor.
            This prevents unnecessary and expensive recalculation
            if the geometry of the setup does not change.

    Returns:
        Demagnetizing field strength in evaluation points.
    """
    if demag_tensor is None:
        useN = False
        demag_tensor = np.zeros(
            shape=(tiles.n, len(pts), 3, 3), dtype=np.float64, order="F"
        )
    else:
        useN = True

    H_out = magtensesource.fortrantopythonio.gethfromtiles(
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
        n=demag_tensor,
        usestoredn=useN,
    )

    return H_out
