import numpy as np

from typing import Optional, List

from magtense import magstatics


class EvaluationPoints:
    '''
    Points at which the demagnetizing field strength shall be calculated.
    A 3-D uniform grid in the center of the cylinder bore is generated.

    Args:
        res: Resolution of the grid of evaluation points.
        r: Radial extension in each direction.
        c: Center point of cylinder bore.
        line: Boolean to produce a 1-D line in each direction.
        coords: Optionally, pre-defined coordinates can be used.
    '''
    def __init__(
        self,
        res: List[int] = [16, 16, 16],
        r: float = 0.002,
        c: List[float] = [0,0,0],
        line: bool = False,
        coords: Optional[np.array] = None
    ) -> None:
        assert len(res) == 3
        l_n = [complex(0, n_res) for n_res in res]
        self._res = res
        self._r = r

        if coords is not None:
            self._coords = coords
            return

        if line:
            _coords = np.vstack((np.mgrid[c[0]-r:c[0]+r:l_n[0], 0:0:1j, 0:0:1j],
                                 np.mgrid[0:0:1j, c[1]-r:c[1]+r:l_n[1], 0:0:1j],
                                 np.mgrid[0:0:1j, 0:0:1j, c[2]-r:c[2]+r:l_n[2]]))
        else:
            zl = c[2] if res[2] == 1 else c[2] - r
            zr = c[2] if res[2] == 1 else c[2] + r
            _coords = np.mgrid[c[0]-r:c[0]+r:l_n[0],
                               c[1]-r:c[1]+r:l_n[1],
                               zl:zr:l_n[2]]

        # Transpose from 'ij' to 'xy' notation
        self._coords = _coords.transpose((0,2,1,3)).reshape(3,-1).T

    @property
    def coordinates(self):
        return self._coords

    @property
    def res(self):
        return self._res

    @property
    def r(self):
        return self._r


class HalbachCylinder:
    """An assembly of cylindrical magnetic pieces in the Halbach configuration.
    Setup is adapted from https://www.sciencedirect.com/science/article/pii/S0304885320307514.
    Optimal angle: https://www.sciencedirect.com/science/article/pii/S0304885309008816
    The resulting magnetic field is pointing in x-direction.

    Args:
        n_layers: Number of cylindrical layers in Halbach setup.
        n_segs: Number of segments each Halbach cylinder is built of.
        n_axial: Number of axial components of each segments.
        B_rem [T]: Remanent magnetic flux density, when H_ext=0 A/m.
        r_bore: Radius of the cylinder bore.
        l: Height of each cylindrical tile.
        r_o: Outer radius of the Halbach cylinder.
        mag_res: Each cylindrical tile is simulated with {mag_res} magnets.
        M_rem: Pre-defined remanent magnetization of each tile.
        easy_axes: Pre-defined easy axis of each tile.
    """
    def __init__(
        self,
        n_layers: int = 1,
        n_segs: int = 16,
        n_axial: int = 1,
        B_rem: float = 1.22,
        r_bore: float = 0.02,
        l: float = 0.2,
        r_o: float = 0.2,
        mag_res: int = 1,
        M_rem: Optional[np.array] = None,
        easy_axes: Optional[np.array] = None,
    ) -> None:
        self.n_hard_tiles = n_layers * mag_res * n_segs * mag_res * n_axial
        self.r_bore = r_bore
        self.n_layers = n_layers
        self.n_segs = n_segs
        self.n_axial = n_axial
        self.mag_res = mag_res
        self.B_rem = B_rem
        self.l_azimuth_opt = []

        # Attributes for shim magnets
        self._pts_shim = None
        self.shim_layers = 0
        self.shim_segs = 0
        self.n_shim = 0
        self.r_i = None
        self.mu_r = 100

        eps = 10e-12

        tiles = magstatics.Tiles(
            n=self.n_hard_tiles,
            tile_type=1,
            mu_r_ea=1.06,
            mu_r_oa=1.17,
            dev_center=[(r_o - r_bore - eps) / (n_axial * mag_res),
                        2 * np.pi / (n_segs * mag_res) - eps * np.pi, l - eps],
            M_rem=M_rem,
            easy_axis=easy_axes
        )

        for z in range(n_layers):
            z0 = - n_layers / 2 * l + z * l + l / 2

            for seg in range(n_segs):
                if easy_axes is None:
                    theta_seg = (2 * np.pi / n_segs) * seg + np.pi / n_segs + eps * np.pi
                    azimuth_opt = 2 * np.arctan2(np.sin(theta_seg), np.cos(theta_seg))

                for res_seg in range(mag_res):
                    theta0 = (2 * np.pi / (n_segs * mag_res)) * (seg * mag_res + res_seg) \
                           + np.pi / (n_segs * mag_res) + eps * np.pi

                    for axial in range(n_axial):
                        for res_ax in range(mag_res):
                            idx = z * n_segs * mag_res * n_axial * mag_res \
                                + seg * mag_res * mag_res * n_axial \
                                + res_seg * mag_res * n_axial + axial * mag_res + res_ax

                            r0 = ((r_o - r_bore) / (n_axial * mag_res)) * (axial * mag_res + res_ax + 0.5) + r_bore
                            tiles.center_pos = ([r0, theta0, z0], idx)
                            if M_rem is None: tiles.M_rem = (B_rem / (4 * np.pi * 1e-7), idx)
                            if easy_axes is None: tiles.set_easy_axis([np.pi/2, azimuth_opt], idx)
                            self.l_azimuth_opt.append(azimuth_opt)
                            color = [0, 0, 0]
                            color[axial - 3 * (axial // 3)] = 1
                            tiles.color = (color, idx)

        self._struct = tiles

    @property
    def tiles(self) -> magstatics.Tiles:
        return self._struct

    @property
    def remanence(self) -> np.array:
        return self.tiles.M_rem[:self.n_hard_tiles]

    @property
    def easy_axes(self) -> np.array:
        return self.tiles.u_ea[:self.n_hard_tiles]

    @property
    def shimming_points(self) -> Optional[EvaluationPoints]:
        return self._pts_shim


    def perturb_config(self, B_rem, azimuth=None, polar=None):
        '''
        Pertubation of the optimal configuration for the finite, segmented Halbach cylinder.
        Norm, azimuth, or polar angle are drawn from a normal distribution.
        This should reflect manufacturing limitations.
        '''
        assert len(B_rem) == self.n_hard_tiles
        self.reset()
        self.tiles.M_rem = B_rem / (4 * np.pi * 1e-7)
        if azimuth is None: azimuth = self.l_azimuth_opt
        if polar is None: polar = np.ones(self.n_hard_tiles) * np.pi / 2
        mag_angles = np.concatenate((np.expand_dims(polar, axis=1), np.expand_dims(azimuth, axis=1)), axis=1)
        
        for idx in range(self.n_hard_tiles):
            self.tiles.set_easy_axis(mag_angles[idx], idx)


    def add_shim_magnets(self, n_layers, n_segs=None, shim_grid=None, mu_r=None, r_i=0.01):
        """
        Shim magnets are placed inside the bore of the Halbach cylinder.
        Realization as soft magnets with constant relative permeability (soft_type=2).
        """
        self.shim_layers = n_layers
        self.r_i = r_i
        if mu_r is not None: self.mu_r = mu_r
        soft_type = 2
        l = 0.1 / n_layers
        self.shim_grid = shim_grid

        shim_dict = {8:[12, 12], 12:[13, 13], 16:[12,14], 20:[14, 14],
                     32:[15, 15], 52:[16, 16], 176:[24, 24], 1148:[50, 50]}

        if shim_grid is None:
            if n_segs is None:
                print('No shim magnet configuration specified!')
                return
            else:
                self.shim_segs = n_segs
                self.shim_grid = shim_dict[n_segs]
        
        else:
            self.shim_segs = list(shim_dict.keys())[list(shim_dict.values()).index(shim_grid)]

        pts_shim = np.zeros([self.shim_segs * n_layers, 3])

        idx = self.n_hard_tiles
        shim_idx = 0
        for z in range(n_layers):
            z0 = - n_layers / 2 * l + z * l + l / 2
            size = (2 * self.r_bore) / np.array(self.shim_grid)

            for grid_x in range(self.shim_grid[0]):
                for grid_y in range(self.shim_grid[1]):
                    off = [(grid_x + 0.5) * size[0] - self.r_bore,
                              (grid_y + 0.5) * size[1] - self.r_bore, z0]
                    dist_c = off[0]**2 + off[1]**2
                    
                    if dist_c < (self.r_bore - np.sqrt(size[0]**2 + size[1]**2))**2 \
                        and dist_c > (r_i + np.sqrt(size[0]**2 + size[1]**2))**2:
                        pts_shim[shim_idx] = np.array(off)
                        self.tiles._add_tiles(1)
                        self.tiles.offset = (off, idx)
                        self.tiles.tile_type = (2, idx)
                        self.tiles.M_rem = (0, idx)
                        self.tiles.size = ([size[0] / 2, size[1] / 2, l / 2], idx)
                        # Line for only square tiles
                        # self.tiles.size = ([min(size) / 2, min(size) / 2, l / 2], idx)
                        self.tiles.set_easy_axis([np.pi/2, 0], idx)
                        self.tiles.mu_r_ea = (self.mu_r, idx)
                        self.tiles.mu_r_oa = (self.mu_r, idx)
                        self.tiles.magnet_type = (soft_type, idx)
                    
                        color = [0.4, 0.4, 0.4]
                        color[z - 3 * (z//3)] = 0.9
                        self.tiles.color = (color, idx)
                        self.tiles.incl_it = (0, idx)
                        idx += 1
                        shim_idx += 1

        self.n_shim = shim_idx
        self._pts_shim = EvaluationPoints([pts_shim.shape[0], 1, 1], self.r_bore, coords=pts_shim)


    def set_shim_matrix(self, shim_mat):
        assert len(shim_mat) == self.n_shim
        self.reset()

        for i in range(self.n_shim):
            if shim_mat[i] == 1: self.tiles.incl_it =(1, self.n_hard_tiles + i)

    
    def set_one_shim_magnet(self, idx):
        assert idx < self.n_hard_tiles
        self.reset(shim_magnets=False)
        self.tiles.incl_it = (1, self.n_hard_tiles + idx)

    
    def reset(self, shim_magnets=True):
        ''' Reset all calculated magnetizations after a change of the Halbach configuration.
        '''
        for idx in range(self.n_hard_tiles):
            self.tiles.M = (0, idx)
            if shim_magnets and idx >= self.n_hard_tiles:
                self.tiles.incl_it = (0, idx)
