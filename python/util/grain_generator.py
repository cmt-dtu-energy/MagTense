import numpy as np
import time
from typing import Optional
from scipy.spatial import KDTree
from scipy.sparse import csr_matrix

class Vertex:
    """
    Represents a 3D coordinate point within the mesh.
    
    Attributes:
        coords (np.array): 1D array of [X, Y, Z] coordinates.
        on_boundary (bool): Flag indicating if the vertex lies on the domain boundary.
    """
    __slots__ = ['coords', 'on_boundary']
    def __init__(self, coords, grid_l, tol=1e-12):
        self.coords = np.array(coords)
        # Check if any coordinate component touches the domain limits (0 or L)
        self.on_boundary = np.any(np.abs(self.coords) < tol) or \
                           np.any(np.abs(self.coords - grid_l) < tol)

class Face:
    """
    Represents a quadrilateral face of a Tile.
    
    Attributes:
        vertices (list): List of 4 Vertex objects defining the corners.
        normal (np.array): Normal vector [nx, ny, nz] pointing outward from the parent tile.
        level (int): Refinement level of the parent tile.
        tile_idx (int): Global index of the tile that owns this face.
        center (np.array): The geometric center of the face.
        bounds (tuple): (min_coords, max_coords) defining the 2D bounding box of the face.
        id (tuple): Geometric identifier used for duplicate pruning.
        plane_key (tuple): Key used to group faces by their common geometric plane.
    """
    __slots__ = ['vertices', 'normal', 'center', 'bounds', 'level', 'tile_idx', 'id', 'plane_key']
    def __init__(self, vertices, normal, level, tile_idx):
        self.vertices = vertices
        self.normal = np.array(normal)
        self.level = level
        self.tile_idx = tile_idx
        
        v_coords = np.array([v.coords for v in vertices])
        self.center = np.mean(v_coords, axis=0)
        self.bounds = (np.min(v_coords, axis=0), np.max(v_coords, axis=0))
        
        # Unique ID for exact twin pruning (same position and size)
        dims = np.round(self.bounds[1] - self.bounds[0], 12)
        self.id = tuple(np.round(self.center, 12)) + tuple(dims)
        
        # Grouping key for faces on the same plane (axis index and rounded coordinate)
        ax = np.argmax(np.abs(self.normal))
        self.plane_key = (ax, np.round(self.center[ax], 12))

class Tile:
    """
    Represents a cubic volume element (cell) in the adaptive mesh.
    
    Attributes:
        pos (np.array): Center position [X, Y, Z].
        dims (np.array): Dimensions [dx, dy, dz].
        level (int): Refinement depth (0 is base grid).
        grain_idx (int): ID of the grain this tile is assigned to.
        faces (list): List of 6 Face objects.
    """
    __slots__ = ['pos', 'dims', 'level', 'grain_idx', 'faces']
    def __init__(self, pos, dims, level):
        self.pos = np.array(pos)
        self.dims = np.array(dims)
        self.level = level
        self.grain_idx = 0
        self.faces = []

class AdaptiveGrainGenerator:
    """
    Generates an adaptive Cartesian mesh based on a Voronoi grain structure.
    Optimized for high-performance plotting and consistent simulation geometry.
    """
    def __init__(self, 
                 n_grains: int = 25, 
                 grid_l: float = 240e-9, 
                 offset_d: float = 3.75e-9, 
                 seed: int = 4, 
                 rng: Optional[np.random.Generator] = None):
        """
        Initializes the generator.
        
        Args:
            n_grains: Number of grains to generate.
            grid_l: Length of the cubic simulation box.
            offset_d: Half-thickness of the intergrain region.
            seed: Random seed (ignored if rng is provided).
            rng: Optional NumPy random generator for reproducibility.
        """
        self.n_grains = n_grains
        self.grid_l = np.array([grid_l]*3 if isinstance(grid_l, (float, int)) else grid_l)
        self.offSetD = offset_d
        
        # Initialize RNG for reproducibility
        self.rng = rng if rng is not None else np.random.default_rng(seed)
        
        # Generate Voronoi seeds
        self.seeds = self.rng.random((n_grains, 3)) * self.grid_l
        self.tree = KDTree(self.seeds)

        # Internal storage
        self.mesh_cart = {}
        self.grid_info = {}
        self.grid_info_plot = {}
        self.mesh_params = {"nGrains": n_grains, "thisGridL": self.grid_l, "offSetD": offset_d}

    def _get_boundary_dist(self, points: np.ndarray):
        """Calculates the distance from points to the nearest grain boundary."""
        dists, idxs = self.tree.query(points, k=2)
        s1, s2 = self.seeds[idxs[:, 0]], self.seeds[idxs[:, 1]]
        # Distance to the bisecting plane between the two closest seeds
        dist_b = (dists[:, 1]**2 - dists[:, 0]**2) / (2 * np.linalg.norm(s1 - s2, axis=1))
        return dist_b, idxs[:, 0]

    def generate(self, res_base: int = 5, num_refinements: int = 3):
        """
        Builds the adaptive mesh by iteratively splitting tiles near grain boundaries.
        """
        start_time = time.time()
        dx, dy, dz = self.grid_l / res_base
        base_size = np.array([dx, dy, dz])
        
        # Recursive tile generation
        tile_list = []
        for cx in np.linspace(dx/2, self.grid_l[0]-dx/2, res_base):
            for cy in np.linspace(dy/2, self.grid_l[1]-dy/2, res_base):
                for cz in np.linspace(dz/2, self.grid_l[2]-dz/2, res_base):
                    self._refine_tile(np.array([cx, cy, cz]), base_size, 0, num_refinements, tile_list)

        # Grain assignment
        pos_arr = np.array([t.pos for t in tile_list])
        dist_b, nearest_grain = self._get_boundary_dist(pos_arr)
        grain_indices = (nearest_grain + 1).astype(int)
        grain_indices[dist_b <= self.offSetD] = self.n_grains + 1
        
        # Setup Face and Vertex objects
        norms = [[-1,0,0], [1,0,0], [0,-1,0], [0,1,0], [0,0,-1], [0,0,1]]
        for i, t in enumerate(tile_list):
            t.grain_idx = int(grain_indices[i])
            for n in norms:
                v_list = []
                ax = np.argmax(np.abs(n))
                others = [a for a in range(3) if a != ax]
                for io, jo in [(-0.5,-0.5), (0.5,-0.5), (0.5,0.5), (-0.5,0.5)]:
                    v_p = np.copy(t.pos)
                    v_p[ax] += n[ax] * t.dims[ax] / 2
                    v_p[others[0]] += io * t.dims[others[0]]
                    v_p[others[1]] += jo * t.dims[others[1]]
                    v_list.append(Vertex(v_p, self.grid_l))
                t.faces.append(Face(v_list, n, t.level, i))

        # Save for simulation use
        self.mesh_cart = {
            "pos_out": pos_arr,
            "dims_out": np.array([t.dims for t in tile_list]),
            "GrainIndex": grain_indices.flatten(),
            "iIn": [np.where(grain_indices == g)[0] for g in range(1, self.n_grains + 2)],
            "tiles": tile_list
        }
        self.mesh_params["resolution"] = [len(tile_list), 1, 1]
        print(f"Mesh generated: {len(tile_list)} tiles in {time.time()-start_time:.2f}s.")

    def _refine_tile(self, pos, size, level, max_level, tile_list):
        """Recursive helper to refine tiles based on proximity to grain boundaries."""
        db, _ = self._get_boundary_dist(pos[np.newaxis, :])
        if level < max_level and db <= (self.offSetD + np.linalg.norm(size)/2):
            new_s = size / 2
            offsets = np.array([[-.5,-.5,-.5], [.5,-.5,-.5], [-.5,.5,-.5], [.5,.5,-.5],
                                [-.5,-.5,.5], [.5,-.5,.5], [-.5,.5,.5], [.5,.5,.5]])
            for off in offsets:
                self._refine_tile(pos + off * size/2, new_s, level + 1, max_level, tile_list)
        else:
            tile_list.append(Tile(pos, size, level))

    def prepare_grid_info(self, ignore_internal: bool = False):
        """
        Processes face geometry to remove overlaps and duplicate twins.
        Generates the standard grid_info and an optional grid_info_plot.
        """
        tiles = self.mesh_cart["tiles"]
        all_faces = [f for t in tiles for f in t.faces]
        tol = 1e-12

        # Group faces by geometric plane
        planes = {}
        for f in all_faces:
            plane_id = (f.plane_key[0], f.plane_key[1])
            if plane_id not in planes: planes[plane_id] = []
            planes[plane_id].append(f)

        # Level-aware pruning: Discard large faces if smaller ones cover them
        final_faces = []
        for p_id, p_faces in planes.items():
            p_kept = []
            for f_check in p_faces:
                is_large_face_split = False
                for f_other in p_faces:
                    if f_other.level > f_check.level:
                        # Check if smaller face f_other is inside larger face f_check
                        if np.all(f_other.bounds[0] >= f_check.bounds[0] - tol) and \
                           np.all(f_other.bounds[1] <= f_check.bounds[1] + tol):
                            is_large_face_split = True
                            break
                
                if not is_large_face_split:
                    # Prune exact mesh twins
                    if not any(f.id == f_check.id for f in p_kept):
                        p_kept.append(f_check)
            final_faces.extend(p_kept)

        # Filter for boundary faces if plotting performance is requested
        plot_faces = [f for f in final_faces if not ignore_internal or any(v.on_boundary for v in f.vertices)]

        def build_struct(f_list):
            n = len(f_list)
            the_signs = csr_matrix(
                (np.ones(n), ([f.tile_idx for f in f_list], np.arange(n))),
                shape=(len(tiles), n)
            )
            return {
                "TheSigns": the_signs,
                "Xf": np.array([f.center[0] for f in f_list]),
                "Yf": np.array([f.center[1] for f in f_list]),
                "Zf": np.array([f.center[2] for f in f_list]),
                "fNormX": np.array([f.normal[0] for f in f_list]),
                "fNormY": np.array([f.normal[1] for f in f_list]),
                "fNormZ": np.array([f.normal[2] for f in f_list])
            }

        self.grid_info = build_struct(final_faces)
        print(f"Grid info ready. Total faces: {len(final_faces)}")
        
        if ignore_internal:
            self.grid_info_plot = build_struct(plot_faces)
            print(f"Plot faces (on boundary): {len(plot_faces)}")