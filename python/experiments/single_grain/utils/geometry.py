"""Conservative adaptive Cartesian filling for analytic 3D domains.

The filler accepts only axis-aligned cubes that can be proven to lie fully
inside the target domain.  For signed-distance domains the proof uses the
inscribed-domain criterion

    sdf(center) >= sqrt(3) * h / 2 + eps,

where ``h`` is the cube side length.  Ambiguous boundary cells are refined or
left unresolved; they are never silently accepted.
"""

from __future__ import annotations

from collections import deque
from dataclasses import asdict, dataclass, field
from enum import IntEnum
from heapq import heappop, heappush
from itertools import count, product
from math import pi, sqrt
from typing import Any, Mapping, Protocol, Sequence

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


class CubeState(IntEnum):
    """Conservative cube classification."""

    OUTSIDE = 0
    INSIDE = 1
    BOUNDARY = 2


@dataclass(frozen=True)
class Cube:
    """A dyadic axis-aligned cube accepted or refined by the filler."""

    center: FloatArray
    h: float
    level: int
    i: int
    j: int
    k: int

    def __post_init__(self) -> None:
        center = _vector3(self.center, "center")
        if not np.isfinite(self.h) or self.h <= 0:
            raise ValueError("h must be positive and finite")
        if self.level < 0:
            raise ValueError("level must be non-negative")
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "h", float(self.h))
        object.__setattr__(self, "level", int(self.level))
        object.__setattr__(self, "i", int(self.i))
        object.__setattr__(self, "j", int(self.j))
        object.__setattr__(self, "k", int(self.k))


@dataclass(frozen=True)
class CubeFillStats:
    """Metadata describing an adaptive fill run."""

    accepted_count: int
    accepted_volume: float
    fill_fraction: float | None
    max_level_reached: int
    classified_cells: int
    rejected_outside_cells: int
    boundary_cells_left_unresolved: int
    termination_reason: str
    boundary_cells_queued: int
    h0: float
    h_min: float
    max_depth: int
    eps: float
    overshoot_policy: str
    queue_policy: str
    complete_layers: bool = False
    target_cells: int | None = None
    target_lower_tolerance: float = 0.0
    target_accepted: bool | None = None
    target_status: str = "not_requested"
    relative_count_error: float | None = None
    layers_used: int = 0
    h_min_reached: bool = False
    min_inside_fraction: float = 1.0
    inside_fraction_samples: int = 1
    selected_grid_shift: tuple[float, float, float] = (0.0, 0.0, 0.0)
    grid_shift_candidates: tuple[Mapping[str, Any], ...] = ()

    def as_dict(self) -> dict[str, Any]:
        """Return a plain-Python representation suitable for metadata."""
        return asdict(self)


@dataclass(frozen=True)
class H0CalibrationRecord:
    """One reusable calibration row for choosing a coarse starting size.

    ``h0_multiplier`` is dimensionless and is applied to the analytic volume
    estimate ``(eta * volume / N_target) ** (1 / 3)``.  The quality fields are
    copied from the fill run that produced the recommendation, so notebooks can
    show why a row was selected rather than treating it as a magic constant.
    """

    shape_key: str
    target_min: int
    target_max: int
    h0_multiplier: float
    fill_fraction: float | None
    achieved_count: int | None
    max_level_reached: int | None
    unresolved_boundary_cells: int | None
    classified_cells: int | None

    def as_dict(self) -> dict[str, Any]:
        """Return a plain-Python representation for tables and metadata."""
        return asdict(self)


H0_REFERENCE_TABLE: Mapping[str, tuple[H0CalibrationRecord, ...]] = {
    # These seed rows are deliberately dimensionless.  They provide a reusable
    # starting point for future notebooks/scripts, while ``calibrate_h0_sweep``
    # can regenerate project-specific rows for a particular shape set.  These
    # bands are intentionally limited to the target ranges covered by the
    # calibration notebooks; larger targets fall back to the analytic estimate
    # unless a new calibration row is added.
    "sphere": (
        H0CalibrationRecord("sphere", 1, 512, 0.5, None, None, None, None, None),
        H0CalibrationRecord("sphere", 513, 4096, 0.5, None, None, None, None, None),
    ),
    "ellipsoid": (
        H0CalibrationRecord("ellipsoid", 1, 512, 0.5, None, None, None, None, None),
        H0CalibrationRecord("ellipsoid", 513, 4096, 0.5, None, None, None, None, None),
    ),
    "cylinder": (
        H0CalibrationRecord("cylinder", 1, 512, 0.5, None, None, None, None, None),
        H0CalibrationRecord("cylinder", 513, 4096, 0.5, None, None, None, None, None),
    ),
    "hexagonal_prism": (
        H0CalibrationRecord("hexagonal_prism", 1, 512, 0.5, None, None, None, None, None),
        H0CalibrationRecord("hexagonal_prism", 513, 4096, 0.5, None, None, None, None, None),
    ),
}


@dataclass(frozen=True)
class PrismMesh:
    """Axis-aligned prism cells and the metadata used to generate them."""

    centers: FloatArray
    dimensions: FloatArray
    levels: IntArray
    shape: str
    shape_metadata: Mapping[str, Any]
    root_bounds: FloatArray
    target_tiles: int | None = None
    refinement_metadata: Mapping[str, Any] = field(default_factory=dict)
    diagnostics: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        centers = np.asarray(self.centers, dtype=float)
        dimensions = np.asarray(self.dimensions, dtype=float)
        levels = np.asarray(self.levels, dtype=np.int64)
        root_bounds = np.asarray(self.root_bounds, dtype=float)

        if centers.ndim != 2 or centers.shape[1] != 3:
            raise ValueError("centers must have shape (number_of_cells, 3)")
        if dimensions.shape != centers.shape:
            raise ValueError("dimensions must have the same shape as centers")
        if levels.shape != (len(centers),):
            raise ValueError("levels must have one value per cell")
        if root_bounds.shape != (2, 3):
            raise ValueError("root_bounds must have shape (2, 3)")
        if np.any(dimensions <= 0):
            raise ValueError("all cell dimensions must be positive")

        object.__setattr__(self, "centers", centers)
        object.__setattr__(self, "dimensions", dimensions)
        object.__setattr__(self, "levels", levels)
        object.__setattr__(self, "root_bounds", root_bounds)

    @property
    def achieved_tiles(self) -> int:
        """Number of cells in the generated mesh."""
        return len(self.centers)

    @property
    def represented_volume(self) -> float:
        """Total volume represented by all retained cells."""
        return float(np.sum(np.prod(self.dimensions, axis=1)))

    @property
    def target_reached(self) -> bool | None:
        """Whether the mesh satisfies the requested target semantics."""
        if self.target_tiles is None:
            return None
        stats = self.refinement_metadata.get("cube_fill_stats", {})
        if isinstance(stats, Mapping) and bool(stats.get("complete_layers", False)):
            accepted = stats.get("target_accepted")
            if accepted is not None:
                return bool(accepted)
        return self.achieved_tiles <= self.target_tiles

    @property
    def target_utilization(self) -> float | None:
        """Fraction of the requested tile target used by the mesh."""
        if self.target_tiles is None:
            return None
        return self.achieved_tiles / self.target_tiles

    @property
    def target_filled(self) -> bool | None:
        """Whether the generated mesh contains exactly the requested tile count."""
        if self.target_tiles is None:
            return None
        return self.achieved_tiles == self.target_tiles

    def to_micromag_kwargs(self) -> dict[str, Any]:
        """Return the grid arguments expected by :class:`MicromagProblem`."""
        return {
            "grid_type": "unstructuredPrisms",
            "grid_pts": self.centers.copy(),
            "grid_abc": self.dimensions.copy(),
            "res": [self.achieved_tiles, 1, 1],
        }


class DomainOracle(Protocol):
    """Domain interface used by :func:`adaptive_cube_fill`.

    Mesh/STL support can be added by implementing this protocol with a BVH-
    backed signed-distance query and, if needed, a custom ``classify_cube``.
    """

    name: str
    volume: float | None
    bounding_box: FloatArray | None

    def sdf(self, point: Sequence[float]) -> float:
        """Return signed distance with positive values inside the domain."""
        ...

    def classify_cube(self, cube: Cube, eps: float) -> CubeState:
        """Classify a cube conservatively."""
        ...


class SignedDistanceDomain:
    """Base class for analytic domains with a valid signed-distance function."""

    name = "signed_distance_domain"
    volume: float | None = None
    bounding_box: FloatArray | None = None

    def sdf(self, point: Sequence[float]) -> float:
        """Return signed distance with positive values inside the domain."""
        raise NotImplementedError

    def classify_cube(self, cube: Cube, eps: float) -> CubeState:
        """Classify using the conservative SDF sphere bound around a cube."""
        radius = sqrt(3.0) * cube.h / 2.0
        value = self.sdf(cube.center)
        if value >= radius + eps:
            return CubeState.INSIDE
        if value <= -radius - eps:
            return CubeState.OUTSIDE
        return CubeState.BOUNDARY


@dataclass(frozen=True)
class Sphere(SignedDistanceDomain):
    """Sphere domain with exact signed distance."""

    center: FloatArray
    radius: float
    name: str = "sphere"

    def __post_init__(self) -> None:
        center = _vector3(self.center, "center")
        if not np.isfinite(self.radius) or self.radius <= 0:
            raise ValueError("radius must be positive and finite")
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "radius", float(self.radius))
        bounds = np.vstack((center - self.radius, center + self.radius))
        object.__setattr__(self, "bounding_box", bounds)
        object.__setattr__(self, "volume", 4.0 * pi * self.radius**3 / 3.0)

    def sdf(self, point: Sequence[float]) -> float:
        """Return exact sphere signed distance, positive inside."""
        return float(self.radius - np.linalg.norm(_vector3(point, "point") - self.center))


@dataclass(frozen=True)
class AxisAlignedBox(SignedDistanceDomain):
    """Axis-aligned box with exact signed distance."""

    center: FloatArray
    dimensions: FloatArray
    name: str = "axis_aligned_box"

    def __post_init__(self) -> None:
        center = _vector3(self.center, "center")
        dimensions = _vector3(self.dimensions, "dimensions", positive=True)
        half = dimensions / 2.0
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "dimensions", dimensions)
        object.__setattr__(self, "bounding_box", np.vstack((center - half, center + half)))
        object.__setattr__(self, "volume", float(np.prod(dimensions)))

    def sdf(self, point: Sequence[float]) -> float:
        """Return exact box signed distance, positive inside."""
        q = np.abs(_vector3(point, "point") - self.center) - self.dimensions / 2.0
        outside = np.linalg.norm(np.maximum(q, 0.0))
        inside = min(float(np.max(q)), 0.0)
        return float(-(outside + inside))

    def classify_cube(self, cube: Cube, eps: float) -> CubeState:
        """Classify boxes exactly by interval containment/separation."""
        lower, upper = self.bounding_box
        cube_lower = cube.center - cube.h / 2.0
        cube_upper = cube.center + cube.h / 2.0
        if np.any(cube_upper < lower - eps) or np.any(cube_lower > upper + eps):
            return CubeState.OUTSIDE
        if np.all(cube_lower >= lower + eps) and np.all(cube_upper <= upper - eps):
            return CubeState.INSIDE
        if np.all(cube_lower >= lower - eps) and np.all(cube_upper <= upper + eps):
            return CubeState.INSIDE
        return CubeState.BOUNDARY


@dataclass(frozen=True)
class Ellipsoid(SignedDistanceDomain):
    """Axis-aligned ellipsoid.

    The ``sdf`` method returns a conservative lower bound on the signed distance
    by scaling the implicit margin with the largest semi-axis.  Cube
    classification uses a safer interval test on the quadratic implicit
    function, so accepted cells are not based on this approximate SDF.
    """

    center: FloatArray
    semi_axes: FloatArray
    name: str = "ellipsoid"

    def __post_init__(self) -> None:
        center = _vector3(self.center, "center")
        semi_axes = _vector3(self.semi_axes, "semi_axes", positive=True)
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "semi_axes", semi_axes)
        object.__setattr__(
            self,
            "bounding_box",
            np.vstack((center - semi_axes, center + semi_axes)),
        )
        object.__setattr__(self, "volume", 4.0 * pi * float(np.prod(semi_axes)) / 3.0)

    def sdf(self, point: Sequence[float]) -> float:
        """Return a conservative signed-distance lower bound, positive inside."""
        relative = (_vector3(point, "point") - self.center) / self.semi_axes
        implicit_margin = 1.0 - float(np.linalg.norm(relative))
        return implicit_margin * float(np.min(self.semi_axes))

    def classify_cube(self, cube: Cube, eps: float) -> CubeState:
        """Classify conservatively using interval bounds of the ellipsoid form."""
        half = cube.h / 2.0
        relative = cube.center - self.center
        nearest = np.maximum(np.abs(relative) - half, 0.0)
        min_q = float(np.sum((nearest / self.semi_axes) ** 2))
        if min_q > 1.0 + eps:
            return CubeState.OUTSIDE

        farthest = np.abs(relative) + half
        max_q = float(np.sum((farthest / self.semi_axes) ** 2))
        if max_q <= 1.0 - eps:
            return CubeState.INSIDE
        if max_q <= 1.0 and eps == 0.0:
            return CubeState.INSIDE
        return CubeState.BOUNDARY


@dataclass(frozen=True)
class Cylinder(SignedDistanceDomain):
    """Finite axis-aligned cylinder with exact signed distance.

    Parameters
    ----------
    center:
        Cylinder center as ``[x, y, z]``.
    radius:
        Radius in the two directions perpendicular to ``axis``.
    length:
        Full cylinder length along ``axis``.
    axis:
        One of ``"x"``, ``"y"``, or ``"z"``.  Only axis-aligned cylinders are
        supported by this conservative Cartesian filler.

    Cube classification is deliberately interval based: a cube is accepted
    only when its full axial interval lies between the end caps and its
    farthest radial corner is still inside the circular cross-section.
    """

    center: FloatArray
    radius: float
    length: float
    axis: str
    name: str = "cylinder"

    def __post_init__(self) -> None:
        center = _vector3(self.center, "center")
        if not np.isfinite(self.radius) or self.radius <= 0:
            raise ValueError("radius must be positive and finite")
        if not np.isfinite(self.length) or self.length <= 0:
            raise ValueError("length must be positive and finite")
        axis_index = _axis_index(self.axis)
        half_extents = np.full(3, float(self.radius))
        half_extents[axis_index] = float(self.length) / 2.0
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "radius", float(self.radius))
        object.__setattr__(self, "length", float(self.length))
        object.__setattr__(self, "axis", self.axis)
        object.__setattr__(self, "_axis_index", axis_index)
        object.__setattr__(
            self,
            "_radial_indices",
            tuple(index for index in range(3) if index != axis_index),
        )
        object.__setattr__(
            self,
            "bounding_box",
            np.vstack((center - half_extents, center + half_extents)),
        )
        object.__setattr__(self, "volume", pi * self.radius**2 * self.length)

    def sdf(self, point: Sequence[float]) -> float:
        """Return exact finite-cylinder signed distance, positive inside."""
        relative = _vector3(point, "point") - self.center
        radial = np.linalg.norm(relative[list(self._radial_indices)]) - self.radius
        axial = abs(relative[self._axis_index]) - self.length / 2.0
        q = np.array([radial, axial], dtype=float)
        outside = np.linalg.norm(np.maximum(q, 0.0))
        inside = min(float(np.max(q)), 0.0)
        return float(-(outside + inside))

    def classify_cube(self, cube: Cube, eps: float) -> CubeState:
        """Classify conservatively with axial and radial interval bounds."""
        half = cube.h / 2.0
        relative = cube.center - self.center
        half_length = self.length / 2.0
        axial_abs = abs(relative[self._axis_index])
        if axial_abs - half > half_length + eps:
            return CubeState.OUTSIDE

        radial_relative = relative[list(self._radial_indices)]
        nearest_radial = np.maximum(np.abs(radial_relative) - half, 0.0)
        if float(np.linalg.norm(nearest_radial)) > self.radius + eps:
            return CubeState.OUTSIDE

        farthest_radial = np.abs(radial_relative) + half
        if (
            axial_abs + half <= half_length - eps
            and float(np.linalg.norm(farthest_radial)) <= self.radius - eps
        ):
            return CubeState.INSIDE
        if (
            eps == 0.0
            and axial_abs + half <= half_length
            and float(np.linalg.norm(farthest_radial)) <= self.radius
        ):
            return CubeState.INSIDE
        return CubeState.BOUNDARY


@dataclass(frozen=True)
class HexagonalPrism(SignedDistanceDomain):
    """Finite regular hexagonal prism with conservative half-space tests.

    Parameters
    ----------
    center:
        Prism center as ``[x, y, z]``.
    side_length:
        Regular-hexagon side length.  For a regular hexagon this is also the
        center-to-vertex radius in the radial plane.
    height:
        Full prism height along ``axis``.
    axis:
        One of ``"x"``, ``"y"``, or ``"z"``.
    rotation_degrees:
        In-plane rotation of the hexagon around ``axis``.

    The radial hexagon is represented as six half-spaces.  For each candidate
    cube, the classifier computes the cube's support interval against every
    half-space normal.  A cube is accepted only if the whole support interval
    is inside every face and inside both axial end caps.
    """

    center: FloatArray
    side_length: float
    height: float
    axis: str
    rotation_degrees: float = 0.0
    name: str = "hexagonal_prism"

    def __post_init__(self) -> None:
        center = _vector3(self.center, "center")
        if not np.isfinite(self.side_length) or self.side_length <= 0:
            raise ValueError("side_length must be positive and finite")
        if not np.isfinite(self.height) or self.height <= 0:
            raise ValueError("height must be positive and finite")
        if not np.isfinite(self.rotation_degrees):
            raise ValueError("rotation_degrees must be finite")
        axis_index = _axis_index(self.axis)
        radial_indices = tuple(index for index in range(3) if index != axis_index)
        rotation = np.deg2rad(self.rotation_degrees)
        vertex_angles = rotation + np.arange(6) * pi / 3.0
        radial_vertices = self.side_length * np.column_stack(
            (np.cos(vertex_angles), np.sin(vertex_angles))
        )
        normal_angles = rotation + (np.arange(6) + 0.5) * pi / 3.0
        face_normals = np.column_stack((np.cos(normal_angles), np.sin(normal_angles)))
        vertices = np.tile(center, (6, 1))
        vertices[:, list(radial_indices)] += radial_vertices
        lower = np.min(vertices, axis=0)
        upper = np.max(vertices, axis=0)
        lower[axis_index] = center[axis_index] - self.height / 2.0
        upper[axis_index] = center[axis_index] + self.height / 2.0
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "side_length", float(self.side_length))
        object.__setattr__(self, "height", float(self.height))
        object.__setattr__(self, "axis", self.axis)
        object.__setattr__(self, "rotation_degrees", float(self.rotation_degrees))
        object.__setattr__(self, "_axis_index", axis_index)
        object.__setattr__(self, "_radial_indices", radial_indices)
        object.__setattr__(self, "_face_normals", face_normals)
        object.__setattr__(self, "_apothem", sqrt(3.0) * self.side_length / 2.0)
        object.__setattr__(self, "_vertices", vertices)
        object.__setattr__(self, "bounding_box", np.vstack((lower, upper)))
        object.__setattr__(
            self,
            "volume",
            3.0 * sqrt(3.0) * self.side_length**2 * self.height / 2.0,
        )

    def sdf(self, point: Sequence[float]) -> float:
        """Return a conservative signed-distance lower bound, positive inside."""
        relative = _vector3(point, "point") - self.center
        radial = relative[list(self._radial_indices)]
        radial_margin = self._apothem - float(np.max(self._face_normals @ radial))
        axial_margin = self.height / 2.0 - abs(relative[self._axis_index])
        return min(radial_margin, axial_margin)

    def classify_cube(self, cube: Cube, eps: float) -> CubeState:
        """Classify via support intervals against all prism half-spaces."""
        half = cube.h / 2.0
        relative = cube.center - self.center
        half_height = self.height / 2.0
        axial_abs = abs(relative[self._axis_index])
        if axial_abs - half > half_height + eps:
            return CubeState.OUTSIDE

        radial = relative[list(self._radial_indices)]
        projections = self._face_normals @ radial
        supports = half * np.sum(np.abs(self._face_normals), axis=1)
        if np.any(projections - supports > self._apothem + eps):
            return CubeState.OUTSIDE

        if (
            axial_abs + half <= half_height - eps
            and np.all(projections + supports <= self._apothem - eps)
        ):
            return CubeState.INSIDE
        if (
            eps == 0.0
            and axial_abs + half <= half_height
            and np.all(projections + supports <= self._apothem)
        ):
            return CubeState.INSIDE
        return CubeState.BOUNDARY


def adaptive_cube_fill(
    domain: DomainOracle,
    bbox: Sequence[Sequence[float]] | None = None,
    N_target: int = 1000,
    h0: float | None = None,
    h_min: float = 0.0,
    max_depth: int = 12,
    eps: float = 1e-12,
    overshoot_policy: str = "closest",
    queue_policy: str = "breadth_first",
    eta: float = 0.75,
    complete_layers: bool = True,
    target_lower_tolerance: float = 0.10,
    grid_shifts: Sequence[Sequence[float]] | str | None = None,
    min_inside_fraction: float = 1.0,
    inside_fraction_samples: int = 3,
) -> tuple[list[Cube], CubeFillStats]:
    """Fill ``domain`` with non-overlapping conservative dyadic cubes.

    Parameters
    ----------
    domain:
        Any object implementing ``DomainOracle``.  Analytic shapes in this
        module provide conservative ``classify_cube`` methods.
    bbox:
        Optional ``[[xmin, ymin, zmin], [xmax, ymax, zmax]]``.  If omitted,
        the domain must provide ``bounding_box``.
    N_target:
        Approximate desired number of accepted cubes.  The algorithm may stop
        below or above this depending on ``overshoot_policy``; containment is
        always more important than exact count.
    h0:
        Starting coarse cube side length.  Larger values produce larger
        interior cubes and more visible boundary refinement.  If ``None``,
        ``h0`` is estimated from domain volume and ``N_target``.
    h_min:
        Smallest side length allowed during boundary refinement.  Boundary
        cubes at or below this size are left unresolved, never accepted.
    max_depth:
        Maximum dyadic refinement level.  Level 0 cubes have side ``h0``,
        level 1 cubes have side ``h0 / 2``, and so on.
    eps:
        Numerical safety tolerance used by conservative classifiers.
    overshoot_policy:
        ``"soft"`` stops after reaching/exceeding ``N_target``; useful for
        visualizing adaptive refinement.  ``"never_exceed"`` rejects
        refinement batches that would cross ``N_target``.  It does not remove
        coarse level-0 cubes that were already conservatively accepted.
        ``"closest"`` accepts refinement steps only when they improve the
        count error.
    queue_policy:
        ``"breadth_first"`` refines coarse boundary cells first. ``"priority"``
        uses a simple larger-cells-first priority queue.
        ``"symmetric_breadth_first"`` and ``"symmetric_priority"`` refine
        mirrored boundary cells as a batch around ``domain.center`` when the
        domain provides one.  Symmetric batches preserve a more balanced mesh
        when a target count stops refinement early, but a batch can overshoot
        ``N_target`` more than a single-parent refinement step.  The symmetric
        priority mode also favors radial/off-axis boundary cells for
        axis-aligned cylinders and hexagonal prisms, because axis-aligned cubes
        already fit most naturally along the shape axis.
    eta:
        Dimensionless fill-efficiency factor for estimating ``h0`` when no
        explicit ``h0`` is provided.  The estimate is
        ``(eta * domain.volume / N_target) ** (1 / 3)``.  Smaller ``eta``
        means a finer starting grid; larger ``eta`` means a coarser one.
    complete_layers:
        When true, refine all boundary cells at a given level before checking
        the target-count criterion.  Queue ordering and symmetry batching do
        not select cells in this mode.  When false, use the legacy target-limited
        transaction logic controlled by ``overshoot_policy``.
    target_lower_tolerance:
        Symmetric fractional target tolerance used in complete-layer mode.
        With the default value, a completed layer is accepted only between
        ``0.9 * N_target`` and ``1.1 * N_target`` cells.  A layer that would
        exceed the upper bound is discarded in full.
    grid_shifts:
        Optional run-level grid shifts.  ``None`` uses the unshifted grid.
        ``"half_step"`` tries all combinations of ``-0.5*h0, 0, 0.5*h0`` in
        x/y/z and selects the best completed run.
    min_inside_fraction:
        Required fraction of sampled cube volume inside the domain.  The
        default ``1.0`` uses the existing strict conservative classifiers.
        Values below one use approximate SDF sub-sampling.
    inside_fraction_samples:
        Number of sample points per coordinate for approximate partial-tile
        acceptance when ``min_inside_fraction < 1.0``.
    """
    _validate_target(N_target)
    bbox_array = _resolve_bbox(domain, bbox)
    h0_value = _resolve_h0(domain, bbox_array, N_target, h0, eta)
    h_min_value = _validate_nonnegative_float(h_min, "h_min")
    max_depth_value = _validate_max_depth(max_depth)
    eps_value = _validate_nonnegative_float(eps, "eps")
    target_lower_tolerance_value = _validate_target_lower_tolerance(target_lower_tolerance)
    min_inside_fraction_value = _validate_min_inside_fraction(min_inside_fraction)
    inside_fraction_samples_value = _validate_sample_count(inside_fraction_samples)
    if overshoot_policy not in {"soft", "never_exceed", "closest"}:
        raise ValueError("overshoot_policy must be 'soft', 'never_exceed', or 'closest'")
    queue_policies = {
        "breadth_first",
        "priority",
        "symmetric_breadth_first",
        "symmetric_priority",
    }
    if queue_policy not in queue_policies:
        raise ValueError(
            "queue_policy must be 'breadth_first', 'priority', "
            "'symmetric_breadth_first', or 'symmetric_priority'"
        )

    shift_candidates = _resolve_grid_shifts(grid_shifts, h0_value)
    runs: list[tuple[list[Cube], CubeFillStats]] = []
    for shift in shift_candidates:
        runs.append(
            _adaptive_cube_fill_single_shift(
                domain=domain,
                bbox_array=bbox_array,
                N_target=N_target,
                h0_value=h0_value,
                h_min_value=h_min_value,
                max_depth_value=max_depth_value,
                eps_value=eps_value,
                overshoot_policy=overshoot_policy,
                queue_policy=queue_policy,
                complete_layers=bool(complete_layers),
                target_lower_tolerance_value=target_lower_tolerance_value,
                min_inside_fraction_value=min_inside_fraction_value,
                inside_fraction_samples_value=inside_fraction_samples_value,
                grid_shift=shift,
            )
        )

    selected_cubes, selected_stats = min(runs, key=lambda item: _fill_run_sort_key(item[1]))
    if len(runs) == 1:
        return selected_cubes, selected_stats

    candidate_summaries = tuple(_grid_shift_summary(stats) for _cubes, stats in runs)
    selected_stats = _copy_stats_with_grid_shift_candidates(
        selected_stats,
        candidate_summaries,
    )
    return selected_cubes, selected_stats


def _adaptive_cube_fill_single_shift(
    *,
    domain: DomainOracle,
    bbox_array: FloatArray,
    N_target: int,
    h0_value: float,
    h_min_value: float,
    max_depth_value: int,
    eps_value: float,
    overshoot_policy: str,
    queue_policy: str,
    complete_layers: bool,
    target_lower_tolerance_value: float,
    min_inside_fraction_value: float,
    inside_fraction_samples_value: int,
    grid_shift: FloatArray,
) -> tuple[list[Cube], CubeFillStats]:
    """Run one fill for one resolved grid shift."""
    accepted: list[Cube] = []
    rejected_outside = 0
    classified = 0
    unresolved = 0
    max_level_reached = 0
    termination_reason = "boundary queue empty"
    h_min_reached = False
    layers_used = 0

    # Level-0 cubes are the coarse grid.  Interior level-0 cubes are accepted
    # immediately and are never split later, which is how the output keeps
    # larger cubes away from the boundary.
    boundary_cells: list[Cube] = []
    boundary_queue = _BoundaryQueue(queue_policy)
    for cube in _coarse_grid(bbox_array, h0_value, grid_shift=grid_shift):
        state = _classify_candidate_cube(
            domain,
            cube,
            eps_value,
            min_inside_fraction_value,
            inside_fraction_samples_value,
        )
        classified += 1
        if state == CubeState.INSIDE:
            accepted.append(cube)
        elif state == CubeState.BOUNDARY:
            if complete_layers:
                boundary_cells.append(cube)
            else:
                boundary_queue.push(cube, domain)
        else:
            rejected_outside += 1

    if complete_layers:
        target_accepted = _target_accepted(
            len(accepted),
            N_target,
            target_lower_tolerance_value,
        )
        if target_accepted:
            termination_reason = "target tolerance reached"

        target_upper_bound = (1.0 + target_lower_tolerance_value) * N_target
        if len(accepted) > target_upper_bound:
            termination_reason = "level-0 count exceeds target tolerance"

        while (
            boundary_cells
            and not target_accepted
            and len(accepted) <= target_upper_bound
        ):
            parent_level = min(cube.level for cube in boundary_cells)
            parent_cells = [cube for cube in boundary_cells if cube.level == parent_level]
            boundary_cells = [cube for cube in boundary_cells if cube.level != parent_level]
            staged_inside: list[Cube] = []
            staged_boundary: list[Cube] = []
            processed_any = False
            processed_parent_count = 0

            for cube in sorted(parent_cells, key=lambda cell: (cell.i, cell.j, cell.k)):
                max_level_reached = max(max_level_reached, cube.level)
                if cube.h <= h_min_value:
                    unresolved += 1
                    h_min_reached = True
                    termination_reason = "h_min reached"
                    continue
                if cube.level >= max_depth_value:
                    unresolved += 1
                    termination_reason = "max_depth reached"
                    continue

                processed_any = True
                processed_parent_count += 1
                for child in _subdivide(cube):
                    state = _classify_candidate_cube(
                        domain,
                        child,
                        eps_value,
                        min_inside_fraction_value,
                        inside_fraction_samples_value,
                    )
                    classified += 1
                    max_level_reached = max(max_level_reached, child.level)
                    if state == CubeState.INSIDE:
                        staged_inside.append(child)
                    elif state == CubeState.BOUNDARY:
                        staged_boundary.append(child)
                    else:
                        rejected_outside += 1

            proposed_count = len(accepted) + len(staged_inside)
            if proposed_count > target_upper_bound:
                unresolved += processed_parent_count
                termination_reason = "next complete layer exceeds target tolerance"
                break

            accepted.extend(staged_inside)
            if processed_any:
                layers_used = max(layers_used, parent_level + 1)
            boundary_cells.extend(staged_boundary)
            target_accepted = _target_accepted(
                len(accepted),
                N_target,
                target_lower_tolerance_value,
            )
            if target_accepted:
                termination_reason = "target tolerance reached"
                break

        unresolved += len(boundary_cells)
        if not target_accepted and not boundary_cells and termination_reason == "boundary queue empty":
            termination_reason = "boundary queue empty"
        return _finalize_fill_stats(
            accepted=accepted,
            domain=domain,
            classified=classified,
            rejected_outside=rejected_outside,
            unresolved=unresolved,
            termination_reason=termination_reason,
            boundary_cells_queued=len(boundary_cells),
            h0_value=h0_value,
            h_min_value=h_min_value,
            max_depth_value=max_depth_value,
            eps_value=eps_value,
            overshoot_policy=overshoot_policy,
            queue_policy=queue_policy,
            complete_layers=True,
            N_target=N_target,
            target_lower_tolerance_value=target_lower_tolerance_value,
            layers_used=layers_used,
            h_min_reached=h_min_reached,
            min_inside_fraction_value=min_inside_fraction_value,
            inside_fraction_samples_value=inside_fraction_samples_value,
            grid_shift=grid_shift,
            max_level_reached=max_level_reached,
        )

    # Only ambiguous boundary cubes enter the queue.  Refinement therefore
    # concentrates resolution where the classifier cannot prove containment.
    while boundary_queue:
        if len(accepted) >= N_target and overshoot_policy == "soft":
            termination_reason = "target reached"
            break

        parent_batch = boundary_queue.pop_batch()
        eligible_parents: list[Cube] = []
        for cube in parent_batch:
            max_level_reached = max(max_level_reached, cube.level)
            if cube.h <= h_min_value:
                unresolved += 1
                h_min_reached = True
                termination_reason = "h_min reached"
            elif cube.level >= max_depth_value:
                unresolved += 1
                termination_reason = "max_depth reached"
            else:
                eligible_parents.append(cube)

        if not eligible_parents:
            continue

        # Stage one parent subdivision as a small transaction.  Under the
        # symmetric queue policies, this transaction is a mirror-orbit batch:
        # all mirrored parents are accepted/rejected/requeued together so the
        # target count cannot refine only one side of a symmetric domain.
        staged_inside: list[Cube] = []
        staged_boundary: list[Cube] = []
        for cube in eligible_parents:
            for child in _subdivide(cube):
                state = _classify_candidate_cube(
                    domain,
                    child,
                    eps_value,
                    min_inside_fraction_value,
                    inside_fraction_samples_value,
                )
                classified += 1
                max_level_reached = max(max_level_reached, child.level)
                if state == CubeState.INSIDE:
                    staged_inside.append(child)
                elif state == CubeState.BOUNDARY:
                    staged_boundary.append(child)
                else:
                    rejected_outside += 1

        new_count = len(accepted) + len(staged_inside)
        if overshoot_policy == "never_exceed" and new_count > N_target:
            unresolved += len(eligible_parents) + len(staged_boundary)
            termination_reason = "target would be exceeded"
            continue

        if overshoot_policy == "closest":
            old_error = abs(N_target - len(accepted))
            new_error = abs(N_target - new_count)
            if new_count > N_target and new_error > old_error:
                unresolved += len(eligible_parents) + len(staged_boundary)
                termination_reason = "target closest count reached"
                continue

        accepted.extend(staged_inside)
        for child in staged_boundary:
            boundary_queue.push(child, domain)

    unresolved += len(boundary_queue)
    if boundary_queue and termination_reason == "boundary queue empty":
        termination_reason = "stopped with unresolved boundary cells"

    return _finalize_fill_stats(
        accepted=accepted,
        domain=domain,
        classified=classified,
        rejected_outside=rejected_outside,
        unresolved=unresolved,
        termination_reason=termination_reason,
        boundary_cells_queued=len(boundary_queue),
        h0_value=h0_value,
        h_min_value=h_min_value,
        max_depth_value=max_depth_value,
        eps_value=eps_value,
        overshoot_policy=overshoot_policy,
        queue_policy=queue_policy,
        complete_layers=False,
        N_target=N_target,
        target_lower_tolerance_value=target_lower_tolerance_value,
        layers_used=max_level_reached,
        h_min_reached=h_min_reached,
        min_inside_fraction_value=min_inside_fraction_value,
        inside_fraction_samples_value=inside_fraction_samples_value,
        grid_shift=grid_shift,
        max_level_reached=max_level_reached,
    )


def _finalize_fill_stats(
    *,
    accepted: list[Cube],
    domain: DomainOracle,
    classified: int,
    rejected_outside: int,
    unresolved: int,
    termination_reason: str,
    boundary_cells_queued: int,
    h0_value: float,
    h_min_value: float,
    max_depth_value: int,
    eps_value: float,
    overshoot_policy: str,
    queue_policy: str,
    complete_layers: bool,
    N_target: int,
    target_lower_tolerance_value: float,
    layers_used: int,
    h_min_reached: bool,
    min_inside_fraction_value: float,
    inside_fraction_samples_value: int,
    grid_shift: FloatArray,
    max_level_reached: int,
) -> tuple[list[Cube], CubeFillStats]:
    """Sort accepted cells and build stats for one fill run."""
    accepted.sort(key=lambda cell: (cell.level, cell.i, cell.j, cell.k))
    accepted_volume = float(sum(cube.h**3 for cube in accepted))
    domain_volume = getattr(domain, "volume", None)
    fill_fraction = None if domain_volume in (None, 0.0) else accepted_volume / float(domain_volume)
    if accepted:
        max_level_reached = max(max_level_reached, max(cube.level for cube in accepted))
    target_accepted = _target_accepted(
        len(accepted),
        N_target,
        target_lower_tolerance_value,
    )

    stats = CubeFillStats(
        accepted_count=len(accepted),
        accepted_volume=accepted_volume,
        fill_fraction=fill_fraction,
        max_level_reached=max_level_reached,
        classified_cells=classified,
        rejected_outside_cells=rejected_outside,
        boundary_cells_left_unresolved=unresolved,
        termination_reason=termination_reason,
        boundary_cells_queued=boundary_cells_queued,
        h0=h0_value,
        h_min=h_min_value,
        max_depth=max_depth_value,
        eps=eps_value,
        overshoot_policy=overshoot_policy,
        queue_policy=queue_policy,
        complete_layers=complete_layers,
        target_cells=N_target,
        target_lower_tolerance=target_lower_tolerance_value,
        target_accepted=target_accepted,
        target_status=_target_status(len(accepted), N_target, target_lower_tolerance_value),
        relative_count_error=(len(accepted) - N_target) / N_target,
        layers_used=layers_used,
        h_min_reached=h_min_reached,
        min_inside_fraction=min_inside_fraction_value,
        inside_fraction_samples=inside_fraction_samples_value,
        selected_grid_shift=tuple(float(value) for value in grid_shift),
    )
    return accepted, stats


class _BoundaryQueue:
    """Boundary refinement queue with optional mirror-orbit batching.

    The non-symmetric policies pop one cube at a time.  The symmetric policies
    group cubes whose centers are mirrored around ``domain.center`` across the
    coordinate planes.  Popping the whole group keeps refinement visually
    balanced when the target count stops the algorithm before the boundary is
    fully resolved.
    """

    def __init__(self, policy: str) -> None:
        self.policy = policy
        self._counter = count()
        self._symmetric = policy.startswith("symmetric_")
        self._priority = policy in {"priority", "symmetric_priority"}
        self._count = 0
        self._groups: dict[tuple[Any, ...], list[Cube]] = {}
        self._items: (
            deque[Cube]
            | deque[tuple[Any, ...]]
            | list[tuple[tuple[float, ...], int, Cube]]
            | list[tuple[tuple[float, ...], int, tuple[Any, ...]]]
        )
        self._items = [] if self._priority else deque()

    def __bool__(self) -> bool:
        return self._count > 0

    def __len__(self) -> int:
        return self._count

    def push(self, cube: Cube, domain: DomainOracle) -> None:
        self._count += 1
        if self._symmetric:
            key = _symmetry_orbit_key(cube, domain)
            self._groups.setdefault(key, []).append(cube)
            if self._priority:
                assert isinstance(self._items, list)
                heappush(self._items, (_queue_priority(cube, domain), next(self._counter), key))
            else:
                assert isinstance(self._items, deque)
                self._items.append(key)
            return

        if not self._priority:
            assert isinstance(self._items, deque)
            self._items.append(cube)
            return
        assert isinstance(self._items, list)
        heappush(self._items, (_queue_priority(cube, domain), next(self._counter), cube))

    def pop_batch(self) -> list[Cube]:
        """Return one refinement transaction, possibly a symmetric orbit."""
        if self._symmetric:
            return self._pop_symmetric_batch()
        cube = self._pop_single()
        self._count -= 1
        return [cube]

    def _pop_single(self) -> Cube:
        if not self._priority:
            assert isinstance(self._items, deque)
            return self._items.popleft()
        assert isinstance(self._items, list)
        return heappop(self._items)[2]

    def _pop_symmetric_batch(self) -> list[Cube]:
        while self._items:
            if self._priority:
                assert isinstance(self._items, list)
                key = heappop(self._items)[2]
            else:
                assert isinstance(self._items, deque)
                key = self._items.popleft()
            batch = self._groups.pop(key, [])
            if batch:
                self._count -= len(batch)
                return sorted(batch, key=lambda cell: (cell.level, cell.i, cell.j, cell.k))
        raise IndexError("pop from an empty boundary queue")


def mesh_from_cubes(
    cubes: Sequence[Cube],
    stats: CubeFillStats,
    domain: DomainOracle,
    bbox: Sequence[Sequence[float]] | None = None,
    *,
    target_tiles: int | None = None,
) -> PrismMesh:
    """Convert accepted cubes to the legacy ``PrismMesh`` container."""
    centers = np.array([cube.center for cube in cubes], dtype=float).reshape(-1, 3)
    dimensions = np.array(
        [[cube.h, cube.h, cube.h] for cube in cubes],
        dtype=float,
    ).reshape(-1, 3)
    levels = np.array([cube.level for cube in cubes], dtype=np.int64)
    indices = np.array([[cube.i, cube.j, cube.k] for cube in cubes], dtype=np.int64).reshape(-1, 3)
    root_bounds = _resolve_bbox(domain, bbox)
    metadata = {
        "cube_fill_stats": stats.as_dict(),
        "termination_reason": stats.termination_reason,
        "level_counts": _level_counts(levels),
        "lattice_indices": indices,
    }
    return PrismMesh(
        centers=centers,
        dimensions=dimensions,
        levels=levels,
        shape=getattr(domain, "name", domain.__class__.__name__),
        shape_metadata=_domain_metadata(domain),
        root_bounds=root_bounds,
        target_tiles=target_tiles,
        refinement_metadata=metadata,
        diagnostics=(
            f"adaptive cube fill accepted {stats.accepted_count} cubes",
            f"refinement stopped: {stats.termination_reason}",
        ),
    )


def fill_domain_mesh(
    domain: DomainOracle,
    *,
    bbox: Sequence[Sequence[float]] | None = None,
    N_target: int,
    h0: float | None = None,
    h_min: float = 0.0,
    max_depth: int = 12,
    eps: float = 1e-12,
    overshoot_policy: str = "closest",
    queue_policy: str = "breadth_first",
    eta: float = 0.75,
    complete_layers: bool = True,
    target_lower_tolerance: float = 0.10,
    grid_shifts: Sequence[Sequence[float]] | str | None = None,
    min_inside_fraction: float = 1.0,
    inside_fraction_samples: int = 3,
) -> PrismMesh:
    """Return a ``PrismMesh`` generated from any ``DomainOracle``.

    This is the generic bridge from the conservative cube-filling API to the
    MagTense-style mesh container.  Use it when you already have a domain
    object, or when experimenting with a custom oracle.  The ``generate_*``
    helpers below are convenience wrappers that build a domain object for a
    common analytic shape and then call this function.
    """
    cubes, stats = adaptive_cube_fill(
        domain=domain,
        bbox=bbox,
        N_target=N_target,
        h0=h0,
        h_min=h_min,
        max_depth=max_depth,
        eps=eps,
        overshoot_policy=overshoot_policy,
        queue_policy=queue_policy,
        eta=eta,
        complete_layers=complete_layers,
        target_lower_tolerance=target_lower_tolerance,
        grid_shifts=grid_shifts,
        min_inside_fraction=min_inside_fraction,
        inside_fraction_samples=inside_fraction_samples,
    )
    return mesh_from_cubes(cubes, stats, domain, bbox, target_tiles=N_target)


def estimate_h0_from_volume(
    domain: DomainOracle,
    N_target: int,
    eta: float = 0.75,
) -> float:
    """Estimate coarse cube side length from analytic volume and target count.

    This is the same volume scaling used when ``adaptive_cube_fill`` receives
    ``h0=None``.  ``eta`` is a dimensionless fill-efficiency factor: values
    below 1 assume only part of the domain volume will be represented by
    accepted cubes at the requested target count.  The estimate is intentionally
    simple; calibration applies dimensionless multipliers around it to find a
    practical starting resolution for a specific shape and target range.
    """
    _validate_target(N_target)
    eta_value = _validate_positive_float(eta, "eta")
    volume = getattr(domain, "volume", None)
    if volume is None:
        raise ValueError("domain must provide volume to estimate h0")
    return _validate_positive_float(
        (eta_value * float(volume) / N_target) ** (1.0 / 3.0),
        "estimated h0",
    )


def calibrate_h0_sweep(
    domain: DomainOracle,
    target_values: Sequence[int],
    h0_multipliers: Sequence[float],
    *,
    bbox: Sequence[Sequence[float]] | None = None,
    h_min_ratio: float = 1.0 / 16.0,
    max_depth: int = 8,
    eps: float = 1e-12,
    overshoot_policy: str = "soft",
    queue_policy: str = "symmetric_priority",
    eta: float = 0.75,
    target_count_min_ratio: float | None = None,
    target_count_max_ratio: float | None = None,
    target_lower_tolerance: float = 0.10,
    complete_layers: bool = True,
    grid_shifts: Sequence[Sequence[float]] | str | None = None,
    min_inside_fraction: float = 1.0,
    inside_fraction_samples: int = 3,
) -> list[dict[str, Any]]:
    """Benchmark several ``h0`` multipliers for one domain.

    Each multiplier is applied to ``estimate_h0_from_volume(domain, target)``.
    ``eps`` is the conservative classification tolerance passed through to the
    fill algorithm.  ``eta`` controls only the analytic reference scale; the
    sweep multipliers then move above or below that scale.  Fill fraction is
    only comparable between rows with similar accepted cell counts, so every
    row records whether it falls inside the target-count band.  The returned
    records are plain dictionaries so notebooks can display them directly or
    convert them to a dataframe.
    """
    if len(target_values) == 0:
        raise ValueError("target_values must not be empty")
    if len(h0_multipliers) == 0:
        raise ValueError("h0_multipliers must not be empty")
    h_min_ratio_value = _validate_nonnegative_float(h_min_ratio, "h_min_ratio")
    tolerance_value = _validate_target_lower_tolerance(target_lower_tolerance)
    min_ratio = (
        1.0 - tolerance_value
        if target_count_min_ratio is None
        else _validate_nonnegative_float(target_count_min_ratio, "target_count_min_ratio")
    )
    max_ratio = (
        float("inf")
        if target_count_max_ratio is None
        else _validate_nonnegative_float(target_count_max_ratio, "target_count_max_ratio")
    )
    if max_ratio < min_ratio:
        raise ValueError(
            "target_count_max_ratio must be greater than or equal to target_count_min_ratio"
        )

    records: list[dict[str, Any]] = []
    shape_key = _shape_key(domain)
    for target in target_values:
        _validate_target(target)
        analytic_h0 = estimate_h0_from_volume(domain, int(target), eta=eta)
        for multiplier in h0_multipliers:
            multiplier_value = _validate_positive_float(multiplier, "h0_multiplier")
            h0_value = analytic_h0 * multiplier_value
            preflight = _level0_calibration_preflight(
                domain,
                bbox,
                h0_value,
                eps,
            )
            if preflight["accepted_count"] > max_ratio * int(target):
                count_ratio = preflight["accepted_count"] / int(target)
                records.append(
                    {
                        "shape": shape_key,
                        "target": int(target),
                        "analytic_h0": analytic_h0,
                        "h0_multiplier": multiplier_value,
                        "h0": h0_value,
                        "accepted_count": preflight["accepted_count"],
                        "count_error": abs(int(target) - preflight["accepted_count"]),
                        "count_ratio": count_ratio,
                        "within_target_band": False,
                        "target_status": "above_target",
                        "target_accepted": True,
                        "relative_count_error": count_ratio - 1.0,
                        "target_count_min_ratio": min_ratio,
                        "target_count_max_ratio": max_ratio,
                        "fill_fraction": None,
                        "level_min": 0 if preflight["accepted_count"] else None,
                        "level_max": 0 if preflight["accepted_count"] else None,
                        "max_level_reached": 0,
                        "unresolved_boundary_cells": preflight["boundary_count"],
                        "classified_cells": preflight["classified_cells"],
                        "rejected_outside_cells": preflight["rejected_outside_cells"],
                        "termination_reason": "level-0 count exceeds target band",
                        "h_min": h0_value * h_min_ratio_value,
                        "max_depth": max_depth,
                        "overshoot_policy": overshoot_policy,
                        "queue_policy": queue_policy,
                        "short_circuited": True,
                        "selection_warning": "increase h0 or increase target cells",
                    }
                )
                continue

            cubes, stats = adaptive_cube_fill(
                domain=domain,
                bbox=bbox,
                N_target=int(target),
                h0=h0_value,
                h_min=h0_value * h_min_ratio_value,
                max_depth=max_depth,
                eps=eps,
                overshoot_policy=overshoot_policy,
                queue_policy=queue_policy,
                eta=eta,
                complete_layers=complete_layers,
                target_lower_tolerance=tolerance_value,
                grid_shifts=grid_shifts,
                min_inside_fraction=min_inside_fraction,
                inside_fraction_samples=inside_fraction_samples,
            )
            levels = [cube.level for cube in cubes]
            count_ratio = stats.accepted_count / int(target)
            final_h_min = min((cube.h for cube in cubes), default=None)
            records.append(
                {
                    "shape": shape_key,
                    "target": int(target),
                    "target_cells": int(target),
                    "analytic_h0": analytic_h0,
                    "h0_multiplier": multiplier_value,
                    "h0": h0_value,
                    "accepted_count": stats.accepted_count,
                    "final_cells": stats.accepted_count,
                    "count_error": abs(int(target) - stats.accepted_count),
                    "count_ratio": count_ratio,
                    "relative_count_error": stats.relative_count_error,
                    "within_target_band": min_ratio <= count_ratio <= max_ratio,
                    "target_status": stats.target_status,
                    "target_accepted": stats.target_accepted,
                    "target_count_min_ratio": min_ratio,
                    "target_count_max_ratio": max_ratio,
                    "fill_fraction": stats.fill_fraction,
                    "level_min": min(levels) if levels else None,
                    "level_max": max(levels) if levels else None,
                    "layers_used": stats.layers_used,
                    "final_smallest_tile_size": final_h_min,
                    "max_level_reached": stats.max_level_reached,
                    "unresolved_boundary_cells": stats.boundary_cells_left_unresolved,
                    "classified_cells": stats.classified_cells,
                    "rejected_outside_cells": stats.rejected_outside_cells,
                    "termination_reason": stats.termination_reason,
                    "h_min": stats.h_min,
                    "max_depth": stats.max_depth,
                    "overshoot_policy": stats.overshoot_policy,
                    "queue_policy": stats.queue_policy,
                    "complete_layers": stats.complete_layers,
                    "target_lower_tolerance": stats.target_lower_tolerance,
                    "selected_grid_shift": stats.selected_grid_shift,
                    "short_circuited": False,
                    "selection_warning": "",
                }
            )
    return records


def select_best_h0_records(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Return the best calibration row for each ``(shape, target)`` pair.

    Rows satisfying the target-count acceptance rule are preferred.  Within
    that valid set, ranking chooses the closest cell count first, then uses
    fill fraction as a tie-breaker.  If no row is within the accepted band, the fallback
    picks the closest accepted-cell count and marks the row with a warning.
    """
    grouped: dict[tuple[str, int], list[Mapping[str, Any]]] = {}
    for record in records:
        key = (str(record["shape"]), int(record["target"]))
        grouped.setdefault(key, []).append(record)

    best_records: list[dict[str, Any]] = []
    for key in sorted(grouped, key=lambda item: (item[0], item[1])):
        group = grouped[key]
        valid = [
            record
            for record in group
            if bool(record.get("target_accepted", record.get("within_target_band", True)))
            and not bool(record.get("short_circuited", False))
        ]
        if valid:
            selected = dict(min(valid, key=_h0_record_sort_key))
            selected["selection_warning"] = ""
        else:
            selected = dict(min(group, key=_h0_fallback_sort_key))
            selected.setdefault("selection_warning", "")
            if not selected["selection_warning"]:
                selected["selection_warning"] = "no rows within target-count band"
        best_records.append(selected)
    return best_records


def recommended_h0(
    domain: DomainOracle,
    N_target: int,
    *,
    shape_key: str | None = None,
    eta: float = 0.75,
) -> float:
    """Return a reusable lookup-table ``h0`` recommendation.

    If ``shape_key`` or ``domain.name`` matches ``H0_REFERENCE_TABLE``, the
    matching row's multiplier is applied to the analytic volume estimate.  For
    unknown/custom domains, the function falls back to the analytic estimate.
    """
    _validate_target(N_target)
    analytic_h0 = estimate_h0_from_volume(domain, N_target, eta=eta)
    resolved_shape = _shape_key(domain if shape_key is None else shape_key)
    for record in H0_REFERENCE_TABLE.get(resolved_shape, ()):
        if record.target_min <= N_target <= record.target_max:
            return analytic_h0 * record.h0_multiplier
    return analytic_h0


def generate_rectangular_prism(
    dimensions: Sequence[float],
    nx: int,
    ny: int,
    nz: int,
    *,
    center: Sequence[float] = (0.0, 0.0, 0.0),
) -> PrismMesh:
    """Generate a regular ``nx`` by ``ny`` by ``nz`` rectangular-prism mesh.

    This helper is intentionally uniform and non-adaptive.  It is kept for
    compatibility with examples that need a simple box grid rather than a
    conservative boundary-fitted fill.
    """
    dimensions_array = _vector3(dimensions, "dimensions", positive=True)
    center_array = _vector3(center, "center")
    if not all(_is_integer(value) and value > 0 for value in (nx, ny, nz)):
        raise ValueError("nx, ny, and nz must be positive integers")
    resolution = np.asarray([nx, ny, nz], dtype=int)
    cell_dimensions = dimensions_array / resolution
    axes = [
        center_array[index]
        + (np.arange(resolution[index]) + 0.5) * cell_dimensions[index]
        - dimensions_array[index] / 2.0
        for index in range(3)
    ]
    grid = np.meshgrid(*axes, indexing="ij")
    centers = np.column_stack([axis.ravel() for axis in grid])
    cell_count = int(np.prod(resolution))
    bounds = np.vstack((center_array - dimensions_array / 2, center_array + dimensions_array / 2))
    return PrismMesh(
        centers=centers,
        dimensions=np.tile(cell_dimensions, (cell_count, 1)),
        levels=np.zeros(cell_count, dtype=np.int64),
        shape="rectangular_prism",
        shape_metadata={
            "center": center_array.copy(),
            "dimensions": dimensions_array.copy(),
            "resolution": tuple(int(value) for value in resolution),
        },
        root_bounds=bounds,
    )


def generate_sphere(
    radius: float,
    target_tiles: int,
    *,
    center: Sequence[float] = (0.0, 0.0, 0.0),
    h0: float | None = None,
    h_min: float = 0.0,
    max_depth: int = 12,
    eps: float = 1e-12,
    overshoot_policy: str = "closest",
    queue_policy: str = "breadth_first",
    complete_layers: bool = True,
    target_lower_tolerance: float = 0.10,
    grid_shifts: Sequence[Sequence[float]] | str | None = None,
    min_inside_fraction: float = 1.0,
    inside_fraction_samples: int = 3,
) -> PrismMesh:
    """Generate an adaptive conservative cube mesh for a sphere.

    Convenience wrapper: constructs ``Sphere(center, radius)`` and forwards the
    fill parameters to ``fill_domain_mesh``.
    """
    domain = Sphere(center=_vector3(center, "center"), radius=radius)
    return fill_domain_mesh(
        domain,
        N_target=target_tiles,
        h0=h0,
        h_min=h_min,
        max_depth=max_depth,
        eps=eps,
        overshoot_policy=overshoot_policy,
        queue_policy=queue_policy,
        complete_layers=complete_layers,
        target_lower_tolerance=target_lower_tolerance,
        grid_shifts=grid_shifts,
        min_inside_fraction=min_inside_fraction,
        inside_fraction_samples=inside_fraction_samples,
    )


def generate_ellipsoid(
    semi_axes: Sequence[float],
    target_tiles: int,
    *,
    center: Sequence[float] = (0.0, 0.0, 0.0),
    h0: float | None = None,
    h_min: float = 0.0,
    max_depth: int = 12,
    eps: float = 1e-12,
    overshoot_policy: str = "closest",
    queue_policy: str = "breadth_first",
    complete_layers: bool = True,
    target_lower_tolerance: float = 0.10,
    grid_shifts: Sequence[Sequence[float]] | str | None = None,
    min_inside_fraction: float = 1.0,
    inside_fraction_samples: int = 3,
) -> PrismMesh:
    """Generate an adaptive conservative cube mesh for an axis-aligned ellipsoid.

    Convenience wrapper: constructs ``Ellipsoid(center, semi_axes)`` and
    forwards the fill parameters to ``fill_domain_mesh``.
    """
    domain = Ellipsoid(
        center=_vector3(center, "center"),
        semi_axes=_vector3(semi_axes, "semi_axes", positive=True),
    )
    return fill_domain_mesh(
        domain,
        N_target=target_tiles,
        h0=h0,
        h_min=h_min,
        max_depth=max_depth,
        eps=eps,
        overshoot_policy=overshoot_policy,
        queue_policy=queue_policy,
        complete_layers=complete_layers,
        target_lower_tolerance=target_lower_tolerance,
        grid_shifts=grid_shifts,
        min_inside_fraction=min_inside_fraction,
        inside_fraction_samples=inside_fraction_samples,
    )


def generate_cylinder(
    radius: float,
    length: float,
    axis: str,
    target_tiles: int,
    *,
    center: Sequence[float] = (0.0, 0.0, 0.0),
    h0: float | None = None,
    h_min: float = 0.0,
    max_depth: int = 12,
    eps: float = 1e-12,
    overshoot_policy: str = "closest",
    queue_policy: str = "breadth_first",
    complete_layers: bool = True,
    target_lower_tolerance: float = 0.10,
    grid_shifts: Sequence[Sequence[float]] | str | None = None,
    min_inside_fraction: float = 1.0,
    inside_fraction_samples: int = 3,
) -> PrismMesh:
    """Generate an adaptive conservative cube mesh for a finite cylinder.

    Convenience wrapper: constructs ``Cylinder(center, radius, length, axis)``
    and forwards the fill parameters to ``fill_domain_mesh``.
    """
    domain = Cylinder(
        center=_vector3(center, "center"),
        radius=radius,
        length=length,
        axis=axis,
    )
    return fill_domain_mesh(
        domain,
        N_target=target_tiles,
        h0=h0,
        h_min=h_min,
        max_depth=max_depth,
        eps=eps,
        overshoot_policy=overshoot_policy,
        queue_policy=queue_policy,
        complete_layers=complete_layers,
        target_lower_tolerance=target_lower_tolerance,
        grid_shifts=grid_shifts,
        min_inside_fraction=min_inside_fraction,
        inside_fraction_samples=inside_fraction_samples,
    )


def generate_hexagonal_prism(
    side_length: float,
    height: float,
    axis: str,
    target_tiles: int,
    *,
    center: Sequence[float] = (0.0, 0.0, 0.0),
    rotation_degrees: float = 0.0,
    h0: float | None = None,
    h_min: float = 0.0,
    max_depth: int = 12,
    eps: float = 1e-12,
    overshoot_policy: str = "closest",
    queue_policy: str = "breadth_first",
    complete_layers: bool = True,
    target_lower_tolerance: float = 0.10,
    grid_shifts: Sequence[Sequence[float]] | str | None = None,
    min_inside_fraction: float = 1.0,
    inside_fraction_samples: int = 3,
) -> PrismMesh:
    """Generate an adaptive conservative cube mesh for a regular hexagonal prism.

    Convenience wrapper: constructs ``HexagonalPrism`` and forwards the fill
    parameters to ``fill_domain_mesh``.
    """
    domain = HexagonalPrism(
        center=_vector3(center, "center"),
        side_length=side_length,
        height=height,
        axis=axis,
        rotation_degrees=rotation_degrees,
    )
    return fill_domain_mesh(
        domain,
        N_target=target_tiles,
        h0=h0,
        h_min=h_min,
        max_depth=max_depth,
        eps=eps,
        overshoot_policy=overshoot_policy,
        queue_policy=queue_policy,
        complete_layers=complete_layers,
        target_lower_tolerance=target_lower_tolerance,
        grid_shifts=grid_shifts,
        min_inside_fraction=min_inside_fraction,
        inside_fraction_samples=inside_fraction_samples,
    )


def _coarse_grid(
    bbox: FloatArray,
    h0: float,
    *,
    grid_shift: Sequence[float] | None = None,
) -> list[Cube]:
    lower, upper = bbox
    shift = np.zeros(3, dtype=float) if grid_shift is None else _vector3(grid_shift, "grid_shift")
    origin = lower + shift
    first = np.floor((lower - origin) / h0).astype(int)
    last = np.ceil((upper - origin) / h0).astype(int)
    counts = np.maximum(1, last - first)
    cubes: list[Cube] = []
    ranges = [range(int(start), int(start + count)) for start, count in zip(first, counts)]
    for i, j, k in product(*ranges):
        lattice_index = np.array([i, j, k], dtype=float)
        center = origin + (lattice_index + 0.5) * h0
        cubes.append(Cube(center=center, h=h0, level=0, i=i, j=j, k=k))
    return cubes


def _subdivide(cube: Cube) -> tuple[Cube, ...]:
    child_h = cube.h / 2.0
    offset = cube.h / 4.0
    children: list[Cube] = []
    for dx, dy, dz in product((0, 1), repeat=3):
        signs = np.array([2 * dx - 1, 2 * dy - 1, 2 * dz - 1], dtype=float)
        children.append(
            Cube(
                center=cube.center + signs * offset,
                h=child_h,
                level=cube.level + 1,
                i=2 * cube.i + dx,
                j=2 * cube.j + dy,
                k=2 * cube.k + dz,
            )
        )
    return tuple(children)


def _symmetry_orbit_key(cube: Cube, domain: DomainOracle) -> tuple[Any, ...]:
    """Return a key shared by cubes mirrored through ``domain.center``.

    The key preserves coordinate order, so it groups mirror images across the
    x/y/z center planes but does not impose rotational symmetry.  If a custom
    domain has no usable center, the cube is its own orbit and the symmetric
    queue gracefully behaves like the corresponding non-symmetric policy.
    """
    center = getattr(domain, "center", None)
    if center is None:
        return ("single", cube.level, cube.i, cube.j, cube.k)
    try:
        domain_center = _vector3(center, "center")
    except ValueError:
        return ("single", cube.level, cube.i, cube.j, cube.k)

    relative = np.abs(cube.center - domain_center)
    # Quantize by the cube size to avoid tiny floating-point differences
    # splitting mirrored cells into separate groups.
    scaled = tuple(int(round(value / cube.h * 1_000_000_000)) for value in relative)
    return ("mirror", cube.level, scaled)


def _queue_priority(cube: Cube, domain: DomainOracle) -> tuple[float, ...]:
    """Priority tuple for boundary cells; lower values are refined first."""
    boundary_distance = abs(float(domain.sdf(cube.center)))
    return (
        -cube.h,
        _off_axis_priority(cube, domain),
        boundary_distance,
        float(cube.level),
    )


def _off_axis_priority(cube: Cube, domain: DomainOracle) -> float:
    """Prefer radial/off-axis boundary cells for axis-aligned prism-like shapes.

    For cylinders and hexagonal prisms, axis-aligned cubes fit naturally along
    the cylinder/prism axis and struggle most around the radial boundary.
    Returning ``radial_margin - axial_margin`` makes radial-boundary cells sort
    before cap-boundary cells.  Other domains use a neutral score.
    """
    if not isinstance(domain, (Cylinder, HexagonalPrism)):
        return 0.0

    relative = cube.center - domain.center
    axis_index = _axis_index(domain.axis)
    axial_extent = domain.length / 2.0 if isinstance(domain, Cylinder) else domain.height / 2.0
    axial_margin = axial_extent - abs(float(relative[axis_index]))
    radial_indices = [index for index in range(3) if index != axis_index]
    radial = relative[radial_indices]

    if isinstance(domain, Cylinder):
        radial_margin = domain.radius - float(np.linalg.norm(radial))
    else:
        radial_margin = domain._apothem - float(np.max(domain._face_normals @ radial))
    return radial_margin - axial_margin


def _classify_candidate_cube(
    domain: DomainOracle,
    cube: Cube,
    eps: float,
    min_inside_fraction: float,
    inside_fraction_samples: int,
) -> CubeState:
    """Classify a candidate cube, using strict or sampled partial acceptance."""
    if min_inside_fraction >= 1.0:
        return domain.classify_cube(cube, eps)

    fraction = _sampled_inside_fraction(domain, cube, eps, inside_fraction_samples)
    if fraction >= min_inside_fraction:
        return CubeState.INSIDE
    if fraction <= 0.0:
        return CubeState.OUTSIDE
    return CubeState.BOUNDARY


def _sampled_inside_fraction(
    domain: DomainOracle,
    cube: Cube,
    eps: float,
    samples_per_axis: int,
) -> float:
    """Approximate the fraction of a cube inside a generic SDF domain."""
    offsets = np.linspace(
        -cube.h / 2.0,
        cube.h / 2.0,
        samples_per_axis + 2,
        dtype=float,
    )[1:-1]
    inside = 0
    total = samples_per_axis**3
    for dx, dy, dz in product(offsets, repeat=3):
        point = cube.center + np.array([dx, dy, dz], dtype=float)
        if domain.sdf(point) >= -eps:
            inside += 1
    return inside / total


def _resolve_grid_shifts(
    grid_shifts: Sequence[Sequence[float]] | str | None,
    h0: float,
) -> tuple[FloatArray, ...]:
    """Return absolute grid-shift candidates for a resolved ``h0``."""
    if grid_shifts is None:
        return (np.zeros(3, dtype=float),)
    if isinstance(grid_shifts, str):
        if grid_shifts != "half_step":
            raise ValueError("grid_shifts must be None, 'half_step', or explicit 3-vectors")
        values = (-0.5 * h0, 0.0, 0.5 * h0)
        return tuple(np.array(shift, dtype=float) for shift in product(values, repeat=3))

    resolved = tuple(_vector3(shift, "grid_shift") for shift in grid_shifts)
    if len(resolved) == 0:
        raise ValueError("grid_shifts must contain at least one shift")
    return resolved


def _fill_run_sort_key(stats: CubeFillStats) -> tuple[int, float, float]:
    """Rank completed fill candidates by target validity, count error, fill."""
    fill_fraction = -float(stats.fill_fraction) if stats.fill_fraction is not None else 0.0
    relative_error = stats.relative_count_error
    error = abs(float(relative_error)) if relative_error is not None else float("inf")
    return (0 if stats.target_accepted else 1, error, fill_fraction)


def _grid_shift_summary(stats: CubeFillStats) -> Mapping[str, Any]:
    """Compact metadata for one grid-shift candidate."""
    return {
        "grid_shift": stats.selected_grid_shift,
        "accepted_count": stats.accepted_count,
        "fill_fraction": stats.fill_fraction,
        "target_accepted": stats.target_accepted,
        "target_status": stats.target_status,
        "relative_count_error": stats.relative_count_error,
        "layers_used": stats.layers_used,
        "termination_reason": stats.termination_reason,
    }


def _copy_stats_with_grid_shift_candidates(
    stats: CubeFillStats,
    candidates: tuple[Mapping[str, Any], ...],
) -> CubeFillStats:
    """Return ``stats`` with run-level grid-shift candidate summaries attached."""
    values = stats.as_dict()
    values["selected_grid_shift"] = tuple(stats.selected_grid_shift)
    values["grid_shift_candidates"] = candidates
    return CubeFillStats(**values)


def _target_accepted(
    accepted_count: int,
    target: int,
    target_lower_tolerance: float,
) -> bool:
    lower = (1.0 - target_lower_tolerance) * target
    upper = (1.0 + target_lower_tolerance) * target
    return lower <= accepted_count <= upper


def _target_status(
    accepted_count: int,
    target: int,
    target_lower_tolerance: float,
) -> str:
    if accepted_count > (1.0 + target_lower_tolerance) * target:
        return "above_tolerance"
    if accepted_count == target:
        return "at_target"
    if accepted_count > target:
        return "within_upper_tolerance"
    if _target_accepted(accepted_count, target, target_lower_tolerance):
        return "within_lower_tolerance"
    return "below_tolerance"


def _shape_key(domain_or_name: DomainOracle | str) -> str:
    """Return the normalized key used by the reusable h0 reference table."""
    if isinstance(domain_or_name, str):
        name = domain_or_name
    else:
        name = getattr(domain_or_name, "name", domain_or_name.__class__.__name__)
    return str(name).lower().replace(" ", "_")


def _level0_calibration_preflight(
    domain: DomainOracle,
    bbox: Sequence[Sequence[float]] | None,
    h0: float,
    eps: float,
) -> dict[str, int]:
    """Classify only the coarse grid for calibration short-circuiting."""
    bbox_array = _resolve_bbox(domain, bbox)
    accepted = 0
    boundary = 0
    outside = 0
    classified = 0
    for cube in _coarse_grid(bbox_array, h0):
        state = domain.classify_cube(cube, eps)
        classified += 1
        if state == CubeState.INSIDE:
            accepted += 1
        elif state == CubeState.BOUNDARY:
            boundary += 1
        else:
            outside += 1
    return {
        "accepted_count": accepted,
        "boundary_count": boundary,
        "rejected_outside_cells": outside,
        "classified_cells": classified,
    }


def _h0_record_sort_key(record: Mapping[str, Any]) -> tuple[int, float, int, int]:
    """Sort key where lower is better according to calibration rules."""
    fill_fraction = record.get("fill_fraction")
    if fill_fraction is None:
        fill_score = float("inf")
    else:
        fill_score = -float(fill_fraction)
    return (
        int(record.get("count_error", 0)),
        fill_score,
        int(record.get("unresolved_boundary_cells", 0)),
        int(record.get("classified_cells", 0)),
    )


def _h0_fallback_sort_key(record: Mapping[str, Any]) -> tuple[int, float, int, int]:
    """Fallback sort key used when every row misses the target-count band."""
    fill_fraction = record.get("fill_fraction")
    fill_score = float("inf") if fill_fraction is None else -float(fill_fraction)
    return (
        int(record.get("count_error", 0)),
        fill_score,
        int(record.get("unresolved_boundary_cells", 0)),
        int(record.get("classified_cells", 0)),
    )


def _resolve_bbox(
    domain: DomainOracle,
    bbox: Sequence[Sequence[float]] | None,
) -> FloatArray:
    source = bbox if bbox is not None else getattr(domain, "bounding_box", None)
    if source is None:
        raise ValueError("bbox is required when the domain does not provide bounding_box")
    bbox_array = np.asarray(source, dtype=float)
    if bbox_array.shape != (2, 3) or not np.all(np.isfinite(bbox_array)):
        raise ValueError("bbox must have shape (2, 3) with finite values")
    if np.any(bbox_array[1] <= bbox_array[0]):
        raise ValueError("bbox upper bounds must be greater than lower bounds")
    return bbox_array


def _resolve_h0(
    domain: DomainOracle,
    bbox: FloatArray,
    N_target: int,
    h0: float | None,
    eta: float,
) -> float:
    if h0 is not None:
        return _validate_positive_float(h0, "h0")
    estimate = estimate_h0_from_volume(domain, N_target, eta=eta)
    span_limit = float(np.max(bbox[1] - bbox[0]))
    return min(_validate_positive_float(estimate, "estimated h0"), span_limit)


def _domain_metadata(domain: DomainOracle) -> Mapping[str, Any]:
    metadata: dict[str, Any] = {}
    for name in (
        "center",
        "radius",
        "length",
        "axis",
        "dimensions",
        "semi_axes",
        "side_length",
        "height",
        "rotation_degrees",
        "volume",
    ):
        if hasattr(domain, name):
            value = getattr(domain, name)
            metadata[name] = value.copy() if isinstance(value, np.ndarray) else value
    if hasattr(domain, "_vertices"):
        metadata["vertices"] = getattr(domain, "_vertices").copy()
    if hasattr(domain, "_apothem"):
        metadata["apothem"] = getattr(domain, "_apothem")
    return metadata


def _level_counts(levels: IntArray) -> dict[int, int]:
    if len(levels) == 0:
        return {}
    unique, counts_array = np.unique(levels, return_counts=True)
    return {int(level): int(count_value) for level, count_value in zip(unique, counts_array)}


def _vector3(value: Sequence[float], name: str, *, positive: bool = False) -> FloatArray:
    result = np.asarray(value, dtype=float)
    if result.shape != (3,) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain three finite values")
    if positive and np.any(result <= 0):
        raise ValueError(f"all {name} values must be positive")
    return result


def _axis_index(axis: str) -> int:
    if axis not in {"x", "y", "z"}:
        raise ValueError("axis must be 'x', 'y', or 'z'")
    return {"x": 0, "y": 1, "z": 2}[axis]


def _validate_target(target_tiles: int) -> None:
    if not _is_integer(target_tiles) or target_tiles <= 0:
        raise ValueError("N_target must be a positive integer")


def _validate_positive_float(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return result


def _validate_nonnegative_float(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be non-negative and finite")
    return result


def _validate_max_depth(value: int) -> int:
    if not _is_integer(value) or value < 0:
        raise ValueError("max_depth must be a non-negative integer")
    return int(value)


def _validate_target_lower_tolerance(value: float) -> float:
    result = _validate_nonnegative_float(value, "target_lower_tolerance")
    if result >= 1.0:
        raise ValueError("target_lower_tolerance must be less than 1.0")
    return result


def _validate_min_inside_fraction(value: float) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0.5 or result > 1.0:
        raise ValueError("min_inside_fraction must satisfy 0.5 <= value <= 1.0")
    return result


def _validate_sample_count(value: int) -> int:
    if not _is_integer(value) or value <= 0:
        raise ValueError("inside_fraction_samples must be a positive integer")
    return int(value)


def _is_integer(value: Any) -> bool:
    return (
        not isinstance(value, (bool, np.bool_))
        and isinstance(value, (int, float, np.integer, np.floating))
        and np.isfinite(value)
        and float(value).is_integer()
    )


def conservative_cube_condition(domain: DomainOracle, cube: Cube, eps: float) -> bool:
    """Return whether ``cube`` satisfies the default conservative SDF test."""
    radius = sqrt(3.0) * cube.h / 2.0
    return bool(domain.sdf(cube.center) >= radius + eps)


def cubes_overlap(a: Cube, b: Cube, tol: float = 1e-12) -> bool:
    """Return True when two axis-aligned cubes overlap with positive volume."""
    a_min = a.center - a.h / 2.0
    a_max = a.center + a.h / 2.0
    b_min = b.center - b.h / 2.0
    b_max = b.center + b.h / 2.0
    return bool(np.all(np.minimum(a_max, b_max) - np.maximum(a_min, b_min) > tol))


def assert_no_cube_overlaps(cubes: Sequence[Cube], tol: float = 1e-12) -> None:
    """Raise ``AssertionError`` if any cubes overlap with positive volume."""
    for first in range(len(cubes)):
        for second in range(first + 1, len(cubes)):
            if cubes_overlap(cubes[first], cubes[second], tol=tol):
                raise AssertionError(f"cubes {first} and {second} overlap")


__all__ = [
    "AxisAlignedBox",
    "Cube",
    "CubeFillStats",
    "CubeState",
    "Cylinder",
    "DomainOracle",
    "Ellipsoid",
    "H0CalibrationRecord",
    "H0_REFERENCE_TABLE",
    "HexagonalPrism",
    "PrismMesh",
    "SignedDistanceDomain",
    "Sphere",
    "adaptive_cube_fill",
    "assert_no_cube_overlaps",
    "calibrate_h0_sweep",
    "conservative_cube_condition",
    "cubes_overlap",
    "estimate_h0_from_volume",
    "fill_domain_mesh",
    "generate_cylinder",
    "generate_ellipsoid",
    "generate_hexagonal_prism",
    "generate_rectangular_prism",
    "generate_sphere",
    "mesh_from_cubes",
    "recommended_h0",
    "select_best_h0_records",
]
