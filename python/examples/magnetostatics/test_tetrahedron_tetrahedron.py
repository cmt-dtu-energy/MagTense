"""
Regression tests for the analytically exact tetrahedron-source to
tetrahedron-target volume-averaged demagnetization tensor.

The tensor under test satisfies ``<H>_{V_target} = N @ M``: it maps the uniform
magnetization of the source tetrahedron to the demagnetizing field averaged
over the volume of the receiving tetrahedron. It is evaluated in closed form
from sixteen analytical triangle-pair Laplace integrals over the outward faces
of the two tetrahedra, with the MagTense normalization ``-1 / (4 pi V_target)``
(no division by the source volume).

Three quantities are compared for each gap:

* ``N_exact``     the new analytical operator,
* ``N_reference`` the supplied analytical reference tensors,
* ``N_gauss``     an independent order-10 tensor-product Gauss-Legendre Duffy
                  average (exactly 1000 points) of the EXISTING MagTense
                  tetrahedron-source to point-target tensor over the target.
"""

import numpy as np
import pytest

from magtense.magstatics import (
    Tiles,
    get_demag_tensor,
    get_demag_tensor_tetrahedron_pair,
)

GAPS = (0.0, 1.0, 5.0, 10.0)

Z_BASE = 0.866025
Y_APEX = 0.816497
Z_APEX = 0.288675

TET0_TOUCHING = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.0, Z_BASE],
        [0.5, Y_APEX, Z_APEX],
    ],
    dtype=np.float64,
)

TET1_TOUCHING = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.0, Z_BASE],
        [0.5, -Y_APEX, Z_APEX],
    ],
    dtype=np.float64,
)


def absolute_vertices(g):
    shift = np.array([0.0, 0.5 * g, 0.0], dtype=np.float64)
    source = TET0_TOUCHING + shift
    target = TET1_TOUCHING - shift
    return source, target


ANALYTICAL = {
    0.0: np.array(
        [
            [-7.52869185174101507e-02, 3.55204657140909045e-17, 2.15614738879441870e-16],
            [3.55204657140909045e-17, 1.50573868067920058e-01, 4.56557709928796125e-18],
            [2.15614738879441870e-16, 4.56557709928796125e-18, -7.52869495505097547e-02],
        ],
        dtype=np.float64,
    ),
    1.0: np.array(
        [
            [-3.19651279863232898e-03, 3.36923256459795206e-17, 6.51077535607359464e-17],
            [3.36923256459795206e-17, 6.39302584138304059e-03, 2.13955342056489863e-17],
            [6.51077535607359464e-17, 2.13955342056489863e-17, -3.19651304275001512e-03],
        ],
        dtype=np.float64,
    ),
    5.0: np.array(
        [
            [-5.92247968139480690e-05, 2.73375812732718110e-15, 3.64363301862866220e-15],
            [2.73375812732718110e-15, 1.18449593929324820e-04, -5.06857723432721162e-16],
            [3.64363301862866220e-15, -5.06857723432721162e-16, -5.92247971133411273e-05],
        ],
        dtype=np.float64,
    ),
    10.0: np.array(
        [
            [-8.31623222274667256e-06, -5.35327175535581492e-15, 2.26408301037721316e-14],
            [-5.35327175535581492e-15, 1.66324643230237244e-05, 1.05882810235410917e-14],
            [2.26408301037721316e-14, 1.05882810235410917e-14, -8.31623221730665258e-06],
        ],
        dtype=np.float64,
    ),
}

# Tolerance on ||N_exact - N_reference|| / ||N_reference||, per gap.
#
# These are not free parameters. Two effects set the achievable agreement and
# both grow with the separation:
#
#  1. The face-pair representation is a difference of sixteen surface integrals
#     whose individual magnitude is set by the face areas over the separation,
#     while the result decays like a dipole field. At g = 10 the sum cancels by
#     about three decades, which costs about three digits in double precision.
#     Recomputing the same routine in quadruple precision changes the g = 10
#     result by ~6e-9 relative and leaves g = 0 unchanged, which is exactly this
#     effect and not an implementation error.
#  2. The supplied reference carries its own error of the same origin. At g = 5
#     and g = 10 its off-diagonal entries are ~3e-15 and ~2e-14 where symmetry
#     demands exactly zero, and at g = 10 its trace departs from zero by 7e-10
#     relative although the exact tensor is trace free for disjoint bodies. The
#     independent 1000-point quadrature below reproduces the deviation of this
#     implementation from the reference at g = 5 and g = 10 to within a few
#     percent, confirming the reference is what differs there.
#
# Measured with ifx 2026.1 at -O3: 2.1e-15, 1.3e-13, 1.5e-10, 3.2e-9. The values
# below keep roughly a factor 20-30 of headroom for other compilers and
# optimisation settings.
REL_TOL_EXACT = {0.0: 5e-14, 1.0: 4e-12, 5.0: 4e-9, 10.0: 8e-8}

# Tolerance on ||N_gauss - N_reference|| / ||N_reference||, per gap. The g = 0
# entry is large by construction: the two tetrahedra share a face, so the
# point-target integrand is singular on part of the target boundary and a smooth
# product rule converges slowly there. That is the whole point of having the
# analytical operator.
REL_TOL_GAUSS = {0.0: 5e-4, 1.0: 1e-11, 5.0: 1e-9, 10.0: 1e-7}

GAUSS_ORDER = 10


def _gauss_duffy_target_points(vertices, order=GAUSS_ORDER):
    """
    Order-``order`` tensor-product Gauss-Legendre Duffy rule on a tetrahedron.

    Returns the ``order**3`` evaluation points and the corresponding weights,
    normalized to sum to one so that the weighted sum is a volume AVERAGE.
    """
    node, weight = np.polynomial.legendre.leggauss(order)
    node = 0.5 * (node + 1.0)
    weight = 0.5 * weight

    r, s, t = np.meshgrid(node, node, node, indexing="ij")
    wr, ws, wt = np.meshgrid(weight, weight, weight, indexing="ij")
    r, s, t = r.ravel(), s.ravel(), t.ravel()

    # dV / V = 6 r^2 s dr ds dt for this Duffy map.
    w = 6.0 * r**2 * s * wr.ravel() * ws.ravel() * wt.ravel()

    lambda0 = 1.0 - r
    lambda1 = r * (1.0 - s)
    lambda2 = r * s * (1.0 - t)
    lambda3 = r * s * t

    pts = (
        lambda0[:, None] * vertices[0]
        + lambda1[:, None] * vertices[1]
        + lambda2[:, None] * vertices[2]
        + lambda3[:, None] * vertices[3]
    )

    # The exact weights already sum to one; renormalize for round-off only.
    return pts, w / w.sum()


def _tetrahedron_point_tensor(source_vertices, pts):
    """The EXISTING MagTense tetrahedron-source to point-target tensor."""
    tile = Tiles(
        n=1,
        tile_type=5,
        offset=[0, 0, 0],
        M_rem=0.0,
        mu_r_ea=1.0,
        mu_r_oa=1.0,
    )
    tile.vertices = source_vertices
    return get_demag_tensor(tile, np.asfortranarray(pts))[0]


def _gauss_target_average(source_vertices, target_vertices):
    pts, w = _gauss_duffy_target_points(target_vertices)
    assert pts.shape == (GAUSS_ORDER**3, 3)
    return np.einsum("q,qij->ij", w, _tetrahedron_point_tensor(source_vertices, pts))


def _relative_error(value, reference):
    return np.linalg.norm(value - reference) / np.linalg.norm(reference)


@pytest.fixture(scope="module")
def results():
    out = {}
    for g in GAPS:
        source, target = absolute_vertices(g)
        out[g] = {
            "exact": get_demag_tensor_tetrahedron_pair(source, target),
            "reference": ANALYTICAL[g],
            "gauss": _gauss_target_average(source, target),
        }
    return out


@pytest.mark.parametrize("g", GAPS)
def test_exact_matches_analytical_reference(results, g):
    """A. The analytical operator reproduces the supplied reference tensors."""
    r = results[g]
    err = _relative_error(r["exact"], r["reference"])
    assert err < REL_TOL_EXACT[g], f"g={g}: relative error {err:.3e}"


@pytest.mark.parametrize("g", GAPS)
def test_gauss_1000_target_average_is_converged(results, g):
    """B. The independent 1000-point target average is converged."""
    r = results[g]
    err = _relative_error(r["gauss"], r["reference"])
    assert err < REL_TOL_GAUSS[g], f"g={g}: relative error {err:.3e}"


def test_analytical_beats_quadrature_for_touching_geometry(results):
    """
    C. Where the comparison is numerically meaningful, the analytical operator
    is closer to the reference than the 1000-point quadrature.

    The touching geometry is the meaningful case: the quadrature has to resolve
    a boundary singularity and is about eleven orders of magnitude worse there.
    For the separated gaps both results sit at or below the accuracy of the
    reference itself, so a strict ordering between them carries no information;
    the assertion there only requires that the analytical operator is not
    materially worse.
    """
    touching = results[0.0]
    err_exact = _relative_error(touching["exact"], touching["reference"])
    err_gauss = _relative_error(touching["gauss"], touching["reference"])
    assert err_exact < err_gauss / 1.0e6, (
        f"touching: exact {err_exact:.3e} vs quadrature {err_gauss:.3e}"
    )

    for g in (1.0, 5.0, 10.0):
        r = results[g]
        err_exact = _relative_error(r["exact"], r["reference"])
        err_gauss = _relative_error(r["gauss"], r["reference"])
        assert err_exact <= max(10.0 * err_gauss, 1.0e-8), (
            f"g={g}: exact {err_exact:.3e} vs quadrature {err_gauss:.3e}"
        )


@pytest.mark.parametrize("g", GAPS)
def test_tensor_structure(results, g):
    """D. Structural invariants of the tensor."""
    N = results[g]["exact"]

    assert N.shape == (3, 3)
    assert np.all(np.isfinite(N)), f"g={g}: non-finite entries"

    scale = np.linalg.norm(N)
    assert scale > 0.0

    # Symmetric to numerical precision.
    assert np.abs(N - N.T).max() <= 1e-14 * scale

    # The two tetrahedra are mirror images about the y = 0 plane and are
    # separated along y, so the mixed components vanish by symmetry.
    off_diagonal = max(abs(N[0, 1]), abs(N[0, 2]), abs(N[1, 2]))
    assert off_diagonal <= REL_TOL_EXACT[g] * scale

    # The bodies are disjoint (they touch at most on a face of zero volume), so
    # the tensor of the external field is trace free.
    assert abs(np.trace(N)) <= REL_TOL_EXACT[g] * scale


def test_touching_faces_are_finite():
    """The g = 0 geometry shares a whole face and must stay finite."""
    source, target = absolute_vertices(0.0)
    N = get_demag_tensor_tetrahedron_pair(source, target)
    assert np.all(np.isfinite(N))
    assert not np.any(np.isnan(N))


@pytest.mark.parametrize("g", GAPS)
def test_vertex_permutation_and_orientation_invariance(g):
    """
    The face normals are made outward explicitly, so neither the order of the
    vertices nor the handedness of the supplied ordering may change the result.
    """
    source, target = absolute_vertices(g)
    reference = get_demag_tensor_tetrahedron_pair(source, target)
    scale = np.linalg.norm(reference)

    # Permuting the vertices permutes the faces and therefore the order in
    # which the sixteen face-pair integrals are summed, so the agreement is
    # limited by the same cancellation budget as the comparison against the
    # reference tensors and is held to the same per-gap tolerance.
    #
    # An odd permutation flips the sign of the signed volume, an even one does
    # not; both have to give the same tensor.
    for perm in ([1, 0, 2, 3], [0, 2, 1, 3], [3, 2, 1, 0], [1, 2, 3, 0]):
        permuted = get_demag_tensor_tetrahedron_pair(source[perm], target[perm])
        assert np.abs(permuted - reference).max() <= REL_TOL_EXACT[g] * scale


def test_translation_invariance():
    """A common rigid translation of both tetrahedra may not change the tensor."""
    source, target = absolute_vertices(1.0)
    reference = get_demag_tensor_tetrahedron_pair(source, target)
    scale = np.linalg.norm(reference)

    shift = np.array([123.25, -47.5, 8.75], dtype=np.float64)
    shifted = get_demag_tensor_tetrahedron_pair(source + shift, target + shift)
    assert np.abs(shifted - reference).max() <= 1e-11 * scale


def _volume(vertices):
    a = vertices[1] - vertices[0]
    b = vertices[2] - vertices[0]
    c = vertices[3] - vertices[0]
    return abs(np.dot(a, np.cross(b, c))) / 6.0


def test_reciprocity_between_unequal_tetrahedra():
    """
    Reciprocity of the double surface integral gives

        V_target * N(source -> target) == V_source * N(target -> source)

    (both tensors being symmetric). Deliberately uses two tetrahedra of
    different size and shape, for which the relation is not trivially the
    symmetry of the geometry.
    """
    source = np.array(
        [[0.0, 0.0, 0.0], [1.3, 0.0, 0.0], [0.2, 0.9, 0.1], [0.4, 0.1, 1.7]],
        dtype=np.float64,
    )
    target = np.array(
        [[2.5, 0.4, 0.3], [3.1, 0.1, -0.2], [2.7, 1.4, 0.0], [2.6, 0.5, 0.8]],
        dtype=np.float64,
    )

    v_source = _volume(source)
    v_target = _volume(target)
    assert not np.isclose(v_source, v_target)

    forward = v_target * get_demag_tensor_tetrahedron_pair(source, target)
    backward = v_source * get_demag_tensor_tetrahedron_pair(target, source)

    # This pair is separated by about twice its own size, so the face-pair sum
    # already cancels by a couple of decades; the tolerance is that budget and
    # not a free parameter (measured 6.4e-11 with ifx 2026.1 at -O3).
    assert np.abs(forward - backward).max() <= 1e-9 * np.linalg.norm(forward)


def test_far_separation_agrees_with_point_target():
    """
    At large separation the target-averaged tensor and the existing
    point-target tensor at the target centroid have to agree, since the source
    field is then nearly constant across the target.
    """
    source, target = absolute_vertices(40.0)
    averaged = get_demag_tensor_tetrahedron_pair(source, target)

    centroid = target.mean(axis=0)[None, :]
    at_centroid = _tetrahedron_point_tensor(source, centroid)[0]

    assert np.abs(averaged - at_centroid).max() <= 1e-3 * np.linalg.norm(averaged)


def test_self_tensor_is_trace_minus_one():
    """
    The coincident case: source and target are the same tetrahedron. The
    volume-averaged self-demagnetization tensor of any body has trace exactly
    -1, which is an independent check of the coincident-face branch that the
    tetrahedral micromagnetic path relies on for its diagonal element.
    """
    for vertices in (
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 0.866025403784439, 0.0],
                [0.5, 0.288675134594813, 0.816496580927726],
            ]
        ),
        np.array([[0.0, 0.0, 0.0], [1.3, 0.0, 0.0], [0.2, 0.9, 0.1], [0.4, 0.1, 1.7]]),
    ):
        N = get_demag_tensor_tetrahedron_pair(vertices, vertices)
        assert np.all(np.isfinite(N))
        assert np.abs(N - N.T).max() <= 1e-14 * np.linalg.norm(N)
        assert abs(np.trace(N) + 1.0) < 1e-12


def test_regular_tetrahedron_self_tensor_is_isotropic():
    """A regular tetrahedron is isotropic enough that its averaged self-demag
    tensor is -I/3, the same as a sphere."""
    vertices = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.866025403784439, 0.0],
            [0.5, 0.288675134594813, 0.816496580927726],
        ]
    )
    N = get_demag_tensor_tetrahedron_pair(vertices, vertices)
    assert np.abs(N + np.eye(3) / 3.0).max() < 1e-12


def test_face_sharing_neighbours_are_finite_and_trace_free():
    """Two tetrahedra that share a whole face, the closest non-overlapping
    configuration a tetrahedral mesh produces."""
    source = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    target = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]])
    N = get_demag_tensor_tetrahedron_pair(source, target)
    assert np.all(np.isfinite(N))
    assert np.abs(N - N.T).max() <= 1e-14 * np.linalg.norm(N)
    assert abs(np.trace(N)) <= 1e-12 * np.linalg.norm(N)


@pytest.mark.parametrize("dz", [2.0, 0.5, 1.0e-3, 1.0e-6, 0.0])
def test_parallel_faces_branch(dz):
    """
    Two congruent tetrahedra stacked along z, so a face of one is parallel to a
    face of the other. The reduction takes its dedicated parallel-plane branch
    here, and it has to stay finite and trace free right down to the two faces
    touching (dz = 0).
    """
    source = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.3, 0.3, 1.0]])
    target = source.copy()
    target[:, 2] += 1.0 + dz

    N = get_demag_tensor_tetrahedron_pair(source, target)
    scale = np.linalg.norm(N)

    assert np.all(np.isfinite(N))
    assert np.abs(N - N.T).max() <= 1e-14 * scale
    assert abs(np.trace(N)) <= 1e-11 * scale


@pytest.mark.parametrize(
    "target",
    [
        # shares only the single vertex (0, 0, 1)
        np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 2.0], [0.0, 1.0, 2.0], [0.0, 0.0, 2.0]]),
        # shares only the edge from (0, 0, 0) to (1, 0, 0)
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]),
    ],
    ids=["vertex-touching", "edge-sharing"],
)
def test_vertex_and_edge_touching_are_finite(target):
    """Lower-dimensional contact, which a tetrahedral mesh produces constantly."""
    source = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    N = get_demag_tensor_tetrahedron_pair(source, target)
    scale = np.linalg.norm(N)
    assert np.all(np.isfinite(N))
    assert abs(np.trace(N)) <= 1e-11 * scale


def test_degenerate_tetrahedron_is_rejected():
    """A flat (zero volume) tetrahedron has to be reported, not silently used."""
    flat = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        dtype=np.float64,
    )
    _, target = absolute_vertices(1.0)

    with pytest.raises(ValueError):
        get_demag_tensor_tetrahedron_pair(flat, target)
