# %%
import os
import sys
import h5py
import numpy as np

from pathlib import Path
from typing import Optional, List
from multiprocessing import Process, cpu_count
from magtense.halbach import HalbachCylinder, EvaluationPoints
from magtense.magstatics import Tiles, grid_config, run_simulation

from tqdm import tqdm

from db_utils import load_demag_tensor, calc_demag_field, eval_shimming
from db_utils import get_shim_mat, gen_s_state, gen_seq


def db_std_prob_4(
    datapath: Path,
    n_seq: int,
    res: List,
    grid_size: List = [500e-9, 500e-9, 3e-9],
    t_steps: int = 500,
    t_per_step: float = 4e-12,
    h_ext_a: List = [0, 360],
    h_ext_n: List = [0, 50],
    seed: int = 0,
    intv: Optional[List] = None,
    name: Optional[str] = None,
    empty: bool = False,
    cuda: bool = False,
) -> None:
    fname = name if name is not None else f"{n_seq}_{t_steps}_{res[0]}_{res[1]}"
    if intv is None:
        intv = [0, n_seq]
    if not empty:
        fname += f"_{intv[0]}_{intv[1]}"
    n_intv = intv[1] - intv[0]

    db = h5py.File(f"{datapath}/{fname}.h5", "w")
    db.create_dataset(
        "sequence", shape=(n_intv, t_steps, 3, res[0], res[1]), dtype="float32"
    )
    db.create_dataset("field", shape=(n_intv, 3), dtype="float32")
    if not empty:
        db.attrs["intv"] = intv

    if empty:
        s_state = gen_s_state(res, grid_size)
        np.save(f"{datapath}/{res[0]}_{res[1]}_s_state.npy", s_state)
        db.attrs["res"] = res
        db.attrs["grid_size"] = grid_size
        db.attrs["t_steps"] = t_steps
        db.attrs["t_per_step"] = t_per_step
        db.attrs["h_ext_angle"] = h_ext_a
        db.attrs["h_ext_norm"] = h_ext_n
        db.attrs["seed"] = seed
        db.close()
        return fname, n_seq

    rng = np.random.default_rng(seed)
    rnd_mat = rng.random(size=(n_seq, 2))

    for i in tqdm(range(n_intv)):
        h_ext = np.zeros(3)
        d = (h_ext_n[1] - h_ext_n[0]) * rnd_mat[i, 0] + h_ext_n[0]
        theta = np.deg2rad((h_ext_a[1] - h_ext_a[0]) * rnd_mat[i, 1] + h_ext_a[0])
        h_ext[:2] = d * np.array([np.cos(theta), np.sin(theta)])
        db["field"][i] = h_ext

        seq = gen_seq(
            m0_state=np.load(f"{datapath}/{res[0]}_{res[1]}_s_state.npy"),
            res=res,
            grid_size=grid_size,
            h_ext=h_ext,
            t_steps=t_steps,
            t_per_step=t_per_step,
            cuda=cuda,
        )

        # Output shape: (t, comp, res_x, res_y)
        db["sequence"][i] = np.moveaxis(
            seq.reshape(t_steps, res[1], res[0], 3).swapaxes(1, 2), -1, 1
        )

    db.close()


def db_halbach(
    datapath: Path,
    n_halbach: int,
    n_mat: int,
    shim_segs: int = 8,
    shim_layers: int = 3,
    field_res: List[int] = [16, 16, 16],
    norm_var: float = 0.05,
    mu_r: float = 100,
    r_probe: float = 0.002,
    symm: bool = True,
    seed_shim: int = 0,
    seed_pertubation: int = 1,
    action: bool = False,
    no_shim: bool = True,
    intv: Optional[List] = None,
    name: Optional[str] = None,
    empty: bool = False,
):
    """
    Dataset creation.
    Name: '{available spots for shim magnets}_{# halbach configs}_{# shim matrices}.h5'

    Args:
        datapath: Path to database.
        n_halbach: Number of different Halbach configuration.
        n_mat: Number of shim matrices per Halbach configuration.
        shim_segs: Number of shim magnets per layer.
        shim_layers: Layers of shim magnets in z-direction.
        field_res: Resolution of measured field in cylindrical bore.
        norm_var: Percentage pertubation of remanent magnetization.
        mu_r: Relative permeability.
        r_probe: Defines size of cube, where the magnetic field is measured.
        symm: If set, then symmetry along xy-plane is assumed.
        seed_shim: Seed for random number generator of shim matrix.
        seed_pertubation: Seed for random number generator of pertubation.
        action: If set, action as single placement of a shim magnet is stored.
        no_shim: If set, information of unshimmed field is stored.
        intv: Data range to iterate over.
        name: Optional filename for database.
        empty: If set, an empty database is created.
    """
    halbach = HalbachCylinder()
    halbach.add_shim_magnets(shim_layers, shim_segs, mu_r=mu_r)
    pts_eval = EvaluationPoints(field_res, r_probe)
    load_demag_tensor(halbach, pts_eval, store=True)
    load_demag_tensor(halbach, halbach.pts_shim, store=True)

    rng_shim = np.random.default_rng(seed_shim)
    spots = shim_segs * (shim_layers // 2 + shim_layers % 2) if symm else halbach.n_shim
    symm_mat = rng_shim.choice([0, 1], size=(n_halbach, n_mat, spots))

    fname = name if name is not None else f"{spots}_{n_halbach}_{n_mat}_{field_res[0]}"
    if intv is None:
        intv = [0, n_halbach]
    n_config = intv[1] - intv[0]
    if not empty:
        fname += f"_{intv[0]}_{intv[1]}"

    db = h5py.File(f"{datapath}/{fname}.h5", libver="latest", mode="w")
    db.create_dataset(
        "halbach_M_rem", shape=(n_config, halbach.n_hard_tiles), dtype="float32"
    )
    db.create_dataset(
        "halbach_ea", shape=(n_config, halbach.n_hard_tiles, 3), dtype="float32"
    )
    db.create_dataset(
        "field",
        shape=(n_config, n_mat, 3, field_res[0], field_res[1], field_res[2]),
        dtype="float32",
    )
    db.create_dataset(
        "shim_field",
        shape=(n_config, n_mat, halbach.pts_shim.coords.shape[0], 3),
        dtype="float32",
    )
    db.create_dataset("p2p", shape=(n_config, n_mat, 1), dtype="float64")
    if not empty:
        db.attrs["intv"] = intv

    if no_shim:
        db.create_dataset(
            "field_no_shim",
            shape=(n_config, 3, field_res[0], field_res[1], field_res[2]),
            dtype="float32",
        )
        db.create_dataset(
            "shim_field_no_shim",
            shape=(n_config, halbach.pts_shim.coords.shape[0], 3),
            dtype="float32",
        )
        db.create_dataset("p2p_no_shim", shape=(n_config, 1), dtype="float64")

    if action:
        db.create_dataset("p2p_next", shape=(n_config, n_mat, 1), dtype="float64")
        rng_action = np.random.default_rng(seed_shim)
        empty_spots = np.where(symm_mat == 0)
        np_action = np.zeros(n_halbach, n_mat)

        for i in range(n_halbach):
            for j in range(n_mat):
                if (i in empty_spots[0]) and (j in empty_spots[1]):
                    np_action[i, j] = rng_action.integers(
                        np.where(symm_mat[i, j] == 0)[0].shape[0]
                    )
                else:
                    np_action[i, j] = -1000

    if empty:
        db.create_dataset("shim_mat", shape=(n_halbach, n_mat, spots), dtype="int32")
        db["shim_mat"][:] = symm_mat

        if action:
            db.create_dataset("action", shape=(n_halbach, n_mat, 1), dtype="int32")
            db["action"][:] = np_action

        db.attrs["n_halbach"] = n_halbach
        db.attrs["n_mat"] = n_mat
        db.attrs["shim_pos"] = halbach.pts_shim.coords
        db.attrs["shim_segs"] = shim_segs
        db.attrs["shim_layers"] = shim_layers
        db.attrs["mu_r"] = mu_r
        db.attrs["r_probe"] = r_probe
        db.attrs["norm_var"] = norm_var
        db.attrs["symmetry"] = symm
        db.attrs["field_res"] = field_res
        db.attrs["seed_shim"] = seed_shim
        db.attrs["seed_pertubation"] = seed_pertubation
        db.close()

        return fname, n_halbach

    rng_per = np.random.default_rng(seed_pertubation)
    B_rem = rng_per.normal(
        loc=halbach.B_rem,
        scale=norm_var * halbach.B_rem,
        size=(n_halbach, halbach.n_hard_tiles),
    )

    with tqdm(total=n_config * n_mat) as pbar:
        for idx_per in range(n_config):
            halbach.perturb_config(B_rem[intv[0] + idx_per])
            db["halbach_M_rem"][idx_per] = halbach.remanence
            db["halbach_ea"][idx_per] = halbach.easy_axes

            if no_shim:
                shim_field_no_shim = calc_demag_field(halbach, halbach.pts_shim)
                B_field_no_shim, p2p_no_shim = eval_shimming(halbach, pts_eval)

                db["shim_field_no_shim"][idx_per] = shim_field_no_shim
                db["field_no_shim"][idx_per] = B_field_no_shim
                db["p2p_no_shim"][idx_per] = p2p_no_shim

            for idx_shim in range(n_mat):
                shim_mat = get_shim_mat(
                    symm_mat[intv[0] + idx_per, idx_shim], shim_layers, shim_segs, symm
                )
                halbach.set_shim_matrix(shim_mat)
                shim_field = calc_demag_field(halbach, halbach.pts_shim)
                B_field, p2p = eval_shimming(halbach, pts_eval)

                db["shim_field"][idx_per, idx_shim] = shim_field
                db["field"][idx_per, idx_shim] = B_field
                db["p2p"][idx_per, idx_shim] = p2p

                if action:
                    if np_action[intv[0] + idx_per, idx_shim] < 0:
                        p2p_next = np_action[intv[0] + idx_per, idx_shim]

                    else:
                        a = np_action[intv[0] + idx_per, idx_shim]
                        halbach.set_one_shim_magnet(a)
                        if symm:
                            a_symm = (
                                halbach.n_shim
                                - ((a // halbach.shim_segs + 1) * halbach.shim_segs)
                                + (a % halbach.shim_segs)
                            )
                            halbach.set_one_shim_magnet(a_symm)
                        _, p2p_next = eval_shimming(halbach, pts_eval)

                    db["p2p_next"][idx_per, idx_shim] = p2p_next

                pbar.update(1)

    db.close()


def db_halbach_field_stats(db_path):
    with h5py.File(db_path, mode="r") as db:
        arr = np.mean(db["field"], axis=(2, 3, 4))

    print(f"Mean: {np.mean(arr, axis=0)}")
    print(f"Std: {np.std(arr, axis=0)}")


def db_magfield(
    datapath: Path,
    n_samples: int,
    res: List[int],
    spots: List = [10, 10, 5],
    area: List = [1, 1, 0.5],
    gap: float = 0.05,
    seed: int = 0,
    intv: Optional[List] = None,
    name: str = "magfield",
    empty: bool = False,
) -> None:
    """
    Generate 3-D magnetic fields of experimental setup.

    Args:
        filepath: Indicates where to store the data sample.
        n_samples: Size of database.
        res: Resolution of magnetic field.
        spots: Available positions in setup.
        area: Size of setup.
        gap: Gap between measurement area and surrounding magnets.
        seed: Seed for random number generator of matrices.
        intv: Data range to iterate over.
        empty: If set, an empty database is created.
    """
    fname = f"{name}_{res[0]}"
    if intv is None:
        intv = [0, n_samples]
    if not empty:
        fname += f"_{intv[0]}_{intv[1]}"
    n_intv = intv[1] - intv[0]

    db = h5py.File(f"{datapath}/{fname}.h5", libver="latest", mode="w")
    out_shape = (n_intv, 3, *res)
    db.create_dataset("field", shape=out_shape, dtype="float32")
    if not empty:
        db.attrs["intv"] = intv

    if empty:
        db.attrs["spots"] = spots
        db.attrs["area"] = area
        db.attrs["gap"] = gap
        db.attrs["seed"] = seed
        db.close()

        return fname, n_samples

    rng = np.random.default_rng(seed)
    tile_size = np.asarray(area) / np.asarray(spots)
    filled_mat = rng.integers(2, size=(n_samples, spots[0], spots[1], spots[2]))
    empty_mat = rng.integers(4, size=(n_samples,))

    for idx in tqdm(range(n_intv)):
        emp_x, emp_y = {0: [4, 5], 1: [4, 6], 2: [3, 6], 3: [3, 7]}[
            empty_mat[idx + intv[0]]
        ]
        s_x = emp_x * tile_size[0]
        s_y = emp_y * tile_size[1]

        filled_pos = [
            [i, j, k]
            for i in range(spots[0])
            for j in range(spots[1])
            for k in range(spots[2])
            if filled_mat[intv[0] + idx][i][j][k] == 1
            and (
                i < emp_x
                or i > emp_y
                or j < emp_x
                or j > emp_y
                or k < emp_x
                or k > emp_y
            )
        ]

        tiles, _ = grid_config(spots, area, filled_pos)

        x_eval = np.linspace(s_x + gap, s_y + gap, res[0])
        y_eval = np.linspace(s_x + gap, s_y + gap, res[1])

        if len(res) == 2:
            xv, yv = np.meshgrid(x_eval, y_eval)
            zv = np.zeros(res[0] * res[1]) + area[2] / 2

        elif len(res) == 3:
            z_eval = np.linspace(s_x + gap, s_y + gap, res[2])
            # # Pixel length in z-direction equal to x-direction
            # s_z = (s_y - s_x) / res[0]
            # z_eval = np.linspace(-s_z, s_z, res[2]) + area[2] / 2
            xv, yv, zv = np.meshgrid(x_eval, y_eval, z_eval)

        else:
            raise ValueError("Only 2-D and 3-D magnetic field can be generated!")

        pts_eval = np.hstack([xv.reshape(-1, 1), yv.reshape(-1, 1), zv.reshape(-1, 1)])
        devnull = open("/dev/null", "w")
        oldstdout_fno = os.dup(sys.stdout.fileno())
        os.dup2(devnull.fileno(), 1)
        _, h_out = run_simulation(tiles, pts_eval)
        os.dup2(oldstdout_fno, 1)

        # Tensor image with shape CxHxWxD [T]
        field = h_out.reshape((*res, 3))
        if len(res) == 3:
            field = field.transpose((3, 0, 1, 2))
        else:
            field = field.transpose((2, 0, 1))

        db["field"][idx] = field * 4 * np.pi * 1e-7

    db.close()


def db_magfield_symm(
    datapath: Path,
    n_samples: int,
    res: List[int],
    spots: List = [10, 10, 5],
    area: List = [1, 1, 0.5],
    gap: float = 0.05,
    seed: int = 0,
    intv: Optional[List] = None,
    name: str = "magfield_symm",
    empty: bool = False,
) -> None:
    """
    Generate 3-D magnetic fields of experimental setup.

    Args:
        filepath: Indicates where to store the data sample.
        n_samples: Size of database.
        res: Resolution of magnetic field.
        spots: Available positions in setup.
        area: Size of setup.
        gap: Gap between measurement area and surrounding magnets.
        seed: Seed for random number generator of matrices.
        intv: Data range to iterate over.
        empty: If set, an empty database is created.
    """
    symm = True
    fname = f"{name}_{res[0]}"
    if intv is None:
        intv = [0, n_samples]
    if not empty:
        fname += f"_{intv[0]}_{intv[1]}"
    n_intv = intv[1] - intv[0]

    db = h5py.File(f"{datapath}/{fname}.h5", libver="latest", mode="w")
    out_shape = (n_intv, 3, *res) if not symm else (n_intv, 2, *res)
    db.create_dataset("field", shape=out_shape, dtype="float32")
    if not empty:
        db.attrs["intv"] = intv

    if symm:
        spots[2] = 1
        area[2] = 10

    if empty:
        db.attrs["spots"] = spots
        db.attrs["area"] = area
        db.attrs["gap"] = gap
        db.attrs["seed"] = seed
        db.close()

        return fname, n_samples

    rng = np.random.default_rng(seed)
    tile_size = np.asarray(area) / np.asarray(spots)

    filled_mat = rng.integers(2, size=(n_samples, spots[0], spots[1], spots[2]))
    mag_angle_mat = rng.random(size=(n_samples, spots[0] * spots[1] * spots[2]))

    for idx in tqdm(range(n_intv)):
        emp_x = spots[0] // 2 - 1
        emp_y = emp_x + 1
        s_x = emp_x * tile_size[0]
        s_y = emp_y * tile_size[1]

        filled_pos = [
            [i, j, k]
            for i in range(spots[0])
            for j in range(spots[1])
            for k in range(spots[2])
            if filled_mat[intv[0] + idx][i][j][k] == 1
            and (
                i < emp_x
                or i > emp_y
                or j < emp_x
                or j > emp_y
                or k < (spots[2] // 2 - 1)
                or k > spots[2] // 2
            )
        ]

        tiles = Tiles(
            n=len(filled_pos),
            size=tile_size,
            tile_type=2,
            M_rem=1.2 / (4 * np.pi * 1e-7),
            color=[1, 0, 0],
            mag_angle=[
                [np.pi / 2, 2 * np.pi * mag_angle_mat[idx + intv[0], i]]
                for i in range(len(filled_pos))
            ],
        )

        for i, pos in enumerate(filled_pos):
            pos = np.asarray(pos)
            if np.greater_equal(pos, spots).any():
                raise ValueError(f"Desired position {pos} is not in the grid!")
            tiles.offset = (np.around((pos + 0.5) * tile_size, decimals=9), i)

        x_eval = np.linspace(s_x + gap, s_y + tile_size[0] - gap, res[0])
        y_eval = np.linspace(s_x + gap, s_y + tile_size[1] - gap, res[1])

        if len(res) == 2:
            xv, yv = np.meshgrid(x_eval, y_eval)
            zv = np.zeros(res[0] * res[1]) + area[2] / 2

        elif len(res) == 3:
            # Pixel length in z-direction equal to x-direction
            s_z = (s_y - s_x) / res[0]
            z_eval = np.linspace(-s_z, s_z, res[2]) + area[2] / 2
            xv, yv, zv = np.meshgrid(x_eval, y_eval, z_eval)

        else:
            raise ValueError("Only 2-D and 3-D magnetic field can be generated!")

        pts_eval = np.hstack([xv.reshape(-1, 1), yv.reshape(-1, 1), zv.reshape(-1, 1)])
        devnull = open("/dev/null", "w")
        oldstdout_fno = os.dup(sys.stdout.fileno())
        os.dup2(devnull.fileno(), 1)
        _, h_out = run_simulation(tiles, pts_eval)
        os.dup2(oldstdout_fno, 1)

        # Tensor image with shape CxHxWxD [T]
        field = h_out.reshape((*res, 3))
        field = (
            field.transpose((3, 0, 1, 2))
            if len(res) == 3
            else field.transpose((2, 0, 1))
        )
        if symm:
            field = field[:2]

        db["field"][idx] = field * 4 * np.pi * 1e-7

    db.close()


def db_single_magnets(
    datapath: Path,
    n_samples: int,
    res: int,
    num_mag: int,
    dim: int = 2,
    M_fixed: bool = False,
    min_pos: float = 0,
    min_magnet_size: float = 0.5,
    max_magnet_size: float = 2,
    height: float = 1,
    area_x: float = 10,
    seed: int = 0,
    intv: Optional[List] = None,
    name: str = "mags",
    empty: bool = False,
) -> None:
    fname = f"{name}_{num_mag}_{res}"
    if M_fixed:
        datapath += "_M_fixed"
    if intv is None:
        intv = [0, n_samples]
    if not empty:
        fname += f"_{intv[0]}_{intv[1]}"
    n_intv = intv[1] - intv[0]

    if empty:
        if Path(f"{datapath}/{fname}.h5").is_file():
            input(f'Overwriting file "{fname}.h5". Press Enter to continue')

    n_params = 6 if dim == 3 else 5
    if M_fixed:
        n_params -= 1

    db = h5py.File(f"{datapath}/{fname}.h5", libver="latest", mode="w")
    db.create_dataset("field", shape=(n_intv, dim, res, res), dtype="float32")
    db.create_dataset("sample", shape=(n_intv, num_mag, n_params), dtype="float32")
    db.create_dataset("info", shape=(n_intv,), dtype="float32")
    if not empty:
        db.attrs["intv"] = intv

    if empty:
        db.attrs["spots"] = M_fixed
        db.attrs["min_pos"] = min_pos
        db.attrs["min_magnet_size"] = min_magnet_size
        db.attrs["max_magnet_size"] = max_magnet_size
        db.attrs["height"] = height
        db.attrs["area_x"] = area_x
        db.attrs["seed"] = seed
        db.close()

        return fname, n_samples

    rng = np.random.default_rng(seed)
    param_mat = rng.random(size=(n_samples, num_mag, n_params))

    x_eval = np.linspace(0, area_x, res) + (area_x / (2 * res))
    xv, yv = np.meshgrid(x_eval, x_eval)
    zv = np.zeros(res * res) + height / 2
    pts_eval = np.hstack([xv.reshape(-1, 1), yv.reshape(-1, 1), zv.reshape(-1, 1)])

    for idx in tqdm(range(n_intv)):
        pos = np.zeros(shape=(num_mag, 3))
        size = np.zeros(shape=(num_mag, 3))
        M_rem = np.zeros(shape=(num_mag,))
        mag_angle = np.zeros(shape=(num_mag, 2))
        n_mag = num_mag

        for i in range(num_mag):
            pos[i] = [
                param_mat[intv[0] + idx, i, 0] * (area_x - 2 * min_pos) + min_pos,
                param_mat[intv[0] + idx, i, 1] * (area_x - 2 * min_pos) + min_pos,
                height / 2,
            ]
            size_x = (
                param_mat[intv[0] + idx, i, 2] * (max_magnet_size - min_magnet_size)
                + min_magnet_size
            )
            size[i] = [size_x, size_x, height]
            M_rem[i] = 1.28 if M_fixed else param_mat[intv[0] + idx, i, 3] + 0.5
            polar_angle = (
                param_mat[intv[0] + idx, i, -1] * np.pi if dim == 3 else np.pi / 2
            )
            mag_angle[i] = [param_mat[intv[0] + idx, i, -2] * 2 * np.pi, polar_angle]

        overlap = True
        while overlap:
            overlap = False
            for i in range(pos.shape[0] - 1):
                for j in range(i + 1, pos.shape[0]):
                    d = np.sqrt(((pos[i, :2] - pos[j, :2]) ** 2).sum())
                    if d < np.sqrt(2) * ((size[i][0] + size[j][0]) / 2):
                        overlap = True
                        pos = np.delete(pos, j, 0)
                        size = np.delete(size, j, 0)
                        M_rem = np.delete(M_rem, j, 0)
                        mag_angle = np.delete(mag_angle, j, 0)
                        n_mag -= 1

        tile = Tiles(
            n=pos.shape[0],
            tile_type=2,
            offset=pos,
            size=size,
            M_rem=M_rem / (4 * np.pi * 1e-7),
            mag_angle=mag_angle,
            color=[1, 0, 0],
        )
        devnull = open("/dev/null", "w")
        oldstdout_fno = os.dup(sys.stdout.fileno())
        os.dup2(devnull.fileno(), 1)
        _, h_out = run_simulation(tile, pts_eval)
        os.dup2(oldstdout_fno, 1)
        h_out = h_out.reshape((res, res, -1)).transpose((2, 0, 1))

        sample = np.zeros(shape=(num_mag, n_params), dtype=np.float32)
        for i in range(pos.shape[0]):
            sample[i, 0] = pos[i][0]
            sample[i, 1] = pos[i][1]
            sample[i, 2] = size[i][0]
            sample[i, 3] = M_rem[i]
            sample[i, 4] = mag_angle[i][0]
            if dim == 3:
                sample[i, 5] = mag_angle[i][1]

        db["field"][idx] = h_out[:dim] * 4 * np.pi * 1e-7
        db["sample"][idx] = sample
        db["info"][idx] = n_mag

    db.close()


def create_db_mp(
    data: str,
    n_workers: Optional[int] = None,
    datapath: Optional[Path] = None,
    **kwargs,
) -> None:
    if datapath is None:
        datapath = Path(__file__).parent.absolute() / ".." / "data"
    if not datapath.exists():
        datapath.mkdir(parents=True)
    kwargs["datapath"] = datapath

    if data == "halbach":
        target = db_halbach
    elif data == "std_prob_4":
        target = db_std_prob_4
    elif data == "magfield":
        target = db_magfield
    elif data == "magfield_symm":
        target = db_magfield_symm
    elif data == "single_magnets":
        target = db_single_magnets
    else:
        raise NotImplementedError()

    db_name, n_tasks = target(**kwargs, empty=True)

    if n_workers is None:
        n_workers = cpu_count()
    intv = n_tasks // n_workers
    if n_tasks % n_workers > 0:
        intv += 1

    l_p = []
    for i in range(n_workers):
        end_intv = min((i + 1) * intv, n_tasks)
        kwargs["intv"] = [i * intv, end_intv]
        p = Process(target=target, kwargs=kwargs)
        p.daemon = True
        p.start()
        l_p.append(p)
        if end_intv == n_tasks:
            break

    try:
        for p in l_p:
            p.join()

    except KeyboardInterrupt:
        for p in l_p:
            p.terminate()

        path = datapath.glob("**/*")
        fnames = [
            x.name
            for x in path
            if x.is_file()
            and x.name[: len(db_name)] == db_name
            and x.name[:-3] != db_name
        ]

        for name in fnames:
            Path(datapath, name).unlink()

        Path(datapath, f"{db_name}.h5").unlink()
        exit(130)

    path = datapath.glob("**/*")
    fnames = [
        x.name
        for x in path
        if x.is_file() and x.name[: len(db_name)] == db_name and x.name[:-3] != db_name
    ]

    with h5py.File(f"{datapath}/{db_name}.h5", mode="a") as db_t:
        for name in fnames:
            print(Path(datapath, name))
            with h5py.File(Path(datapath, name), mode="r") as db_s:
                intv = db_s.attrs["intv"]
                for key in db_s.keys():
                    db_t[key][intv[0] : intv[1]] = db_s[key]
            Path(datapath, name).unlink()

    print("Database created")


# %%

if __name__ == "__main__":
    # db_kwargs = {
    #     'n_halbach': 1000,
    #     'n_mat': 1,
    #     'shim_segs': 16,
    #     'shim_layers': 3,
    #     'name': 'test_pre_3000f_1r',
    #     'seed_shim' : 51,
    #     'seed_pertubation': 52
    # }

    # create_db_mp('halbach', n_workers=15, **db_kwargs)

    # db_kwargs = {
    #     'res': [16, 4, 1],
    #     'grid_size': [500e-9, 500e-9, 3e-9],
    #     'n_seq': 4,
    # }

    # create_db_mp('std_prob_4', n_workers=1, **db_kwargs)

    db_kwargs = {
        "n_samples": 30000,
        "res": [64, 64, 64],
        "spots": [10, 10, 10],
        "area": [1, 1, 1],
    }
    create_db_mp("magfield", n_workers=15, **db_kwargs)

    # db_kwargs = {
    #     'n_samples': 10,
    #     'res': 32,
    #     'num_mag': 3,
    #     'dim': 3
    # }
    # create_db_mp('single_magnets', n_workers=5, **db_kwargs)

# %%
