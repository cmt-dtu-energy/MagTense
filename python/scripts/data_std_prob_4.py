#%%
import h5py
import numpy as np
import matplotlib.pyplot as plt

from magtense.micromag import MicromagProblem
from typing import Optional, List
from pathlib import Path
from math import radians


def d_theta_to_coords(d, theta):
    theta_rad = radians(theta)
    return np.array([d * np.cos(theta_rad), d * np.sin(theta_rad)])

def gen_s_state(
    res: List,
    cuda: bool = False,
    show: bool = False,
) -> np.ndarray:
    problem_ini = MicromagProblem(res, [500e-9, 125e-9, 3e-9], alpha=4.42e3, cuda=cuda)
    problem_ini.m0[:] = 1 / np.sqrt(3)

    Hyst_dir = 1 / (4 * np.pi * 1e-7) *  np.array([1, 1, 1])
    HextFct = lambda t: np.expand_dims(1e-9 - t, 0).T * Hyst_dir * np.expand_dims(t < 1e-9, 0).T

    t_out, M_out, _, _, _, _, _ = problem_ini.run_simulation(100e-9, 200, HextFct, 2000)
    M_sq_ini = np.squeeze(M_out, axis=2)

    if show:
        plt.clf()
        plt.plot(t_out, np.mean(M_sq_ini[:, :, 0], axis=1), 'rx')
        plt.plot(t_out, np.mean(M_sq_ini[:, :, 1], axis=1), 'gx')
        plt.plot(t_out, np.mean(M_sq_ini[:, :, 2], axis=1), 'bx')
        plt.show()

        plt.clf()
        plt.figure(figsize=(8, 2), dpi=80)
        s_state = M_sq_ini.reshape(200, res[1], res[0], 3)[-1]
        plt.quiver(s_state[:,:,0], s_state[:,:,1], pivot='mid')
        plt.show()

    return M_sq_ini[-1]


def gen_seq(
    m0_state,
    res,
    grid: List = [500e-9, 500e-9, 3e-9],
    H_ext: List = [0, 0, 0],
    t_steps: int = 500,
    t_per_step: float = 4e-12,
    cuda: bool = False,
    show: bool = False
) -> np.ndarray:
    problem = MicromagProblem(res, grid, alpha=4.42e3, gamma=2.211e5, cuda=cuda)
    problem.m0 = m0_state
    
    t_end = t_per_step * t_steps
    Hyst_dir = 1 / (4 * np.pi * 1e-7) *  np.array(H_ext) / 1000
    HextFct = lambda t: np.expand_dims(t > -1, 0).T * Hyst_dir

    t_out, M_out, _, _, _, _, _ = problem.run_simulation(t_end, t_steps, HextFct, 2000)
    M_sq = np.squeeze(M_out, axis=2)

    if show:
        plt.plot(t_out, np.mean(M_sq[:, :, 0], axis=1), 'rx')
        plt.plot(t_out, np.mean(M_sq[:, :, 1], axis=1), 'gx')
        plt.plot(t_out, np.mean(M_sq[:, :, 2], axis=1), 'bx')
        plt.show()

    return M_sq


def gen_std_prob_4(
    fname,
    res,
    grid: List = [500e-9, 500e-9, 3e-9],
    n_seq: int = 4,
    t_steps: int = 500,
    t_per_step: float = 4e-12,
    H_ext_angle: List = [0, 360],
    H_ext_norm: List = [0, 50],
    cuda: bool = False,
    seed: Optional[int] = None,
) -> None:
    rng = np.random.default_rng(seed)
    datapath = Path(__file__).parent.absolute() / '..' / 'data'
    if not datapath.exists(): datapath.mkdir()
    db = h5py.File(f'{datapath}/{fname}.h5', 'w')

    for i in range(n_seq):
        print(f'Generating sequence: {i+1}/{n_seq}')
        # TODO Load from database
        s_state = gen_s_state(res)

        # Random 2-D external field
        field = np.zeros(3)
        d = (H_ext_norm[1] - H_ext_norm[0]) * rng.random() + H_ext_norm[0]
        theta = (H_ext_angle[1] - H_ext_angle[0]) * rng.random() + H_ext_angle[0]
        field[:2] = d_theta_to_coords(d, theta)

        seq = gen_seq(
            m0_state=s_state,
            res=res,
            grid=grid,
            H_ext=field,
            t_steps=t_steps,
            t_per_step=t_per_step,
            cuda=cuda,
            show=False,
        )

        # Bring to shape T x X x Y x D
        seq = seq.reshape(t_steps, res[1], res[0], 3).swapaxes(1,2).moveaxis()
        seq = np.moveaxis(seq, -1, 1)

        g = db.create_group(str(i))
        g.create_dataset('sequence', seq)
        g.create_dataset('field', field)

    db.close()

# %%
