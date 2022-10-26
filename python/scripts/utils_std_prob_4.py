#%%
import numpy as np
import matplotlib.pyplot as plt

from magtense.micromag import MicromagProblem
from typing import List


def gen_s_state(
    res: List,
    grid_size: List,
    cuda: bool = False,
    show: bool = False,
) -> np.ndarray:
    problem_ini = MicromagProblem(
        res=res,
        grid_L=grid_size,
        m0=1/np.sqrt(3),
        alpha=4.42e3,
        cuda=cuda
    )

    h_ext = np.array([1, 1, 1]) / (4 * np.pi * 1e-7)
    h_ext_fct = lambda t: np.expand_dims(np.where(t < 1e-9, 1e-9 - t, 0), axis=1) * h_ext
    
    t_out, M_out, _, _, _, _, _ = problem_ini.run_simulation(100e-9, 200, h_ext_fct, 2000)
    M_sq_ini = np.squeeze(M_out, axis=2)

    if show:
        plt.clf()
        plt.plot(t_out, np.mean(M_sq_ini[:, :, 0], axis=1), 'rx')
        plt.plot(t_out, np.mean(M_sq_ini[:, :, 1], axis=1), 'gx')
        plt.plot(t_out, np.mean(M_sq_ini[:, :, 2], axis=1), 'bx')
        plt.show()

        plt.clf()
        plt.figure(figsize=(8, 2), dpi=80)
        s_state = np.reshape(M_sq_ini[-1], (res[1], res[0], 3))
        plt.quiver(s_state[:,:,0], s_state[:,:,1], pivot='mid')
        plt.show()

    return M_sq_ini[-1]


def gen_seq(
    m0_state: np.ndarray,
    res: List,
    grid_size: List,
    h_ext: List = [0, 0, 0],
    t_steps: int = 500,
    t_per_step: float = 4e-12,
    cuda: bool = False,
    show: bool = False
) -> np.ndarray:
    problem = MicromagProblem(
        res=res,
        grid_L=grid_size,
        m0=m0_state,
        alpha=4.42e3,
        gamma=2.211e5,
        cuda=cuda
    )
    
    t_end = t_per_step * t_steps
    h_ext = np.array(h_ext) / 1000 / (4 * np.pi * 1e-7)
    h_ext_fct = lambda t: np.expand_dims(t > -1, axis=1) * h_ext

    t_out, M_out, _, _, _, _, _ = problem.run_simulation(t_end, t_steps, h_ext_fct, 2000)
    M_sq = np.squeeze(M_out, axis=2)

    if show:
        plt.plot(t_out, np.mean(M_sq[:, :, 0], axis=1), 'rx')
        plt.plot(t_out, np.mean(M_sq[:, :, 1], axis=1), 'gx')
        plt.plot(t_out, np.mean(M_sq[:, :, 2], axis=1), 'bx')
        plt.show()

    return M_sq

# %%
