#%%
import numpy as np

from magtense.micromag import MicromagProblem
from magtense.utils import plot_M_avg_seq, plot_M_thin_film


def std_prob_4(res=[64,16,1], NIST_field=1, cuda=False, show=True):
    mu0 = 4 * np.pi * 1e-7
    grid_L = [500e-9, 125e-9, 3e-9]

    ### Magnetization to s-state
    problem_ini = MicromagProblem(res, grid_L, m0=1/np.sqrt(3), alpha=4.42e3, cuda=cuda)
    h_ext = np.array([1, 1, 1]) / mu0
    h_ext_fct = lambda t: np.expand_dims(np.where(t < 1e-9, 1e-9 - t, 0), axis=1) * h_ext

    _, M_out, _, _, _, _, _ = problem_ini.run_simulation(100e-9, 200, h_ext_fct, 2000)
    M_sq_ini = np.squeeze(M_out, axis=2)

    ### Time-dependent solver
    problem_dym = MicromagProblem(res, grid_L, m0=M_sq_ini[-1], alpha=4.42e3, gamma=2.21e5, cuda=cuda)

    # Two applied external fields of std problem 4
    if NIST_field == 1:
        h_ext_nist = np.array([-24.6, 4.3, 0])
    elif NIST_field == 2:
        h_ext_nist = np.array([-35.5, -6.3, 0])
    else:
        raise NotImplementedError()

    h_ext_fct = lambda t: np.expand_dims(t > -1, axis=1) * (h_ext_nist / 1000 / mu0)
    t_dym, M_out, _, _, _, _, _ = problem_dym.run_simulation(1e-9, 200, h_ext_fct, 2000)

    if show:
        M_sq_dym = np.squeeze(M_out, axis=2)
        plot_M_avg_seq(t_dym, M_sq_dym)
        plot_M_thin_film(M_sq_dym[0], res, 'Start state')
        plot_M_thin_film(M_sq_dym[-1], res, 'Final state')
    
# %%

if __name__ == '__main__':
    std_prob_4(NIST_field=1, show=True, cuda=True)
