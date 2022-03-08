import numpy as np
import matplotlib.pyplot as plt

from magtense import magtense, micromag_problem


def transformer_data(res=[36,36,1], use_CUDA=False, use_CVODE=False, show=True):
    rng = np.random.default_rng(0)
    mu0 = 4 * np.pi * 1e-7

    problem = micromag_problem.DefaultMicroMagProblem(res)
    problem.set_use_CUDA(use_CUDA)
    problem.set_use_CVODE(use_CVODE)

    problem.grid_L = [500e-9,500e-9,3e-9]
    mu0 = 4*np.pi*1e-7
    problem.alpha = 4.42e3
    problem.gamma = 2.21e5
    problem.dem_appr = micromag_problem.get_micromag_demag_approx(None)
    problem.setTimeDis = 10
    problem.set_time(np.linspace(0,2e-9,500))

    Hyst_dir = 1 / mu0 * np.array([-25, 5, 0]) / 1000
    HextFct = lambda t: np.expand_dims(t > -1, 0).T * Hyst_dir
    problem.set_Hext(HextFct, np.linspace(0,2e-9,2000))
    
    for i in range(problem.m0.shape[0]):
        v = 2 * rng.random((3)) - 1
        problem.m0[i] = 1 / np.linalg.norm(v) * v

    t, M, _, _, _, _, _ = magtense.run_micromag_simulation(problem)
    

    if show:
        plt.plot(t, np.mean(M[:,:,0,0], axis=1), 'rx')
        plt.plot(t, np.mean(M[:,:,0,0], axis=1), 'gx')
        plt.plot(t, np.mean(M[:,:,0,0], axis=1), 'bx')
        plt.show()

        # TODO Plot start and end state

if __name__ == '__main__':
    transformer_data(show=True)
