import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

from magtense.micromag import MicromagProblem


def std_prob_4(res=[64,16,1], NIST_field=1, cuda=False, show=True):
    mu0 = 4 * np.pi * 1e-7
    grid_L = [500e-9, 125e-9, 3e-9]

    ### Magnetization to s-state
    problem_ini = MicromagProblem(res, grid_L, alpha=4.42e3, cuda=cuda)    
    problem_ini.m0[:] = 1 / np.sqrt(3)

    Hyst_dir = 1 / mu0 *  np.array([1, 1, 1])
    HextFct = lambda t: np.expand_dims(1e-9 - t, 0).T * Hyst_dir * np.expand_dims(t < 1e-9, 0).T

    _, M_out, pts, _, _, _, _ = problem_ini.run_simulation(100e-9, 200, HextFct, 2000)
    M_sq_ini = np.squeeze(M_out, axis=2)

    ### Time-dependent solver
    problem_dym = MicromagProblem(res, grid_L, alpha=4.42e3, gamma=2.21e5, cuda=cuda)
    # s-state as initial magnetization
    problem_dym.m0 = M_out[-1, :, 0, :]

    # Two applied fields of std problem 4
    if NIST_field == 1:
        ext_field = [-24.6, 4.3, 0]
    elif NIST_field == 2:
        ext_field = [-35.5, -6.3, 0]
    else:
        raise NotImplementedError()

    Hyst_dir = (np.array(ext_field) / 1000) / mu0
    HextFct = lambda t: np.expand_dims(t > -1, 0).T * Hyst_dir
    t_dym, M_out, pts, _, _, _, _ = problem_dym.run_simulation(1e-9, 200, HextFct, 2000)
    M_sq_dym = np.squeeze(M_out, axis=2)

    if show:
        plt.plot(t_dym, np.mean(M_sq_dym[:,:,0], axis=1), 'rx')
        plt.plot(t_dym, np.mean(M_sq_dym[:,:,1], axis=1), 'gx')
        plt.plot(t_dym, np.mean(M_sq_dym[:,:,2], axis=1), 'bx')
        plt.show()

        for k, m in enumerate([M_sq_ini[0,:,:], M_sq_ini[-1,:,:], M_sq_dym[-1,:,:]]):
            fig = go.Figure(
                data=go.Cone(
                    x=pts[:,0],
                    y=pts[:,1],
                    z=pts[:,2],
                    u=m[:,0], v=m[:,1], w=m[:,2]
                ),
                layout_title_text=f'State{k}',
            )
            # fig.update_layout(scene=dict(aspectratio=dict(x=res[0]/10, y=res[1]/10, z=res[2]/10)))
            fig.show()


if __name__ == '__main__':
    std_prob_4(NIST_field=1, show=False, cuda=False)
