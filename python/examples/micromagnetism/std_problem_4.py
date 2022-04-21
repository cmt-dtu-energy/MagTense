import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

from magtense import magtense, micromag_problem


def std_prob_4(res=[64,16,1], NIST_field=1, use_CUDA=False, show=True):
    mu0 = 4*np.pi*1e-7

    ### Magnetization to s-state
    # Setup problem
    problem_ini = micromag_problem.DefaultMicroMagProblem(res)
    problem_ini.dem_appr = micromag_problem.get_micromag_demag_approx(None)
    problem_ini.set_use_CUDA(use_CUDA)

    problem_ini.grid_L = [500e-9, 125e-9, 3e-9]

    problem_ini.alpha = 4.42e3
    problem_ini.gamma = 0

    # Initial magnetization
    problem_ini.m0[:] = (1 / np.sqrt(3))

    # Time-dependent applied field
    t_end = 100e-9
    Hyst_dir = 1 / mu0 *  np.array([1, 1, 1])
    HextFct = lambda t: np.expand_dims(1e-9 - t, 0).T * Hyst_dir * np.expand_dims(t < 1e-9, 0).T
    problem_ini.set_Hext(HextFct, np.linspace(0, t_end, 2000))

    timesteps = 200
    M_ini_out = np.zeros(shape=(timesteps, problem_ini.m0.shape[0], 1, 3))
    t_ini_out = np.zeros(shape=(timesteps))

    for n_t in range(timesteps//50):
        dt = t_end / ( timesteps // 50 )
        problem_ini.set_time(np.linspace(n_t * dt, (n_t + 1) * dt , 50))    

        t, M, pts, _, _, _, _ = magtense.run_micromag_simulation(problem_ini)
        problem_ini.m0 = M[-1,:,0,:]
        M_ini_out[n_t*50:(n_t+1)*50] = M.copy()
        t_ini_out[n_t*50:(n_t+1)*50] = t.copy()

    M_sq_ini = np.squeeze(M_ini_out, axis=2)

    ### Time-dependent solver
    # Setup problem
    problem_dym = micromag_problem.DefaultMicroMagProblem(res)
    problem_dym.dem_appr = micromag_problem.get_micromag_demag_approx(None)
    problem_dym.set_use_CUDA(use_CUDA)

    problem_dym.grid_L = [500e-9, 125e-9, 3e-9]

    problem_dym.alpha = 4.42e3
    problem_dym.gamma = 2.21e5

    t_end = 1e-9

    # S-state as initial magnetization
    problem_dym.m0 = problem_ini.m0

    # Two applied fields of standard problem 4
    if NIST_field == 1:
        Hyst_dir = 1 / mu0 * np.array([-24.6, 4.3, 0]) / 1000
    if NIST_field == 2:
        Hyst_dir = 1 / mu0 * np.array([-35.5, -6.3, 0]) / 1000

    HextFct = lambda t: np.expand_dims(t > -1, 0).T * Hyst_dir
    problem_dym.set_Hext(HextFct, np.linspace(0, t_end, 2000))

    timesteps = 200
    M_dym_out = np.zeros(shape=(timesteps, problem_dym.m0.shape[0], 1, 3))
    t_dym_out = np.zeros(shape=(timesteps))

    for n_t in range(timesteps//50):
        dt = t_end / ( timesteps // 50 )
        problem_dym.set_time(np.linspace(n_t * dt, (n_t + 1) * dt , 50))    

        t, M, pts, _, _, _, _ = magtense.run_micromag_simulation(problem_dym)
        problem_dym.m0 = M[-1,:,0,:]
        M_dym_out[n_t*50:(n_t+1)*50] = M.copy()
        t_dym_out[n_t*50:(n_t+1)*50] = t.copy()
    
    M_sq_dym = np.squeeze(M_dym_out, axis=2)


    if show:
        plt.plot(t_dym_out, np.mean(M_dym_out[:,:,0,0], axis=1), 'rx')
        plt.plot(t_dym_out, np.mean(M_dym_out[:,:,0,1], axis=1), 'gx')
        plt.plot(t_dym_out, np.mean(M_dym_out[:,:,0,2], axis=1), 'bx')
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
    std_prob_4(NIST_field=1, show=True, use_CUDA=True)
