#%%
import numpy as np
import matplotlib.pyplot as plt
import pathlib

from magtense import magtense, micromag_problem


def std_prob_3(res=[10,10,10], L_loop=np.linspace(8,9,10), use_CUDA=False, 
               show=True, show_details=False, save=False):
    problem = micromag_problem.DefaultMicroMagProblem(res)
    problem.dem_appr = micromag_problem.get_micromag_demag_approx(None)
    problem.set_use_CUDA(use_CUDA)

    mu0 = 4*np.pi*1e-7
    problem.gamma = 0
    problem.Ms = 1000e3
    problem.K0 = 0.1 * 0.5 * mu0 * problem.Ms**2
    problem.A0 = 1.74532925199e-10
    problem.u_ea = np.zeros(shape=(np.prod(res), 3), dtype=np.float64, order='F')
    problem.u_ea[:,2] = 1
    lex = np.sqrt(problem.A0 / (0.5 * mu0 * problem.Ms**2))
    problem.setTimeDis = 10
    Hext_fct = lambda t: np.atleast_2d(t).T * [0, 0, 0]
    L_loop = [8.5]

    # Time-dependent alpha parameter, to ensure faster convergence
    problem.alpha = 1e3
    E_arr = np.zeros(shape=(4, len(L_loop), 2))

    ax = plt.figure().add_subplot(projection='3d')

    for i in range(len(L_loop)):
        print(f'ITERATION: {i} / {len(L_loop)}')
        for j in range(2):
            # Initial magnetization
            if j == 0:
                print('Flower state')
                problem.m0[:,0:2] = 0
                problem.m0[:,2] = 1
                t_end = 10e-9
                
            elif j == 1:
                print('Vortex state')
                xv = np.linspace(-1, 1, res[0])
                yv = np.linspace(-1, 1, res[1])
                zv = np.linspace(-1, 1, res[2])
                [x,y,z] = np.meshgrid(xv, yv, zv, indexing='ij')
                xvec =  np.sin(np.arctan2(z, x))
                yvec = -np.cos(np.arctan2(z, x))
                problem.m0[:,0] = xvec.swapaxes(0,2).reshape((-1))
                problem.m0[:,2] = yvec.swapaxes(0,2).reshape((-1))
                problem.m0 = problem.m0 / np.tile(np.expand_dims(np.sqrt(np.sum(problem.m0**2, axis=1)), axis=1), (1, 3))
                t_end = 200e-9

            # Time grid on which to solve the problem
            problem.set_time(np.linspace(0, t_end, 50))

            # Time-dependent applied field
            problem.set_Hext(Hext_fct, np.linspace(0, t_end, 2))
            
            problem.grid_L = np.array([lex, lex, lex]) * L_loop[i]

            t, M, pts, H_exc, H_ext, H_dem, H_ani = magtense.run_micromag_simulation(problem)

            M_sq = np.squeeze(M, axis=2)

            if show:
                for k, m in enumerate([M_sq[0,:,:], M_sq[-1,:,:]]):
                    plt.clf()
                    ax.quiver(pts[:,0], pts[:,1], pts[:,2], m[:,0], m[:,1], m[:,2])
                    plt.title('Fortran starting magnetization')
                    plt.savefig(f'{pathlib.Path(__file__).parent.resolve()}/state{j}_{k}.png')

                plt.clf()
                plt.plot(t, np.mean(M[:,:,0,0], axis=1), 'rx')
                plt.plot(t, np.mean(M[:,:,0,1], axis=1), 'gx')
                plt.plot(t, np.mean(M[:,:,0,2], axis=1), 'bx')
                plt.savefig(f'{pathlib.Path(__file__).parent.resolve()}/mag_components.png')
            
            # Calculate the energy terms
            E_exc = np.sum((1/2) * (M[:,:,0,0] * H_exc[:,:,0,0] \
                + M[:,:,0,1] * H_exc[:,:,0,1] \
                + M[:,:,0,2] * H_exc[:,:,0,2]), axis=1)

            E_ext = np.sum(M[:,:,0,0] * H_ext[:,:,0,0] \
                + M[:,:,0,1] * H_ext[:,:,0,1] \
                + M[:,:,0,2] * H_ext[:,:,0,2], axis=1)

            E_dem = np.sum((1/2) * (M[:,:,0,0] * H_dem[:,:,0,0] \
                + M[:,:,0,1] * H_dem[:,:,0,1] \
                + M[:,:,0,2] * H_dem[:,:,0,2]), axis=1)
            
            E_ani = np.sum((1/2) * (M[:,:,0,0] * H_ani[:,:,0,0] \
                + M[:,:,0,1] * H_ani[:,:,0,1] \
                + M[:,:,0,2] * H_ani[:,:,0,2]), axis=1)
            
            E_arr[:,i,j] = mu0 * np.array([E_exc[-1], E_ext[-1], E_dem[-1], E_ani[-1]])

            if show_details:
                plt.clf()
                plt.plot(t, mu0 * E_exc - mu0 * E_exc[0], '.')
                plt.plot(t, mu0 * E_ext - mu0 * E_ext[0], '.')
                plt.plot(t, mu0 * E_dem - mu0 * E_dem[0], '.')
                plt.plot(t, mu0 * E_ani - mu0 * E_ani[0], '.')
                plt.xlabel('Time [s]')
                plt.ylabel('Energy [-]')
                plt.legend({'E_{exc}', 'E_{ext}', 'E_{dem}', 'E_{ani}'}, 'Location', 'East')

    if show:
        plt.clf()
        plt.plot(L_loop, np.sum(E_arr[:,:,0], axis=0), '.')
        plt.plot(L_loop, np.sum(E_arr[:,:,1], axis=0), 'o')
        plt.xlabel('L [l_ex]')
        plt.ylabel('E [-]')
        plt.savefig(f'{pathlib.Path(__file__).parent.resolve()}/prob3.png')

    # print(f'Energy intersection: {np.interp(np.sum(E_arr[:,:,0], axis=0) - np.sum(E_arr[:,:,1], axis=0), L_loop, 0)})')

    return problem, (t, M, pts, H_exc, H_ext, H_dem, H_ani), E_arr, L_loop

#%%

if __name__ == '__main__':
    std_prob_3(show=True)