import numpy as np
import matplotlib.pyplot as plt

from magtense.magtense import run_micromag_simulation
from magtense.micro_mag_problem import DefaultMicroMagProblem, get_micromag_demag_approx

def std_prob_3(res=[10,10,10], L_loop=np.linspace(8,9,10), use_CUDA=True, 
               show=True, show_details=False, save=False):
    problem = DefaultMicroMagProblem(res)
    problem.dem_appr = get_micromag_demag_approx('none')
    problem.set_use_CUDA(use_CUDA)

    mu0 = 4*np.pi*1e-7
    problem.gamma = 0
    problem.Ms = 1000e3
    problem.K0 = 0.1 * 1/2 * mu0 * problem.Ms**2
    problem.A0 = 1.74532925199e-10
    problem.u_ea = np.zeros(np.prod(res), 3)
    problem.u_ea[:,2] = 1
    lex = np.sqrt(problem.A0 / (1/2 * mu0 * problem.Ms**2))
    problem.setTimeDis = 10
    Hext_fct = lambda t: np.atleast_2d(t).T * [0, 0, 0]

    # Time-dependent alpha parameter, to ensure faster convergence
    problem.alpha = 1e3

    for i in range(len(L_loop)):
        print(f'ITERATION :{i} / {len(L_loop)}')
        for j in range(2):
            # Initial magnetization
            if j == 0:
                print('Flower state')
                problem.m0[:,0:2] = 0
                problem.m0[:,2] = 1
                t_end = 10e-9
                
            elif j == 2:
                print('Vortex state')
                xv = np.linspace(-1, 1, res[0])
                yv = np.linspace(-1, 1, res[1])
                zv = np.linspace(-1, 1, res[2])
                [x,y,z] = np.meshgrid(xv, yv, zv)
                xvec =  np.sin(np.arctan2(z, x))
                yvec = -np.cos(np.arctan2(z, x))
                problem.m0[:,0] = xvec[:]
                problem.m0[:,1] = yvec[:]
                problem.m0 = problem.m0 / np.tile(np.sqrt(sum(problem.m0**2, 2)), (1, 3))
                t_end = 200e-9

            # Time grid on which to solve the problem
            problem.set_time(np.linspace(0, t_end, 50))

            # Time-dependent applied field
            problem.set_Hext(Hext_fct, np.linspace(0, t_end, 2))
            
            problem.grid_L = [lex, lex, lex] * L_loop[i]

            t, M, pts, H_exc, H_ext, H_dem, H_ani = run_micromag_simulation(problem)
            E_arr = np.zeros(shape=solution.M.shape)

            if show:
                pass
                # TODO Plot
                # M_1 = np.squeeze(solution.M(1,:,:)); figure(figure3); subplot(2,2,1); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_1(:,1),M_1(:,2),M_1(:,3)); axis equal; title('Fortran starting magnetization')
                # M_end = np.squeeze(solution.M(end,:,:)); figure(figure3); subplot(2,2,3); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_end(:,1),M_end(:,2),M_end(:,3)); axis equal; title('Fortran ending magnetization')

                # plt.clf()
                # plt.plot(solution.t, np.mean(solution.M[:,:,0], 2), 'rx')
                # plt.plot(solution.t, np.mean(solution.M[:,:,1], 2), 'gx')
                # plt.plot(solution.t, np.mean(solution.M[:,:,2], 2), 'bx')
            
            # Calculate the energy terms
            E_exc = sum((1/2) * (solution.M[:,:,0] * solution.H_exc[:,:,0] \
                + solution.M[:,:,1] * solution.H_exc[:,:,1] \
                + solution.M[:,:,2] * solution.H_exc[:,:,2]), 2)

            E_ext = sum(solution.M[:,:,0] * solution.H_ext[:,:,0] \
                + solution.M[:,:,1] * solution.H_ext[:,:,1] \
                + solution.M[:,:,2] * solution.H_ext[:,:,2], 2)

            E_dem = sum((1/2) * (solution.M[:,:,0] * solution.H_dem[:,:,0] \
                + solution.M[:,:,1] * solution.H_dem[:,:,1] \
                + solution.M[:,:,2] * solution.H_dem[:,:,2]), 2)
            
            E_ani = sum((1/2) * (solution.M[:,:,0] * solution.H_ani[:,:,0] \
                + solution.M[:,:,1] * solution.H_ani[:,:,1] \
                + solution.M[:,:,2] * solution.H_ani[:,:,2]), 2)
            
            E_arr[:,i,j] = mu0 * [E_exc[-1], E_ext[-1], E_dem[-1], E_ani[-1]]

            if show_details:
                #TODO Plot
                pass
                # plt.clf()
                # plt.plot(solution.t, mu0 * E_exc - mu0 * E_exc[0], '.')
                # plt.plot(solution.t, mu0 * E_ext - mu0 * E_ext[0], '.')
                # plt.plot(solution.t, mu0 * E_dem - mu0 * E_dem[0], '.')
                # plt.plot(solution.t, mu0 * E_ani - mu0 * E_ani[0], '.')
                # plt.xlabel('Time [s]')
                # plt.ylabel('Energy [-]')
                # plt.legend({'E_{exc}','E_{ext}','E_{dem}','E_{ani}'},'Location','East')

    if show:
        plt.clf()
        plt.plot(L_loop, sum(E_arr[:,:,0], 1), '.')
        plt.plot(L_loop, sum(E_arr[:,:,0], 1), '.')
        plt.xlabel('L [l_ex]')
        plt.ylabel('E [-]')

    # print(f'Energy intersection: {interp1(sum(E_arr[:,:,0], 1) - sum(E_arr[:,:,1], 1), L_loop, 0)})')

    return problem, solution, E_arr, L_loop
