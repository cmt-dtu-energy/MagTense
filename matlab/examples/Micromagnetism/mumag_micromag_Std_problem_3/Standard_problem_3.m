function [elapsedTime,problem,solution,E_arr,L_loop] = Standard_problem_3( resolution, use_CUDA, ShowTheResult, ShowTheResultDetails, L_loop )

clearvars -except  resolution use_CUDA ShowTheResult SaveTheResult ShowTheResultDetails RunMatlab L_loop 
close all

if ~exist('resolution','var')
    resolution = [10,10,10];
end
if ~exist('use_CUDA','var')
    use_CUDA = true;
end
if ~exist('ShowTheResult','var')
    ShowTheResult = 1;
end
if ~exist('ShowTheResultDetails','var')
    ShowTheResultDetails = 0;
end
if ~exist('L_loop','var')
    L_loop = linspace(8,9,10);
end

if use_CUDA
    dir_m = 'GPU';
else
    dir_m = 'CPU';
end

if (ShowTheResult)
    figure10= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig10 = axes('Parent',figure10,'Layer','top','FontSize',16); hold on; grid on; box on
end

if (ShowTheResultDetails)
    RunMatlab = 0;
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
    figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold on; grid on; box on
    figure3= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig3 = axes('Parent',figure3,'Layer','top','FontSize',16); hold on; grid on; box on
end
mu0 = 4*pi*1e-7;

addpath('../../../MEX_files');
addpath('../../../util');

%% --------------------------------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------- MAGTENSE ---------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Setup the problem for the initial configuration
%takes the size of the grid as arguments (nx,ny,nz) and a function handle
%that produces the desired field (if not present zero applied field is inferred)
problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

problem.dem_appr = getMicroMagDemagApproximation('none');
problem = problem.setUseCuda( use_CUDA );

problem.gamma = 0;
Ms = 1000e3;
K0 = 0.1*1/2*mu0*Ms^2;
problem.Ms = Ms*ones(prod(resolution),1);
problem.K0 = K0*ones(prod(resolution),1);
problem.A0 = 1.74532925199e-10;

problem.u_ea = zeros( prod(resolution), 3 );
problem.u_ea(:,3) = 1;
lex = sqrt(problem.A0/(1/2*mu0*Ms^2));
problem.setTimeDis = int32(10);
HextFct = @(t) (t)' .* [0,0,0];

%time-dependent alpha parameter, to ensure faster convergence
problem.alpha = 1e3;

tic
for i = 1:length(L_loop)  
    disp(['ITERATION :',num2str(i),'/',num2str(length(L_loop))])
    for j = 1:2
        switch j
            %initial magnetization
            case 1
                disp('Flower state') ;
                problem.m0(:) = 0;
                problem.m0(:,3) = 1 ;
                
                t_end = 10e-9;  
                
            case 2
                disp('Vortex state') ;
                [x,y,z]=ndgrid(linspace(-1,1,resolution(1)),linspace(-1,1,resolution(2)),linspace(-1,1,resolution(3)));
                xvec =  sin(atan2(z,x));
                yvec = -cos(atan2(z,x));
                problem.m0(:,1) = xvec(:) ;
                problem.m0(:,3) = yvec(:) ;
                problem.m0 = problem.m0./repmat(sqrt(sum(problem.m0.^2,2)),1,3);
                
                t_end = 200e-9;
        end

        %time grid on which to solve the problem
        problem = problem.setTime( linspace(0,t_end,50) );

        %time-dependent applied field
        problem = problem.setHext( HextFct, linspace(0,t_end,2) );
        
        problem.grid_L = [lex,lex,lex]*L_loop(i);%m

        solution = struct();
        prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran

        solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

        if (ShowTheResultDetails)
            M_1 = squeeze(solution.M(1,:,:)); figure(figure3); subplot(2,1,1); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_1(:,1),M_1(:,2),M_1(:,3)); axis equal; title('Fortran starting magnetization')
            M_end = squeeze(solution.M(end,:,:)); figure(figure3); subplot(2,1,2); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_end(:,1),M_end(:,2),M_end(:,3)); axis equal; title('Fortran ending magnetization')

            plot(fig1,solution.t,mean(solution.M(:,:,1),2),'rx'); 
            plot(fig1,solution.t,mean(solution.M(:,:,2),2),'gx'); 
            plot(fig1,solution.t,mean(solution.M(:,:,3),2),'bx'); 
        end
        
        %--- Calculate the energy terms
        E_exc = sum((1/2)*(solution.M(:,:,1).*solution.H_exc(:,:,1) + solution.M(:,:,2).*solution.H_exc(:,:,2) + solution.M(:,:,3).*solution.H_exc(:,:,3)),2) ;        
        E_ext = sum(      (solution.M(:,:,1).*solution.H_ext(:,:,1) + solution.M(:,:,2).*solution.H_ext(:,:,2) + solution.M(:,:,3).*solution.H_ext(:,:,3)),2) ;
        E_dem = sum((1/2)*(solution.M(:,:,1).*solution.H_dem(:,:,1) + solution.M(:,:,2).*solution.H_dem(:,:,2) + solution.M(:,:,3).*solution.H_dem(:,:,3)),2) ;
        E_ani = sum((1/2)*(solution.M(:,:,1).*solution.H_ani(:,:,1) + solution.M(:,:,2).*solution.H_ani(:,:,2) + solution.M(:,:,3).*solution.H_ani(:,:,3)),2) ;
                   
        %--- Normalize the energy to the volume elements
        dV = prod(problem.grid_L./resolution);
        E_exc = mu0*E_exc*Ms*dV;  % [J]
        E_ext = mu0*E_ext*Ms*dV;  % [J]
        E_dem = mu0*E_dem*Ms*dV;  % [J]
        E_ani = mu0*E_ani*Ms*dV;  % [J]
        
        %--- The anistropy energy must be added the constant volume term of the integral of K0 over the volume 
        %--- See "Nonlinear Magnetization Dynamics in Thin-films and Nanoparticles" by Massimiliano d'Aquino, Eq. (1.43)
        E_ani = E_ani + K0*prod(problem.grid_L);
        
        %--- Calcute the energy densities
        E_exc_dens = E_exc/prod(problem.grid_L);  % [J/m^3]
        E_ext_dens = E_ext/prod(problem.grid_L);  % [J/m^3]
        E_dem_dens = E_dem/prod(problem.grid_L);  % [J/m^3]
        E_ani_dens = E_ani/prod(problem.grid_L);  % [J/m^3]
        
        Km = (1/2)*(Ms^2)*mu0;  % [J/m^3] 
        %--- See Eq. (29) from "Standard Problems in Micromagnetics", Donald Porter and Michael Donahue (2018)
        
        %--- Reduced energies, i.e. normalized by Km
        E_exc_red = E_exc_dens./Km;
        E_ext_red = E_ext_dens./Km;
        E_dem_red = E_dem_dens./Km;
        E_ani_red = E_ani_dens./Km;
        
        E_arr(:,i,j) = [E_dem_red(end) E_exc_red(end) E_ani_red(end) E_ext_red(end)];
        
        if (ShowTheResultDetails)
            plot(fig2,solution.t,mu0*E_exc-mu0*E_exc(1),'.')
            plot(fig2,solution.t,mu0*E_ext-mu0*E_ext(1),'.')
            plot(fig2,solution.t,mu0*E_dem-mu0*E_dem(1),'.')
            plot(fig2,solution.t,mu0*E_ani-mu0*E_ani(1),'.')
            xlabel(fig2,'Time [s]')
            ylabel(fig2,'Energy [-]')
            legend(fig2,{'E_{exc}','E_{ext}','E_{dem}','E_{ani}'},'Location','East');
        end
    end
end
elapsedTime = toc

if (ShowTheResult)
    plot(fig10,L_loop,sum(E_arr(:,:,1),1),'.')
    plot(fig10,L_loop,sum(E_arr(:,:,2),1),'.')
    xlabel(fig10,'L [l_{ex}]')
    ylabel(fig10,'E [-]')
    legend(fig10,'Flower state','Vortex state');
end

end