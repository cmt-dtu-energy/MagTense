function [elapsedTime,problem,solution,E_arr,L_loop] = Standard_problem_3( resolution, L_loop, options )

arguments
    resolution (1,3) {mustBeInteger}                        = [10,10,10];       %--- [nx,ny,nz] of the grid
    L_loop (1,:) {mustBeNumeric}                            = linspace(8,9,10); %--- The side length values of the simulation cube
    options.use_CUDA {mustBeNumericOrLogical}               = true              %--- Use CUDA for the calculations
    options.ShowTheResult {mustBeNumericOrLogical}          = true              %--- Show the result
    options.ShowTheResultDetails {mustBeNumericOrLogical}   = false             %--- Show the magnetization at each L
    options.use_CVODE {mustBeNumericOrLogical}              = false;            %--- Use CVODE for the numerical time evolution
    options.use_minimizer {mustBeNumericOrLogical}          = true;            %--- Relax with the energy minimizer instead of integrating the LL equation in time
end

if (options.ShowTheResult)
    figure10= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig10 = axes('Parent',figure10,'Layer','top','FontSize',16); hold on; grid on; box on
end

if (options.ShowTheResultDetails)
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
    figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold on; grid on; box on
end

mu0 = 4*pi*1e-7;

addpath('../../../MEX_files');
addpath('../../../util');

%% --------------------------------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------- MAGTENSE ---------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Setup the problem for the initial configuration

     %--- Define parameters
    alpha = 1e3;
    gamma = 0;
    Ms = 1000e3;
    K0 = 0.1*1/2*mu0*Ms^2;
    A0 = 1.74532925199e-10;

    %--- Setup the problem
    problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
    problem = problem.setMicroMagDemagApproximation('none');
    problem = problem.setUseCuda( options.use_CUDA );
    problem = problem.setUseCVODE( options.use_CVODE );
    %--- The equilibrium is found either by integrating the LL equation in time at zero field over
    %--- the window set with setTime (the default 'Dynamic' solver of this example) or by the
    %--- energy minimizer ('Minimizer'), which ignores the time window and stops when the largest
    %--- torque is below problem.min_tol. The energies come back in solution.E in both cases.
    if options.use_minimizer
        problem = problem.setMicroMagSolver( 'Minimizer' );
    end
    problem.ReturnHall = int32(1);

    %--- Save the parameters
    problem.alpha = alpha;
    problem.gamma = gamma;
    problem.Ms = Ms*ones(prod(resolution),1);
    problem.K0 = K0*ones(prod(resolution),1);
    problem.A0 = A0;
    problem.u_ea = zeros( prod(resolution), 3 );
    problem.u_ea(:,3) = 1;

    %--- Information on the grid
    lex = sqrt(problem.A0/(1/2*mu0*Ms^2));
    
    problem.setTimeDis = int32(10);
    HextFct = @(t) (t)' .* [0,0,0];

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
                problem.m0(:,2) = 0.*xvec(:) ;
                problem.m0(:,3) = yvec(:) ;
                problem.m0 = problem.m0./repmat(sqrt(sum(problem.m0.^2,2)),1,3);
                
                t_end = 200e-9;
        end

        %time grid on which to solve the problem
        problem = problem.setTime( linspace(0,t_end,50) );

        %applied field (zero). The minimizer treats every row of the field table as a separate
        %field to relax at, so it gets a single row.
        if options.use_minimizer
            problem = problem.setHext( HextFct, 0 );
        else
            problem = problem.setHext( HextFct, linspace(0,t_end,2) );
        end
        
        problem.grid_L = [lex,lex,lex]*L_loop(i);%m
        tile_volumes = repmat(prod(problem.grid_L)/(prod(resolution)),prod(resolution),1);

        solution = struct();
        prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran

        solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
        
        [Mx,My,Mz,mx,my,mz] = computeMagneticMomentGeneralMesh(solution.M,tile_volumes) ;

        if (options.ShowTheResultDetails)
            if (j == 1)
                s1 = 1;
                s2 = 3;
            else
                s1 = 2;
                s2 = 4;
            end
            M_1 = squeeze(solution.M(1,:,:)); figure(figure2); subplot(2,2,s1); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_1(:,1),M_1(:,2),M_1(:,3)); axis equal; title('Starting magnetization')
            M_end = squeeze(solution.M(end,:,:)); figure(figure2); subplot(2,2,s2); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_end(:,1),M_end(:,2),M_end(:,3)); axis equal; title('Ending magnetization')

            %--- The time evolution of the magnetization
            plot(fig1,solution.t,Mx,'rd'); 
            plot(fig1,solution.t,My,'gd'); 
            plot(fig1,solution.t,Mz,'bd'); 
            xlabel(fig1,'Time [ns]')
            ylabel(fig1,'Reduced magnetization, m_i [-]')
            legend(fig1,'<m_x>','<m_y>','<m_z>')         
        end
        
        %--- The energy terms. solution.E holds them in J as (time, field, term) with the terms in
        %--- the order exchange, external, demag, anisotropy, evaluated in Fortran from the fields
        %--- the solver ran on. Divided by Km*V they are the reduced energies of the mumag problem
        %--- (this is what computeMagneticEnergy computed from the returned fields before).
        Km = 1/2*mu0*Ms^2;
        E_red = squeeze(solution.E(end,1,:)) / (Km*prod(problem.grid_L));
        E_arr(:,i,j) = [E_red(3) E_red(1) E_red(4) E_red(2)];   % dem, exc, ani, ext - the order this example has always used

        n_feval(i,j) = sum(solution.n_feval);
        if options.use_minimizer
            disp(['   Minimizer: ' num2str(sum(solution.n_feval)) ' field evaluations, ' num2str(sum(solution.min_iter)) ...
                  ' iterations, status ' num2str(solution.min_status(end)) ', E/(Km V) = ' num2str(sum(E_red))])
        else
            disp(['   LL relaxation: ' num2str(sum(solution.n_feval)) ' field evaluations, E/(Km V) = ' num2str(sum(E_red))])
        end
    end
end
elapsedTime = toc
disp(['Total field evaluations: ' num2str(sum(n_feval(:)))])

if (options.ShowTheResult)
    plot(fig10,L_loop,sum(E_arr(:,:,1),1),'.')
    plot(fig10,L_loop,sum(E_arr(:,:,2),1),'.')
    xlabel(fig10,'L [l_{ex}]')
    ylabel(fig10,'E [-]')
    legend(fig10,'Flower state','Vortex state');
end

end