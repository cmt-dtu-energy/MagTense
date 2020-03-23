%% Benchmark MagTense, both Matlab and Fortran versions, in terms of computation time using Standard problem #4

%--- Use field X from the NIST example
NIST_field = 2;
res_arr = 1:9;

%--- Loop over using the code in both CPU and CUDA mode
for j = 1:2
    if j == 1
        j_l = false;
        j_s = 'no_CUDA';
        dir_f = 'CPU';
        dir_m = 'CPU';
        disp('Running no CUDA')
    else
        j_l = true;
        j_s = 'with_CUDA';
        dir_f = 'CUDA';
        dir_m = 'GPU';
         disp('Running with CUDA')
    end
    
    dem_thres = 0;
    dem_string = 'none';

    %--- Loop over the resolution of the problem
    for i = res_arr

        clearvars -except i j_l j_s dem_thres dem_string dir_f dir_m NIST_field res_arr
        
        %--- Vary the resolution
        resolution = [i*20,i*5,1];

        mu0 = 4*pi*1e-7;

        addpath('../../MEX_files');
        addpath('../../util');

        %% Setup the problem for the initial configuration
        tic
        problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

        problem = problem.setUseCuda( j_l );

        problem.alpha = 4.42e3;
        problem.gamma = 0;
        problem.dem_thres = dem_thres;
        problem.dem_appr = getMicroMagDemagApproximation(dem_string);

        %initial magnetization
        problem.m0(:) = 1/sqrt(3);

        %time grid on which to solve the problem
        problem = problem.setTime( linspace(0,100e-9,200) );
        problem.setTimeDis = int32(10);
        HystDir = 1/mu0*[-24.6,4.3,0]/1000 ;

        %time-dependent applied field
        HextFct = @(t) (1e-9-t)' .* HystDir .* (t<1e-9)';

        problem = problem.setHext( HextFct );

        solution = struct();
        %convert the class obj to a struct so it can be loaded into fortran
        prob_struct = struct(problem);

        solution = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
%         figure; M_end = squeeze(solution.M(end,:,:)); quiver(solution.pts(:,1),solution.pts(:,2),M_end(:,1),M_end(:,2)); axis equal; title('Starting state - Fortran')
        elapsedTime_part1 = toc;

        problem_initial = problem;
        solution_initial = solution;

        %% Setup problem for the time-dependent solver
        tic
        problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
        problem.alpha = 4.42e3 ;
        problem.gamma = 2.21e5 ;
        problem.dem_thres = dem_thres;
        problem.dem_appr = getMicroMagDemagApproximation(dem_string);
        problem = problem.setTime( linspace(0,1e-9,200) ); %
        problem.setTimeDis = int32(10);

        problem = problem.setUseCuda( j_l );

        if (NIST_field == 1)
            %field 1
            HystDir = 1/mu0*[-24.6,4.3,0]/1000 ;
        end
        if (NIST_field == 2)
            %field 2
            HystDir = 1/mu0*[-35.5,-6.3,0]/1000 ;
        end

        HextFct = @(t) (t>-1)' .*HystDir;
        problem = problem.setHext( HextFct );

        problem.m0(:) = solution.M(end,:,:);

        %convert the class obj to a struct so it can be loaded into fortran
        solution_t = struct();
        prob_struct = struct(problem);

        solution_t = MagTenseLandauLifshitzSolver_mex( prob_struct, solution_t );
        elapsedTime_part2 = toc

        memory_t = memory;
        gpuDevice_t = gpuDevice;

        problem_t = problem;

        mkdir(['Field_' num2str(NIST_field) '\' dir_f ])
        save(['Field_' num2str(NIST_field) '\' dir_f '\Solution_' num2str(i) '_x_20_5_1_' j_s '.mat'],'problem_initial','solution_initial','problem_t','solution_t','elapsedTime_part1','elapsedTime_part2','memory_t','gpuDevice_t')

             
        %% Run the Matlab version of MagTense
        
        %--- Fix MySim.SimulationName 
%         figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
        addpath('..\..\micromagnetism')
        Matlab_model_params.nGrid = resolution';
        Matlab_model_params.Field_dir = HystDir;
        Matlab_model_params.SaveTheResult = 1;
        Matlab_model_params.ShowTheResult = 0;
        Matlab_model_params.use_gpuArray = j_l;
        Matlab_model_params.use_single = 1;
        Matlab_model_params.DirectoryFilename = ['Field_' num2str(NIST_field)' '\Matlab_simulations_' dir_m '\Matlab_resolution_' num2str(i) '_x_20_5_1'];
        Script_3D_Std_Problem_4([],Matlab_model_params);
    end
end


