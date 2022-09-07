function [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym] = Standard_problem_4( NIST_field, resolution, options )
%STANDARD_PROBLEM_4 
%A function script to setup and simulate mumag standard problem 4
%
%Syntax:
%------
%   Standard_problem_4()
%   [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym] = Standard_problem_4( NIST_field, resolution, options)
%   [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym] = Standard_problem_4( NIST_field, resolution), use_CUDA, ShowTheResult )
%   [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym] = Standard_problem_4( NIST_field, resolution), use_CUDA, ShowTheResult )
%   function_with_docstring_following_the_template()
%   output_1 = function_with_docstring_following_the_template(var_1)
%   [output_1, output_2] = function_with_docstring_following_the_template(var_1, var_2)
%   [__] = function_with_docstring_following_the_template(__, options=val)
%   [__, complex_matrix] = function_with_docstring_following_the_template(__, complex_matrix=val)
%
%Description of syntax:
%------
%   Standard_problem_4() 
%       Uses the default parameters to solve mumag problem 4 and displays the results on screen
%
%   Standard_problem_4( NIST_field, resolution, options)
%       Takes 1 or 2 input argument which specifies the applied field and the resolution of the problem. Additional options can also be specified
%
%   [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym] = Standard_problem_4( NIST_field, resolution, options)
%       As above but also returns the computation times, the problem setup file and the solution for both the initial and dynamical part of the problem
%
%Input arguments:
%------
%   NIST_field : Either 1 or 2
%       Determines if the first or second applied field specified in the mumag problem description is used (Default value is 1)
%
%   resolution : Array of size 1x3 (Default value is [36,9,1])
%       The resolution of the prismal mesh used to solve the problem
%
%Options:
%-------
%   use_CUDA : Interpreted as a logical - Default is true
%       Determines if CUDA is used for the computation.
%
%   ShowTheResult : Interpreted as a logical - Default is true
%       Determines if the results are plotted or not.
%
%Output arguments:
%-------
%   elapsedTime_part1 : Double
%      The time takes to compute the initial part of the problem
%
%   elapsedTime_part2 : Double
%      The time takes to compute the dynamic part of the problem
%
%   problem_ini : Struct
%      A struct containing the MagTense problem setup for the initial part of the mumag standard problem 4
% 
%   solution_ini : Struct
%      A struct containing the MagTense solution for the initial part of the mumag standard problem 4
%
%   problem_dym : Struct
%      A struct containing the MagTense problem setup for the dynamic part of the mumag standard problem 4
% 
%   solution_dym : Struct
%      A struct containing the MagTense solution for the dynamic part of the mumag standard problem 4
%
%Detailed description:
%-------
%   The script setups up and runs the mumag standard problem 4 for a prismal mesh.
%
%Version: 1.0.0
%Author:  Rasmus Bjørk
%Date:    2022.06.21
%
%See also: Benchmark_using_Standard_problem_4

arguments
    NIST_field (1,1) {mustBeInteger} = 1                %--- Use either field 1 or field 2 from the mumag example
    resolution (1,3) {mustBeInteger} = [36,9,1];
    options.use_CUDA {mustBeNumericOrLogical} = true 
    options.ShowTheResult {mustBeNumericOrLogical} = true
    options.use_CVODE {mustBeNumericOrLogical} = false;
end

if (options.ShowTheResult)
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
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
problem_ini = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem_ini.grid_L = [500e-9,125e-9,3e-9];%m
problem_ini.nThreads = int32(8);

problem_ini = problem_ini.setMicroMagDemagApproximation('none');
problem_ini = problem_ini.setUseCuda( options.use_CUDA );
problem_ini = problem_ini.setUseCVODE( options.use_CVODE );

% loadFile = 'N_out_test.dat';
% if exist(loadFile,'file')
%     problem_ini = problem_ini.setLoadNFilename( loadFile );
% end
% problem_ini = problem_ini.setReturnNFilename( loadFile );

problem_ini.Ms = 8e5*ones(prod(resolution),1);
problem_ini.K0 = 0*zeros(prod(resolution),1);
problem_ini.A0 = 1.3e-11;

problem_ini.alpha = 4.42e3;
problem_ini.gamma = 0;

%initial magnetization
problem_ini.m0(:,1) = 1/sqrt(3);
problem_ini.m0(:,2) = 1/sqrt(3);
problem_ini.m0(:,3) = 1/sqrt(3);

% problem_ini.u_ea(:,:) = 0;
    
%time grid on which to solve the problem
problem_ini = problem_ini.setTime( linspace(0,100e-9,200) );
problem_ini.setTimeDis = int32(100);
HystDir = 1/mu0*[1,1,1] ;

%time-dependent applied field
HextFct = @(t) (1e-9-t)' .* HystDir .* (t<1e-9)';
problem_ini = problem_ini.setHext( HextFct, linspace(0,100e-9,2000) );


problem_ini = problem_ini.setConvergenceCheckTime( linspace(0,40e-9,2) );
% problem_ini.conv_tol = 1e-6;


%% Create Preset
solution_ini = struct();
%convert the class obj to a struct so it can be loaded into fortran
prob_struct = struct(problem_ini);

tic
solution_ini = problem_ini.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_ini );
elapsedTime_part1 = toc
if (options.ShowTheResult)
    figure; M_end = squeeze(solution_ini.M(end,:,:)); quiver(solution_ini.pts(:,1),solution_ini.pts(:,2),M_end(:,1),M_end(:,2)); axis equal; title('Starting state - Fortran')
end

%% Setup problem for the time-dependent solver
problem_dym = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem_dym.alpha = 4.42e3 ;
problem_dym.gamma = 2.21e5 ;
problem_dym = problem_dym.setUseCuda( options.use_CUDA );
problem_dym = problem_dym.setUseCVODE( options.use_CVODE );
problem_dym = problem_dym.setMicroMagDemagApproximation('none');
problem_dym = problem_dym.setTime( linspace(0,1e-9,200) );
problem_dym.setTimeDis = int32(10);
problem_dym.Ms = 8e5*ones(prod(resolution),1);
problem_dym.K0 = 0*zeros(prod(resolution),1);
problem_dym.A0 = 1.3e-11;
problem_dym.nThreads = problem_ini.nThreads;

if (NIST_field == 1)
    %field 1
    HystDir = 1/mu0*[-24.6,4.3,0]/1000 ;
end
if (NIST_field == 2)
    %field 2
    HystDir = 1/mu0*[-35.5,-6.3,0]/1000 ;
end

HextFct = @(t) (t>-1)' .*HystDir;
problem_dym = problem_dym.setHext( HextFct, linspace(0,1e-9,2000) );

problem_dym.m0(:) = solution_ini.M(end,:,:);

%convert the class obj to a struct so it can be loaded into fortran
solution_dym = struct();
prob_struct = struct(problem_dym);

tic
solution_dym = problem_dym.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_dym );
elapsedTime_part2 = toc

if (options.ShowTheResult)
    plot(fig1,solution_dym.t,mean(solution_dym.M(:,:,1),2),'rx'); 
    plot(fig1,solution_dym.t,mean(solution_dym.M(:,:,2),2),'gx'); 
    plot(fig1,solution_dym.t,mean(solution_dym.M(:,:,3),2),'bx'); 
    % figure; hold all; for i=2:4; plot(problem.Hext(:,1),problem.Hext(:,i),'.'); end;
end


%% --------------------------------------------------------------------------------------------------------------------------------------
%% --------------------------------------------------------------------  mumag -----------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Compare with published solutions available from mumag webpage
if (options.ShowTheResult)
    t=linspace(0,1,1000);
    data = load(['Published_solutions_field' num2str(NIST_field)]);
    mean_avg=data.mean_avg; err=data.err;
    colours = [[1 0 0];[0 1 0];[0 0 1]];
    weak_colours = colours + ~colours*0.75;
    fill_ts=[t,fliplr(t)];  
    for j=1:3
        std_errors{j}(1:2,:)=[mean_avg(1,:,j)+err(1,:,j);mean_avg(1,:,j)-err(1,:,j)];
        interval = [std_errors{j}(1,:),fliplr(std_errors{j}(2,:))];
        fill(fig1,1e-9*fill_ts,interval,weak_colours(j,:),'linestyle','none')
        plot(fig1,1e-9*t,mean_avg(1,:,j),'color',colours(j,:))
    end

    legend(fig1,'Fortran Mx','Fortran My','Fortran Mz','Matlab Mx','Matlab My','Matlab Mz','\mu{}mag \sigma{}(Mx)','\mu{}mag <Mx>','\mu{}mag \sigma{}(My)','\mu{}mag <My>','\mu{}mag \sigma{}(Mz)','\mu{}mag <Mz>');
    ylabel(fig1,'<M_i>/M_s')
    xlabel(fig1,'Time [ns]')

    xlim(fig1,[0 1e-9])
    figure(figure1)
end

end