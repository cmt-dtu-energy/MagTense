function [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym,rel_int_error,mesh,GridInfo] = Standard_problem_4_unstructured_cart( mumag_field, options )
%STANDARD_PROBLEM_4_UNSTRUCTURED_CART 
%A function script to setup and simulate mumag standard problem 4 with an unstructured cartesian mesh loaded from a file
%
%Syntax:
%------
%   Standard_problem_4_unstructured_cart()
%   [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym,mesh,GridInfo] = Standard_problem_4( NIST_field, options)
%
%Description of syntax:
%------
%   Standard_problem_4_unstructured_cart() 
%       Uses the default parameters to solve mumag problem 4 and displays the results on screen
%
%   Standard_problem_4_unstructured_cart( NIST_field, options)
%       Takes 1 input argument which specifies the applied field. Additional options can also be specified
%
%   [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym,mesh,GridInfo] = Standard_problem_4_unstructured_cart( NIST_field, options)
%       As above but also returns the computation times, the problem setup file and the solution for both the initial and dynamical part of the problem as well as the mesh information
%
%Input arguments:
%------
%   NIST_field : Either 1 or 2
%       Determines if the first or second applied field specified in the mumag problem description is used (Default value is 1)
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
%   mesh : Struct
%      Contains information on the mesh used
%
%   GridInfo : Struct
%      Contains information on the mesh used
%
%Detailed description:
%-------
%   The script setups up and runs the mumag standard problem 4 for a prismal mesh.
%
%Version: 1.0.1
%Author:  Rasmus Bjørk
%Date:    2025.03.28
%
%See also: Standard_problem_4

arguments
    mumag_field (1,1) {mustBeInteger}                = 1             %--- Use either field 1 or field 2 from the mumag example
    options.use_CUDA {mustBeNumericOrLogical}       = true          %--- Use CUDA for the calculations
    options.ShowTheResult {mustBeNumericOrLogical}  = true          %--- Show the result
    options.use_CVODE {mustBeNumericOrLogical}      = false;        %--- Use CVODE for the numerical time evolution
end

mu0 = 4*pi*1e-7;

addpath('../../../MEX_files');
addpath('../../../util');

%% --------------------------------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------- MAGTENSE ---------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Setup the problem for the initial configuration
% load('Std_prob_4_unstructured_mesh_TEST.mat');
load('Std_prob_4_unstructured_mesh_grains_6_res_80_20_ref_2.mat');
% load('Std_prob_4_unstructured_mesh_grains_6_res_100_25_ref_3.mat');
cartesianUnstructuredMeshPlot(mesh.pos_out,mesh.dims_out,GridInfo,mesh.iIn);

%--- Setup the problem
resolution = [length(mesh.pos_out) 1 1];
disp(['Prisms N_grid = ' num2str(prod(resolution))])
problem_ini = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem_ini = problem_ini.setMicroMagGridType('unstructuredPrisms');

%--- Information on the grid
problem_ini.grid_pts    = mesh.pos_out;
problem_ini.grid_abc    = mesh.dims_out;

% Set specific flags on options
problem_ini = problem_ini.setMicroMagDemagApproximation('none');
problem_ini = problem_ini.setUseCuda( options.use_CUDA );
problem_ini = problem_ini.setUseCVODE( options.use_CVODE );

% Material properties
problem_ini.alpha = 4.42e3;
problem_ini.gamma = 0;
problem_ini.Ms = 8e5*ones(prod(resolution),1);
problem_ini.K0 = 0*zeros(prod(resolution),1);
problem_ini.A0 = 1.3e-11*ones(prod(resolution),1);
problem_ini.exch_weigh = 8;
problem_ini = problem_ini.setMicroMagExchMethod( 'DirectLaplacianNeumann'  );
problem_ini = problem_ini.setMicroMagExchInterpn( 'Extended'  ); 

%--- Calculate the exchange matrix
InteractionMatrices.GridInfo = GridInfo;
InteractionMatrices.X = GridInfo.Xel ;
InteractionMatrices.Y = GridInfo.Yel ;
InteractionMatrices.Z = GridInfo.Zel ;

[D2X,D2Y] = computeDifferentialOperatorsFromMesh_DirectLap(GridInfo,'extended',8,"DirectLaplacianNeumann", problem_ini.A0./( mu0 * problem_ini.Ms));
% [D2X,D2Y] = computeDifferentialOperatorsFromMesh_DirectLap(GridInfo,'extended',8,"DirectLaplacianNeumann");
InteractionMatrices.A2 = D2X + D2Y ;

%--- Convert the exchange matrix to sparse
problem_ini = problem_ini.setExchangeMatrixSparse( InteractionMatrices.A2 );

% Initial magnetization
problem_ini.m0(:,1) = 1/sqrt(3);
problem_ini.m0(:,2) = 1/sqrt(3);
problem_ini.m0(:,3) = 1/sqrt(3);
    
% Time points at which to return the solution
problem_ini = problem_ini.setTime( linspace(0,100e-9,200) );
problem_ini.setTimeDis = int32(100);

% The applied field as function of time
HystDir = 1/mu0*[1,1,1] ;
HextFct = @(t) (1e-9-t)' .* HystDir .* (t<1e-9)';
problem_ini = problem_ini.setHext( HextFct, linspace(0,100e-9,2000) );

%% Solve the initial configuration
% Convert the class obj to a struct so it can be loaded into fortran
solution_ini = struct();
prob_struct = struct(problem_ini);

tic
% solution_ini = problem_ini.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_ini );
[solution_ini, GridInfoFortran] = problem_ini.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_ini );
% [GridInfoFortran] = problem_ini.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_ini );

elapsedTime_part1 = toc
if (options.ShowTheResult)
    figure; 
    M_end = squeeze(solution_ini.M(end,:,:)); 
    quiver(solution_ini.pts(:,1),solution_ini.pts(:,2),M_end(:,1),M_end(:,2)); 
    axis equal; 
    title('Starting state of dynamical simulation')
end

%% Setup problem for the time-dependent solver'
% Use the initial problem to setup the dynamical part of the simulations
problem_dym = problem_ini;

% Calculate to 1 ns and save the results in 200 steps
problem_dym = problem_dym.setTime( linspace(0,1e-9,200) );
problem_dym.setTimeDis = int32(10);

% Material properties
problem_dym.gamma = 2.21e5 ;

% The external field applied
if (mumag_field == 1)
    %field 1
    HystDir = 1/mu0*[-24.6,4.3,0]/1000 ;
end
if (mumag_field == 2)
    %field 2
    HystDir = 1/mu0*[-35.5,-6.3,0]/1000 ;
end
HextFct = @(t) (t>-1)' .*HystDir;
problem_dym = problem_dym.setHext( HextFct, linspace(0,1e-9,2000) );

% Set the starting state to be that found in the initial part of the problem
problem_dym.m0(:) = solution_ini.M(end,:,:);

%% Solve the dynamic configuration
% Convert the class obj to a struct so it can be loaded into fortran
solution_dym = struct();
prob_struct = struct(problem_dym);

tic
solution_dym = problem_dym.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_dym );
elapsedTime_part2 = toc

[Mx,My,Mz,mx,my,mz] = computeMagneticMomentGeneralMesh(solution_dym.M,GridInfo.Volumes) ;

%% --------------------------------------------------------------------------------------------------------------------------------------
%% --------------------------------------------------------------------  mumag ----------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Compare with published solutions available from mumag webpage
t=1e-9*linspace(0,1,1000);
M_mumag = load(['../../../../documentation/examples_mumag_validation/Validation_standard_problem_4/Field_' num2str(mumag_field) '_mumag_mean_solution.txt']);

% Interpolate the MagTense solution to the mumag solution and calculate the relative error in percent
rel_int_error(1) = calculate_relative_integral_error(t,M_mumag(:,1),solution_dym.t,Mx);
rel_int_error(2) = calculate_relative_integral_error(t,M_mumag(:,3),solution_dym.t,My);
rel_int_error(3) = calculate_relative_integral_error(t,M_mumag(:,5),solution_dym.t,Mz);

%% --------------------------------------------------------------------------------------------------------------------------------------
%% ---------------------------------------------------------------  Plot the results ----------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
if (options.ShowTheResult)
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on

    %--- Plot the MagTense magnetization
    plot(fig1,solution_dym.t,Mx,'rx');
    plot(fig1,solution_dym.t,My,'gx');
    plot(fig1,solution_dym.t,Mz,'bx');

    %--- Plot the mumag solutions
    colours = [[1 0 0];[0 1 0];[0 0 1]];
    weak_colours = colours + ~colours*0.75;
    fill_ts=[t,fliplr(t)];  
    for j=1:3
        std_errors(1:2,:)=[M_mumag(:,(j-1)*2+1)+M_mumag(:,j*2), M_mumag(:,(j-1)*2+1)-M_mumag(:,j*2)]';
        interval = [std_errors(1,:),fliplr(std_errors(2,:))];
        plot(fig1,t,M_mumag(:,(j-1)*2+1),'color',colours(j,:))
        fill(fig1,fill_ts,interval,weak_colours(j,:),'linestyle','none')
    end

    legend(fig1,'MagTense M_x','MagTense M_y','MagTense M_z','\mu{}mag <M_x>','\mu{}mag \sigma{}(M_x)','\mu{}mag <M_y>','\mu{}mag \sigma{}(M_y)','\mu{}mag <M_z>','\mu{}mag \sigma{}(M_z)','Location','eastoutside');
    ylabel(fig1,'<M_i>/M_s')
    xlabel(fig1,'Time [ns]')
    xlim(fig1,[0 1e-9])
    figure(figure1)
end

end