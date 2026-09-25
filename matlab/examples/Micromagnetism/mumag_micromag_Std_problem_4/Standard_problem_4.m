function [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym,rel_int_error,GridInfo] = Standard_problem_4( mumag_field, resolution, options )
%STANDARD_PROBLEM_4
%A function script to setup and simulate mumag standard problem 4
%
%Syntax:
%------
%   Standard_problem_4()
%   [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym,rel_int_error,GridInfo] = Standard_problem_4( mumag_field, resolution, options)
%
%Description of syntax:
%------
%   Standard_problem_4()
%       Uses the default parameters to solve mumag problem 4 and displays the results on screen
%
%   Standard_problem_4( mumag_field, resolution, options)
%       Takes 1 or 2 input argument which specifies the applied field and the resolution of the problem. Additional options can also be specified
%
%   [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym,rel_int_error,GridInfo] = Standard_problem_4( mumag_field, resolution, options)
%       As above but also returns the computation times, the problem setup file and the solution for both the initial and dynamical part of the problem
%
%Input arguments:
%------
%   mumag_field : Either 1 or 2
%       Determines if the first or second applied field specified in the mumag problem description is used (Default value is 1)
%
%   resolution : Array of size 1x3 (Default value is [36,9,1])
%       The resolution of the uniform grid used to solve the problem. Ignored for the unstructured mesh
%
%Options:
%-------
%   mesh_type : Either 'uniform' or 'unstructuredPrisms' - Default is 'uniform'
%       The type of mesh to run on. 'uniform' is a regular grid of resolution(1) x resolution(2) x resolution(3) cells.
%       'unstructuredPrisms' is a grid of unstructured Cartesian prisms read from the text file mesh_file (one row per
%       prism: centre [x,y,z] and side lengths [a,b,c] in metres). The exchange operator that MagTense builds from the
%       mesh in the first stage is handed to the second stage with setExchangeMatrixCOO, so the mesh is analysed only once.
%
%   mesh_file : Path to a text file
%       The unstructured Cartesian mesh, used with mesh_type = 'unstructuredPrisms'.
%
%   use_CUDA : Interpreted as a logical - Default is true
%       Determines if CUDA is used for the computation.
%
%   use_CVODE : Interpreted as a logical - Default is false
%       Determines if CVODE is used for the numerical time evolution.
%
%   use_AvgN : Interpreted as a logical - Default is true
%       Use the averaged prism tensor for the demag field.
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
%   rel_int_error : Array
%      The integrated difference between the NIST published solutions and the MagTense computed solution, relative to the
%      integral of the published solution and in percent, for the three components of the average magnetization.
%
%   GridInfo : Struct
%      The information on the mesh that MagTense built, including the cell volumes and the exchange matrix. Only filled
%      for the unstructured mesh.
%
%Detailed description:
%-------
%   The script setups up and runs the mumag standard problem 4 on a uniform grid or on an unstructured prism mesh.
%
%Version: 1.1.0
%Author:  Rasmus Bjørk
%Date:    2026.09.22

arguments
    mumag_field (1,1) {mustBeInteger}                = 1            %--- Use either field 1 or field 2 from the mumag example
    resolution (1,3) {mustBeInteger}                = [36,9,1];     %--- [nx,ny,nz] of the uniform grid (ignored for the unstructured mesh)
    options.mesh_type {mustBeMember(options.mesh_type,{'uniform','unstructuredPrisms'})} = 'uniform' %--- The type of mesh to run on
    options.mesh_file (1,:) char                    = '../../../../documentation/examples_mumag_validation/Validation_standard_problem_4/Std_prob_4_unstructured_mesh_grains_6_res_80_20_ref_2.txt' %--- The unstructured Cartesian mesh
    % options.mesh_file (1,:) char                    = '../../../../documentation/examples_mumag_validation/Validation_standard_problem_4/Std_prob_4_unstructured_mesh_grains_6_res_100_25_ref_3.txt' %--- The finer unstructured Cartesian mesh
    options.use_CUDA {mustBeNumericOrLogical}       = true          %--- Use CUDA for the calculations
    options.ShowTheResult {mustBeNumericOrLogical}  = true          %--- Show the result
    options.use_CVODE {mustBeNumericOrLogical}      = false;        %--- Use CVODE for the numerical time evolution
    options.use_AvgN {mustBeNumericOrLogical}       = true;         %--- Use the averaged prism tensor for demag field calculations
end

mu0 = 4*pi*1e-7;

addpath('../../../MEX_files');
addpath('../../../util');

%% --------------------------------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------- MAGTENSE ---------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Setup the mesh
switch options.mesh_type
    case 'uniform'
        % Construct a default problem, with a grid with size (nx,ny,nz)
        problem_ini = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

    case 'unstructuredPrisms'
        [data] = load(options.mesh_file);
        pos_out  = data(:,1:3);
        dims_out = data(:,4:6);

        resolution = [length(pos_out) 1 1];
        disp(['Prisms N_grid = ' num2str(prod(resolution))])
        problem_ini = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
        problem_ini = problem_ini.setMicroMagGridType('unstructuredPrisms');

        %--- Information on the grid
        problem_ini.grid_pts    = pos_out;
        problem_ini.grid_abc    = dims_out;
end
problem_ini.grid_L = [500e-9,125e-9,3e-9]; %m

%% Setup the problem for the initial configuration
% Set specific flags on options
problem_ini = problem_ini.setMicroMagDemagApproximation('none');
problem_ini = problem_ini.setUseCuda( options.use_CUDA );
problem_ini = problem_ini.setUseCVODE( options.use_CVODE );
problem_ini.useAvgN = int32(options.use_AvgN);

% Material properties
problem_ini.alpha = 4.42e3;
problem_ini.gamma = 0;
problem_ini.Ms = 8e5*ones(prod(resolution),1);
problem_ini.K0 = 0*zeros(prod(resolution),1);
problem_ini.A0 = 1.3e-11*ones(prod(resolution),1);

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
if strcmp(options.mesh_type,'uniform')
    solution_ini = problem_ini.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_ini );
    GridInfo = struct();
    tile_volumes = [];   % all tiles have the same volume, so computeMagneticMomentGeneralMesh takes the plain mean
else
    %--- The second output is the GridInfo that MagTense built from the mesh. Its Volumes are what
    %--- the magnetic moment has to be weighted by, and its exchange matrix is reused below
    [solution_ini, GridInfo] = problem_ini.MagTenseLandauLifshitzSolver_mex( prob_struct, solution_ini );
    tile_volumes = GridInfo.Volumes;
end
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

% Pass the exchange matrix MagTense built from the unstructured mesh on to the dynamic problem,
% so the mesh is analysed only once
if strcmp(options.mesh_type,'unstructuredPrisms')
    problem_dym = problem_dym.setExchangeMatrixCOO( GridInfo.ExchMat_nr, GridInfo.ExchMat_nc ...
                                        , GridInfo.ExchMat_r, GridInfo.ExchMat_c, GridInfo.ExchMat_v );
end

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

[Mx,My,Mz,mx,my,mz] = computeMagneticMomentGeneralMesh(solution_dym.M,tile_volumes);

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
    if strcmp(options.mesh_type,'unstructuredPrisms')
        cartesianUnstructuredMeshPlot(pos_out,dims_out,GridInfo);
    end

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
    title(fig1,['Standard problem 4, Field ' num2str(mumag_field) ', mesh: ' options.mesh_type])
    figure(figure1)
end

end
