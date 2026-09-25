function [elapsedTime,problem,solution,E_arr,L_loop] = Standard_problem_3( resolution, L_loop, options )
%STANDARD_PROBLEM_3
% Mumag standard problem 3: the flower/vortex energy crossover of a cube as a function of its
% edge length in units of the exchange length.
%
% The problem can be run on three kinds of mesh, selected with options.mesh_type:
%
%   'uniform'            a regular grid of resolution(1) x resolution(2) x resolution(3) cells.
%                        This is the default and the only mesh type that uses the resolution
%                        argument.
%   'unstructuredPrisms' a grid of unstructured Cartesian prisms read from the text file
%                        options.mesh_file (one row per prism: centre [x,y,z] and side lengths
%                        [a,b,c] in metres). The mesh is scaled to each edge length in L_loop.
%   'tetrahedron'        a tetrahedral mesh generated at each edge length with CreateTetraMesh,
%                        with options.mesh_res_param cells per edge length. The nodes and the
%                        connectivity are handed to MagTense with setMicroMagGridTetrahedron,
%                        and MagTense analyses the mesh and builds the exchange operator itself.
%                        Requires the PDE Toolbox, for generateMesh in CreateTetraMesh.
%
% For the two unstructured meshes the mesh is the only thing given to MagTense: it computes the
% cell volumes, the exchange operator and the energies itself. The energies of all mesh types
% come back in solution.E, evaluated in Fortran with the cell volumes of the mesh.
%
% E_arr(:,i,j) holds the reduced energies [demag, exchange, anisotropy, external] at the edge
% length L_loop(i) for the flower (j=1) and the vortex (j=2) state.

arguments
    resolution (1,3) {mustBeInteger}                        = [10,10,10];       %--- [nx,ny,nz] of the uniform grid (ignored for the other mesh types)
    L_loop (1,:) {mustBeNumeric}                            = linspace(8,9,10); %--- The side length values of the simulation cube [l_ex]
    options.mesh_type {mustBeMember(options.mesh_type,{'uniform','unstructuredPrisms','tetrahedron'})} = 'uniform' %--- The type of mesh to run on
    options.mesh_file (1,:) char                            = '../../../../documentation/examples_mumag_validation/Validation_standard_problem_3/Std_prob_3_unstructured_cartesian_grains_9_mesh_4_ref_2.txt' %--- The unstructured Cartesian mesh
    % options.mesh_file (1,:) char                            = '../../../../documentation/examples_mumag_validation/Validation_standard_problem_3/Std_prob_3_unstructured_cartesian_grains_9_mesh_7_ref_2.txt' %--- The unstructured Cartesian mesh
    options.mesh_res_param (1,1) {mustBePositive}           = 5                 %--- Tetrahedra per edge length
    options.use_CUDA {mustBeNumericOrLogical}               = true              %--- Use CUDA for the calculations
    options.ShowTheResult {mustBeNumericOrLogical}          = true              %--- Show the result
    options.ShowTheResultDetails {mustBeNumericOrLogical}   = false             %--- Show the magnetization at each L
    options.use_CVODE {mustBeNumericOrLogical}              = false;            %--- Use CVODE for the numerical time evolution
    options.use_minimizer {mustBeNumericOrLogical}          = true;             %--- Relax with the energy minimizer instead of integrating the LL equation in time
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
%% Define the material parameters. These are common to all mesh types
alpha = 1e3;
gamma = 0;
Ms = 1000e3;
K0 = 0.1*1/2*mu0*Ms^2;
A0 = 1.74532925199e-10;

%--- The exchange length and the energy scale of the mumag problem
lex = sqrt(A0/(1/2*mu0*Ms^2));
Km  = 1/2*mu0*Ms^2;

HextFct = @(t) (t)' .* [0,0,0];

%% Load the unstructured Cartesian mesh
%--- The mesh is stored for a single edge length and is scaled to each value in L_loop. The edge
%--- length it was made at is recovered from the total volume of the prisms, which is the cube
if strcmp(options.mesh_type,'unstructuredPrisms')
    data = load(options.mesh_file);
    pos_mesh  = data(:,1:3);
    dims_mesh = data(:,4:6);
    L_mesh    = sum(prod(dims_mesh,2))^(1/3);
    disp(['Loaded ' num2str(size(pos_mesh,1)) ' prisms meshed at L = ' num2str(L_mesh/lex) ' l_ex'])
end

tic
for i = 1:length(L_loop)
    disp(['ITERATION :',num2str(i),'/',num2str(length(L_loop))])
    thisGridL = [lex,lex,lex]*L_loop(i); %m

    %% Setup the mesh
    %--- Each branch creates the problem for its mesh and returns the cell centres, pts, which are
    %--- used below to set up the vortex state
    switch options.mesh_type
        case 'uniform'
            ntot = prod(resolution);
            problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

            [x,y,z] = ndgrid(linspace(-1,1,resolution(1)),linspace(-1,1,resolution(2)),linspace(-1,1,resolution(3)));
            pts = [x(:) y(:) z(:)];

        case 'unstructuredPrisms'
            ntot = size(pos_mesh,1);
            disp(['Prisms N_grid = ' num2str(ntot)])
            problem = DefaultMicroMagProblem(ntot,1,1);
            problem = problem.setMicroMagGridType('unstructuredPrisms');

            %--- Scale the stored mesh to the current edge length
            problem.grid_pts = pos_mesh*L_loop(i)/(L_mesh/lex);
            problem.grid_abc = dims_mesh*L_loop(i)/(L_mesh/lex);
            pts = problem.grid_pts;

        case 'tetrahedron'
            %--- Create the tetrahedral mesh, with a cell size that scales with the cube
            model = CreateTetraMesh(thisGridL, lex*L_loop(i)/options.mesh_res_param);
            ntot = size(model.Mesh.Elements,2);
            disp(['Tetra N_grid = ' num2str(ntot)])
            problem = DefaultMicroMagProblem(ntot,1,1);

            %--- Information on the grid. This is the whole of it: MagTense runs the mesh
            %--- analysis and builds the exchange operator from the mesh
            problem = problem.setMicroMagGridTetrahedron(model.Mesh.Nodes, model.Mesh.Elements);
            pts = problem.grid_pts;
    end
    problem.grid_L = thisGridL;

    %% Setup the problem
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
    problem.Ms = Ms*ones(ntot,1);
    problem.K0 = K0*ones(ntot,1);
    problem.A0 = A0;
    problem.u_ea = zeros( ntot, 3 );
    problem.u_ea(:,3) = 1;

    problem.setTimeDis = int32(10);

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
                %--- A vortex in the xz-plane about the centre of the cube. The centre is taken
                %--- from the extent of the cell centres, so the mesh can sit anywhere
                center = (min(pts,[],1) + max(pts,[],1))/2;
                xvec =  sin(atan2(pts(:,3)-center(3),pts(:,1)-center(1)));
                yvec = -cos(atan2(pts(:,3)-center(3),pts(:,1)-center(1)));
                problem.m0(:,1) = xvec(:) ;
                problem.m0(:,2) = 0.*xvec(:) ;
                problem.m0(:,3) = yvec(:) ;
                problem.m0 = problem.m0./repmat(sqrt(sum(problem.m0.^2,2)),1,3);

                t_end = 200e-9;
        end

        %time grid on which to solve the problem
        problem = problem.setTime( linspace(0,t_end,50) );

        %For the LL relaxation a convergence check at every output time lets the integration stop
        %as soon as the magnetization is stationary instead of running to t_end
        problem = problem.setConvergenceCheckTime( linspace(0,t_end,50) );
        problem.conv_tol = 1e-6;

        %applied field (zero). The minimizer treats every row of the field table as a separate
        %field to relax at, so it gets a single row.
        if options.use_minimizer
            problem = problem.setHext( HextFct, 0 );
        else
            problem = problem.setHext( HextFct, linspace(0,t_end,2) );
        end

        solution = struct();
        prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran

        %--- For an unstructured mesh the second output is the GridInfo that MagTense built from
        %--- the mesh. Its Volumes are what the magnetic moment has to be weighted by
        if strcmp(options.mesh_type,'uniform')
            solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
            tile_volumes = repmat(prod(problem.grid_L)/ntot,ntot,1);
        else
            [solution, GridInfoFortran] = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
            tile_volumes = GridInfoFortran.Volumes;
        end

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
        %--- the solver ran on and the cell volumes of the mesh. Divided by Km*V they are the
        %--- reduced energies of the mumag problem.
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

    %--- Show the unstructured Cartesian mesh once
    if (options.ShowTheResult && i == 1 && strcmp(options.mesh_type,'unstructuredPrisms'))
        cartesianUnstructuredMeshPlot(problem.grid_pts,problem.grid_abc,GridInfoFortran);
    end
end
elapsedTime = toc
disp(['Total field evaluations: ' num2str(sum(n_feval(:)))])

if (options.ShowTheResult)
    plot(fig10,L_loop,sum(E_arr(:,:,1),1),'.','MarkerSize',20)
    plot(fig10,L_loop,sum(E_arr(:,:,2),1),'.','MarkerSize',20)
    xlabel(fig10,'L [l_{ex}]')
    ylabel(fig10,'E [-]')
    legend(fig10,'Flower state','Vortex state');
    title(fig10,['Mesh: ' options.mesh_type])
end

end
