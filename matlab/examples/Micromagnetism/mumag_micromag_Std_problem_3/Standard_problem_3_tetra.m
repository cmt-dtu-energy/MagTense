function [elapsedTime,problem,solution,E_arr_v,L_loop] = Standard_problem_3_tetra( options )
%STANDARD_PROBLEM_3_TETRA
% Mumag standard problem 3 on a tetrahedral mesh.
%
% The cube is meshed with tetrahedra and the flower and the vortex state are relaxed at a
% sequence of edge lengths, in units of the exchange length. The energy of the two states
% crosses at the edge length that the standard problem asks for.
%
% The mesh is the only thing given to MagTense: setMicroMagGridTetrahedron passes the nodes and
% the connectivity, and MagTense analyses the mesh and assembles the exchange operator itself,
% exactly as Standard_problem_3_unstructured_cart does by giving the centres and sizes of a grid
% of unstructured prisms. Nothing has to be computed in Matlab beforehand.
%
% This replaces the older approach, in which the mesh analysis and the differential operators
% were run in Matlab and the assembled exchange matrix was handed to MagTense through
% setExchangeMatrixSparse. That route still works and is still the way to use an exchange
% operator that MagTense cannot build itself, but it is no longer needed for an ordinary
% tetrahedral problem.
%
% Requires the PDE Toolbox, for generateMesh in CreateTetraMesh.

arguments
    options.use_CUDA {mustBeNumericOrLogical}             = true    %--- Use CUDA for the calculations
    options.ShowTheResult {mustBeNumericOrLogical}        = true    %--- Show the energy crossing
    options.ShowTheResultDetails {mustBeNumericOrLogical} = true   %--- Show the magnetisation at each L
    options.use_CVODE {mustBeNumericOrLogical}            = false   %--- Use CVODE for the time evolution
    options.L_loop                                        = linspace(8,9,10)  %--- Edge lengths [l_ex]
    options.mesh_res_param                                = 5       %--- Cells per edge length
end

L_loop = options.L_loop;

if (options.ShowTheResultDetails)
    figure1 = figure('PaperType','A4','Visible','on','PaperPositionMode','auto');
    fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
    figure3 = figure('PaperType','A4','Visible','on','PaperPositionMode','auto');
end

mu0 = 4*pi*1e-7;

addpath('../../../MEX_files');
addpath('../../../util');

%% --------------------------------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------- MAGTENSE ---------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
tic
for i = 1:length(L_loop)
    disp(['ITERATION :',num2str(i),'/',num2str(length(L_loop))])

    %--- Define parameters
    alpha = 1e3;
    gamma = 0;
    Ms = 1000e3;
    K0 = 0.1*1/2*mu0*Ms^2;
    A0 = 1.74532925199e-10;
    lex = sqrt(A0/(1/2*mu0*Ms^2));
    thisGridL = [lex,lex,lex]*L_loop(i);    %m
    mesh_res = lex*L_loop(i)/options.mesh_res_param;

    %--- Create the tetrahedral mesh
    model = CreateTetraMesh(thisGridL,mesh_res);

    %--- Setup the problem
    resolution = [size(model.Mesh.Elements,2) 1 1];
    disp(['Tetra N_grid = ' num2str(prod(resolution))])
    problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
    problem = problem.setUseCuda( options.use_CUDA );
    problem = problem.setUseCVODE( options.use_CVODE );
    problem = problem.setMicroMagDemagApproximation('none');
    problem.ReturnHall = int32(1);

    %--- Information on the grid. This is the whole of it: MagTense runs the mesh analysis and
    %--- builds the exchange operator from the mesh, so no exchange matrix is passed in
    problem = problem.setMicroMagGridTetrahedron(model.Mesh.Nodes, model.Mesh.Elements);

    %--- Save the parameters
    problem.alpha = alpha;
    problem.gamma = gamma;
    problem.Ms = Ms;
    problem.K0 = K0;
    problem.A0 = A0;
    problem.u_ea = zeros( prod(resolution), 3 );
    problem.u_ea(:,3) = 1;

    problem.setTimeDis = int32(10);
    HextFct = @(t) (t)' .* [0,0,0];

    for j = 1:2
        switch j
            case 1
                disp('Flower state');
                problem.m0(:) = 0;
                problem.m0(:,3) = 1;
                t_end = 10e-9;

            case 2
                disp('Vortex state');
                xvec =  sin(atan2(problem.grid_pts(:,3),problem.grid_pts(:,1)));
                yvec = -cos(atan2(problem.grid_pts(:,3),problem.grid_pts(:,1)));
                problem.m0(:,1) = xvec(:);
                problem.m0(:,2) = 0.*xvec(:);
                problem.m0(:,3) = yvec(:);
                problem.m0 = problem.m0./repmat(sqrt(sum(problem.m0.^2,2)),1,3);
                t_end = 200e-9;
        end

        %--- Time grid on which to solve the problem
        problem = problem.setTime( linspace(0,t_end,50) );

        %--- Time-dependent applied field
        problem = problem.setHext( HextFct, linspace(0,t_end,2) );

        problem.grid_L = thisGridL;

        problem.nThreads = int32(2);    %Threads used by OpenMP for building the demag tensor

        solution = struct();
        prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran

        %--- The second output is the GridInfo that MagTense built from the mesh. Its Volumes are
        %--- what the magnetic moment of an unstructured mesh has to be weighted by
        [solution, GridInfo] = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

        [Mx,My,Mz,mx,my,mz] = computeMagneticMomentGeneralMesh(solution.M,GridInfo.Volumes);

        if (options.ShowTheResultDetails)
            M_1 = squeeze(solution.M(1,:,:));
            figure(figure3); subplot(2,2,2);
            quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_1(:,1),M_1(:,2),M_1(:,3));
            axis equal; title('Starting magnetization')
            M_end = squeeze(solution.M(end,:,:));
            figure(figure3); subplot(2,2,4);
            quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_end(:,1),M_end(:,2),M_end(:,3));
            axis equal; title('Ending magnetization')

            plot(fig1,solution.t,Mx,'rd');
            plot(fig1,solution.t,My,'gd');
            plot(fig1,solution.t,Mz,'bd');
        end

        %--- Calculate the energy terms
        E_exc = sum((1/2)*(mx(:,:).*solution.H_exc(:,:,1) + my(:,:).*solution.H_exc(:,:,2) + mz(:,:).*solution.H_exc(:,:,3)),2);
        E_ext = sum(      (mx(:,:).*solution.H_ext(:,:,1) + my(:,:).*solution.H_ext(:,:,2) + mz(:,:).*solution.H_ext(:,:,3)),2);
        E_dem = sum((1/2)*(mx(:,:).*solution.H_dem(:,:,1) + my(:,:).*solution.H_dem(:,:,2) + mz(:,:).*solution.H_dem(:,:,3)),2);
        E_ani = sum((1/2)*(mx(:,:).*solution.H_ani(:,:,1) + my(:,:).*solution.H_ani(:,:,2) + mz(:,:).*solution.H_ani(:,:,3)),2);
        E_arr_v(:,i,j) = mu0*[E_exc(end) E_ext(end) E_dem(end) E_ani(end)];
    end
end
elapsedTime = toc

if (options.ShowTheResult)
    figure10 = figure('PaperType','A4','Visible','on','PaperPositionMode','auto');
    fig10 = axes('Parent',figure10,'Layer','top','FontSize',16); hold on; grid on; box on
    plot(fig10,L_loop,sum(E_arr_v(:,:,1),1)/sum(E_arr_v(:,1,1),1),'rd')
    plot(fig10,L_loop,sum(E_arr_v(:,:,2),1)/sum(E_arr_v(:,1,1),1),'bd')
    xlabel(fig10,'L [l_{ex}]')
    ylabel(fig10,'E [a.u.]')
    legend(fig10,'Flower state','Vortex state');
end

end
