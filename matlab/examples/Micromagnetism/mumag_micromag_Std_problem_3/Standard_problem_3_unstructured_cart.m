function [elapsedTime,problem,solution,E_arr] = Standard_problem_3_unstructured_cart( options )

arguments
    options.use_CUDA {mustBeNumericOrLogical}               = true              %--- Use CUDA for the calculations
    options.ShowTheResult {mustBeNumericOrLogical}          = true              %--- Show the result
    options.ShowTheResultDetails {mustBeNumericOrLogical}   = false             %--- Show the magnetization at each L
    options.use_CVODE {mustBeNumericOrLogical}              = false;            %--- Use CVODE for the numerical time evolution
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
%--- Load the unstructured mesh
%--- The parameters Ms, K0 and A0 are used to generate the mesh and thus they are loaded as well
load('Std_prob_3_unstructured_cartesian_grains_9_mesh_4_ref_2.mat')
% load('Std_prob_3_unstructured_cartesian_grains_9_mesh_7_ref_2.mat')

%% Setup the problem for the initial configuration
for i = 1:length(L_loop)  
    disp(['ITERATION :',num2str(i),'/',num2str(length(L_loop))])
    
    mesh = mesh_arr(i);
    GridInfo = GridInfo_arr(i);
    if (i == 1)
        cartesianUnstructuredMeshPlot(mesh.pos_out,mesh.dims_out,GridInfo,mesh.iIn);
    end

    %--- Define parameters
    alpha = 1e3;
    gamma = 0;

    %--- Setup the problem
    resolution = [length(mesh.pos_out) 1 1];
    disp(['Prisms N_grid = ' num2str(prod(resolution))])
    problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
    problem = problem.setUseCuda( options.use_CUDA );
    problem = problem.setUseCVODE( options.use_CVODE );
    problem = problem.setMicroMagGridType('unstructuredPrisms');

    %--- Save the parameters
    problem.alpha = alpha;
    problem.gamma = gamma;
    problem.Ms = Ms;
    problem.K0 = K0;
    problem.A0 = A0;
    problem.u_ea = zeros( prod(resolution), 3 );
    problem.u_ea(:,3) = 1;

    %--- Information on the grid
    problem.grid_pts    = mesh.pos_out;
    problem.grid_abc    = mesh.dims_out;
    lex = sqrt(problem.A0/(1/2*mu0*Ms^2));
    thisGridL = [lex,lex,lex]*L_loop(i); %m

    %--- Calculate the exchange matrix
    InteractionMatrices.GridInfo = GridInfo;
    InteractionMatrices.X = GridInfo.Xel ;
    InteractionMatrices.Y = GridInfo.Yel ;
    InteractionMatrices.Z = GridInfo.Zel ;

    [D2X,D2Y,D2Z] = computeDifferentialOperatorsFromMesh_DirectLap(GridInfo,'extended',6,"DirectLaplacianNeumann");
    InteractionMatrices.A2 = D2X + D2Y + D2Z ;

    %--- Convert the exchange matrix to sparse
    problem = problem.setExchangeMatrixSparse( InteractionMatrices.A2 );

    problem.setTimeDis = int32(10);
    HextFct = @(t) (t)' .* [0,0,0];

    tic

    for j = 1:2
        switch j
            %initial magnetization
            case 1
                disp('Flower state') ;
                problem.m0(:) = 0;
                problem.m0(:,3) = 1 ;
                
                t_end = 300e-9;  
                
            case 2
                disp('Vortex state') ;
                xvec =  sin(atan2(problem.grid_pts(:,3)-thisGridL(3)/2,problem.grid_pts(:,1)-thisGridL(1)/2));
                yvec = -cos(atan2(problem.grid_pts(:,3)-thisGridL(3)/2,problem.grid_pts(:,1)-thisGridL(1)/2));
                problem.m0(:,1) = xvec(:) ;
                problem.m0(:,2) = 0.*xvec(:) ;
                problem.m0(:,3) = yvec(:) ;
                problem.m0 = problem.m0./repmat(sqrt(sum(problem.m0.^2,2)),1,3);
                
                t_end = 300e-9;
        end

        %time grid on which to solve the problem
        problem = problem.setTime( linspace(0,t_end,50) );

        %time-dependent applied field
        problem = problem.setHext( HextFct, linspace(0,t_end,2) );
        
        problem.grid_L = [lex,lex,lex]*L_loop(i);%m

        solution = struct();
        prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran

        solution = problem.MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
        
        [Mx,My,Mz,mx,my,mz] = computeMagneticMomentGeneralMesh(solution.M,GridInfo.Volumes) ;

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
        
        toc

        %--- Calculate the energy terms            
        [E_dem_red, E_exc_red, E_ani_red, E_ext_red] = computeMagneticEnergy(solution,GridInfo.Volumes,problem,Ms);
        E_arr(:,i,j) = [E_dem_red(end) E_exc_red(end) E_ani_red(end) E_ext_red(end)];
    end
end
elapsedTime = toc

if (options.ShowTheResult)
    plot(fig10,L_loop,sum(E_arr(:,:,1),1),'.')
    plot(fig10,L_loop,sum(E_arr(:,:,2),1),'.')
    xlabel(fig10,'L [l_{ex}]')
    ylabel(fig10,'E [-]')
    legend(fig10,'Flower state','Vortex state');
end

end