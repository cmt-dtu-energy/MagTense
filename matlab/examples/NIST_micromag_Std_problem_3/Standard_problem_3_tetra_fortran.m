function [elapsedTime,problem,solution,E_arr,L_loop] = Standard_problem_3_tetra_fortran( resolution, use_CUDA, ShowTheResult, SaveTheResult, ShowTheResultDetails, RunMatlab, L_loop )

clearvars -except  resolution use_CUDA ShowTheResult SaveTheResult ShowTheResultDetails RunMatlab L_loop 
close all

if ~exist('resolution','var')
    resolution = [9,9,9];
end
if ~exist('use_CUDA','var')
    use_CUDA = true;
end
if ~exist('ShowTheResult','var')
    ShowTheResult = 1;
end
if ~exist('SaveTheResult','var')
    SaveTheResult = 0;
end
if ~exist('ShowTheResultDetails','var')
    ShowTheResultDetails = 1;
end
if ~exist('RunMatlab','var')
    RunMatlab = 1;
end
if ~exist('L_loop','var')
    L_loop = 8.5;%linspace(8,9,10);
end

if use_CUDA
    dir_m = 'GPU';
else
    dir_m = 'CPU';
end

% if (ShowTheResult)
%     figure10= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig10 = axes('Parent',figure10,'Layer','top','FontSize',16); hold on; grid on; box on
% end

if (ShowTheResultDetails)
    RunMatlab = 1;
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
%     figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold on; grid on; box on
    figure3= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig3 = axes('Parent',figure3,'Layer','top','FontSize',16); hold on; grid on; box on
end
mu0 = 4*pi*1e-7;

addpath('../../MEX_files');
addpath('../../util');
addpath('../../micromagnetism');

%% --------------------------------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------- MAGTENSE ---------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Setup the problem for the initial configuration
%takes the size of the grid as arguments (nx,ny,nz) and a function handle
%that produces the desired field (if not present zero applied field is inferred)
load('InteractionMatrices_test.mat') ; % Voronoi & Tetra
resolution = [length(model.Mesh.Elements(:,:)) 1 1];

problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
problem.dem_appr = getMicroMagDemagApproximation('none');
problem = problem.setUseCuda( use_CUDA );
problem.grid_type = getMicroMagGridType('tetrahedron');

GridInfo = TetrahedralMeshAnalysis(model) ;
problem.grid_pts = [GridInfo.Xel, GridInfo.Yel, GridInfo.Zel] ;
problem.grid_ele = int32(model.Mesh.Elements(:,:)) ;
problem.grid_nod = model.Mesh.Nodes(:,:) ;      
problem.grid_n = [int32(length(problem.grid_ele)) int32(1) int32(1)];
problem.grid_nnod = int32(length(problem.grid_nod));
% InteractionMatrices.A2 = 0*InteractionMatrices.A2;
% problem.exch_mat = InteractionMatrices.A2;
[v,c,rs,re] = ConvertToCSR(InteractionMatrices.A2);
problem.exch_nval = int32(numel(v));
problem.exch_nrow = int32(numel(rs));
problem.exch_val  = single(v);
problem.exch_rows = int32(rs);
problem.exch_rowe = int32(re);
problem.exch_col  = int32(c);

problem.gamma = 0;
problem.Ms = 1000e3 ;
problem.K0 = 0.1*1/2*mu0*problem.Ms^2 ;
problem.A0 = 1.74532925199e-10;
problem.u_ea = zeros( prod(problem.grid_n), 3 );
problem.u_ea(:,3) = 1;
lex = sqrt(problem.A0/(1/2*mu0*problem.Ms^2));
problem.setTimeDis = int32(10);
HextFct = @(t) (t)' .* [0,0,0];

%time-dependent alpha parameter, to ensure faster convergence
problem.alpha = 1e3;

tic
for i = 1:length(L_loop)  
    disp(['ITERATION :',num2str(i),'/',num2str(length(L_loop))])
    for j = 2%1:2
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
                xvec(round(resolution(1)/2),:,round(resolution(3)/2)) = 1/sqrt(2);
                yvec(round(resolution(1)/2),:,round(resolution(3)/2)) = -1/sqrt(2);
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

        solution = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

        if (ShowTheResultDetails)
            M_1 = squeeze(solution.M(1,:,:)); figure(figure3); subplot(2,2,2); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_1(:,1),M_1(:,2),M_1(:,3)); axis equal; title('Fortran starting magnetization')
            M_end = squeeze(solution.M(end,:,:)); figure(figure3); subplot(2,2,4); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_end(:,1),M_end(:,2),M_end(:,3)); axis equal; title('Fortran ending magnetization')

            plot(fig1,solution.t,mean(solution.M(:,:,1),2),'rd'); 
            plot(fig1,solution.t,mean(solution.M(:,:,2),2),'gd'); 
            plot(fig1,solution.t,mean(solution.M(:,:,3),2),'bd'); 
        end
        
        %--- Calculate the energy terms
        E_exc = sum((1/2)*(solution.M(:,:,1).*solution.H_exc(:,:,1) + solution.M(:,:,2).*solution.H_exc(:,:,2) + solution.M(:,:,3).*solution.H_exc(:,:,3)),2) ;        
        E_ext = sum(      (solution.M(:,:,1).*solution.H_ext(:,:,1) + solution.M(:,:,2).*solution.H_ext(:,:,2) + solution.M(:,:,3).*solution.H_ext(:,:,3)),2) ;
        E_dem = sum((1/2)*(solution.M(:,:,1).*solution.H_dem(:,:,1) + solution.M(:,:,2).*solution.H_dem(:,:,2) + solution.M(:,:,3).*solution.H_dem(:,:,3)),2) ;
        E_ani = sum((1/2)*(solution.M(:,:,1).*solution.H_ani(:,:,1) + solution.M(:,:,2).*solution.H_ani(:,:,2) + solution.M(:,:,3).*solution.H_ani(:,:,3)),2) ;
        E_arr(:,i,j) = mu0*[E_exc(end) E_ext(end) E_dem(end) E_ani(end)];
        
%% --------------------------------------------------------------------------------------------------------------------------------------
%% -------------------------------------------------------------------- MATLAB ----------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Run the Matlab version of the micromagnetism code
        if (RunMatlab)
            problem = problem.setUseCuda( false );
            addpath('..\..\micromagnetism')
            problem = problem.setSolverType( 'UseDynamicSolver' );
            problem.DirectoryFilename = [''];
            problem.SimulationName = 'Matlab_DynamicsStdProb3';

            thisGridL = [lex,lex,lex]*L_loop(i);%m
            resolution = [size(model.Mesh.Elements,2),1,1];

            %% Setup the problem for the initial configuration
            %takes the size of the grid as arguments (nx,ny,nz) and a function handle
            %that produces the desired field (if not present zero applied field is inferred)
            problem_tetra = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
            problem_tetra.grid_L = thisGridL;% NOT STANDARD
            problem_tetra.RecomputeInteractionMatrices = 1 ;
            problem_tetra.ExternalMesh = 1 ;
            problem_tetra.MeshType = 'Tetra' ;
            problem_tetra.ExternalMeshFileName = 'InteractionMatrices_test.mat' ;

            problem_tetra.dem_appr = problem.dem_appr;
            problem_tetra = problem_tetra.setUseCuda( false );
            problem_tetra.setTimeDis = int32(10);

            problem_tetra.gamma = problem.gamma;
            problem_tetra.Ms = problem.Ms;
            problem_tetra.K0 = problem.K0;
            problem_tetra.A0 = problem.A0;
            problem_tetra.alpha = problem.alpha;
            problem_tetra.m0 = problem.m0;

            problem_tetra.u_ea = problem.u_ea;
            problem_tetra = problem_tetra.setTime( linspace(0,t_end,50) );

            %time-dependent applied field
            problem_tetra = problem_tetra.setHext( HextFct, linspace(0,t_end,2) );
            problem_tetra = problem_tetra.setSolverType( 'UseDynamicSolver' );
            problem_tetra.DirectoryFilename = [''];
            problem_tetra.SimulationName = 'Matlab_DynamicsStdProb3';
            disp(['Tetra N_grid = ' num2str(prod(resolution))])
            [SigmaSol1,~,InteractionMatrices] = ComputeTheSolution(problem_tetra);
            [Mx,My,Mz] = ComputeMagneticMomentGeneralMesh(SigmaSol1,InteractionMatrices.GridInfo.Volumes) ;

            plot(fig1, problem_tetra.t,Mx,'rx')
            plot(fig1, problem_tetra.t,My,'gx')
            plot(fig1, problem_tetra.t,Mz,'bx')
            h_l =legend(fig1,{'Uniform M_x','Uniform M_y','Uniform M_z','Tetra M_x','Tetra M_y','Tetra M_z'},'Location','East');
            set(h_l,'Fontsize',8)

            for k=1:size(SigmaSol1,1)
                Sigma = SigmaSol1(k,:).' ;
                NN = round(numel(Sigma)/3) ;
                SigmaX_tetra(:,k) = Sigma(0*NN+[1:NN]) ;
                SigmaY_tetra(:,k) = Sigma(1*NN+[1:NN]) ;
                SigmaZ_tetra(:,k) = Sigma(2*NN+[1:NN]) ;
            end

            figure(figure3); subplot(2,2,1); quiver3(InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:),SigmaX_tetra(:,1),SigmaY_tetra(:,1),SigmaZ_tetra(:,1)); axis equal;  title('Matlab starting magnetization')
            figure(figure3); subplot(2,2,3); quiver3(InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:),SigmaX_tetra(:,end),SigmaY_tetra(:,end),SigmaZ_tetra(:,end)); axis equal;  title('Matlab ending magnetization')
            xlabel('x'); ylabel('y'); zlabel('z');

        end
    end
end
elapsedTime = toc

end