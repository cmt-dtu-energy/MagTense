function [elapsedTime,problem,solution,E_arr,L_loop] = Standard_problem_3( resolution, use_CUDA, ShowTheResult, SaveTheResult, ShowTheResultDetails, RunMatlab, L_loop )

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
if ~exist('SaveTheResult','var')
    SaveTheResult = 0;
end
if ~exist('ShowTheResultDetails','var')
    ShowTheResultDetails = 0;
end
if ~exist('RunMatlab','var')
    RunMatlab = 0;
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
    RunMatlab = 1;
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
    figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold on; grid on; box on
    figure3= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig3 = axes('Parent',figure3,'Layer','top','FontSize',16); hold on; grid on; box on
end
mu0 = 4*pi*1e-7;

addpath('../../MEX_files');
addpath('../../util');

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
problem.Ms = 1000e3 ;
problem.K0 = 0.1*1/2*mu0*problem.Ms^2 ;
problem.A0 = 1.74532925199e-10;
problem.u_ea = zeros( prod(resolution), 3 );
problem.u_ea(:,3) = 1;
lex = sqrt(problem.A0/(1/2*mu0*problem.Ms^2));
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
            M_1 = squeeze(solution.M(1,:,:)); figure(figure3); subplot(2,2,1); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_1(:,1),M_1(:,2),M_1(:,3)); axis equal; title('Fortran starting magnetization')
            M_end = squeeze(solution.M(end,:,:)); figure(figure3); subplot(2,2,3); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_end(:,1),M_end(:,2),M_end(:,3)); axis equal; title('Fortran ending magnetization')

            plot(fig1,solution.t,mean(solution.M(:,:,1),2),'rx'); 
            plot(fig1,solution.t,mean(solution.M(:,:,2),2),'gx'); 
            plot(fig1,solution.t,mean(solution.M(:,:,3),2),'bx'); 
        end
        
        %--- Calculate the energy terms
        E_exc = sum((1/2)*(solution.M(:,:,1).*solution.H_exc(:,:,1) + solution.M(:,:,2).*solution.H_exc(:,:,2) + solution.M(:,:,3).*solution.H_exc(:,:,3)),2) ;        
        E_ext = sum(      (solution.M(:,:,1).*solution.H_ext(:,:,1) + solution.M(:,:,2).*solution.H_ext(:,:,2) + solution.M(:,:,3).*solution.H_ext(:,:,3)),2) ;
        E_dem = sum((1/2)*(solution.M(:,:,1).*solution.H_dem(:,:,1) + solution.M(:,:,2).*solution.H_dem(:,:,2) + solution.M(:,:,3).*solution.H_dem(:,:,3)),2) ;
        E_ani = sum((1/2)*(solution.M(:,:,1).*solution.H_ani(:,:,1) + solution.M(:,:,2).*solution.H_ani(:,:,2) + solution.M(:,:,3).*solution.H_ani(:,:,3)),2) ;
        E_arr(:,i,j) = mu0*[E_exc(end) E_ext(end) E_dem(end) E_ani(end)];

        if (ShowTheResultDetails)
            plot(fig2,solution.t,mu0*E_exc-mu0*E_exc(1),'.')
            plot(fig2,solution.t,mu0*E_ext-mu0*E_ext(1),'.')
            plot(fig2,solution.t,mu0*E_dem-mu0*E_dem(1),'.')
            plot(fig2,solution.t,mu0*E_ani-mu0*E_ani(1),'.')
            xlabel(fig2,'Time [s]')
            ylabel(fig2,'Energy [-]')
            legend(fig2,{'E_{exc}','E_{ext}','E_{dem}','E_{ani}'},'Location','East');
        end
        
%% --------------------------------------------------------------------------------------------------------------------------------------
%% -------------------------------------------------------------------- MATLAB ----------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------------------------------------
%% Run the Matlab version of the micromagnetism code
        if (RunMatlab)
            addpath('..\..\micromagnetism')
            problem = problem.setSolverType( 'UseDynamicSolver' );
            problem.DirectoryFilename = [''];
            problem.SimulationName = 'Matlab_DynamicsStdProb3';
            [SigmaSol1,~,InteractionMatrices] = ComputeTheSolution(problem);

            for k=1:size(SigmaSol1,1)
                Sigma = SigmaSol1(k,:).' ;
                NN = round(numel(Sigma)/3) ;
                SigmaX(:,k) = Sigma(0*NN+[1:NN]) ;
                SigmaY(:,k) = Sigma(1*NN+[1:NN]) ;
                SigmaZ(:,k) = Sigma(2*NN+[1:NN]) ;
                SigmaN(:,k) = sqrt(SigmaX(:,k).^2+SigmaY(:,k).^2+SigmaZ(:,k).^2) ;
                Mx(k) = mean(SigmaX(:,k)./SigmaN(:,k)) ;
                My(k) = mean(SigmaY(:,k)./SigmaN(:,k)) ;
                Mz(k) = mean(SigmaZ(:,k)./SigmaN(:,k)) ;
            end
disp(['Max ERROR: ',num2str(max(abs(SigmaN(:)-1)))]) ;
            if (ShowTheResultDetails)
                plot(fig1, problem.t,Mx,'ro')
                plot(fig1, problem.t,My,'go')
                plot(fig1, problem.t,Mz,'bo')
                xlabel(fig1,'Time [s]')
                ylabel(fig1,'Magnetization [-]')
                legend(fig1,{'Fortran M_x','Fortran M_y','Fortran M_z','Matlab M_x','Matlab M_y','Matlab M_z'});

                figure(figure3); subplot(2,2,2); quiver3(InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:),SigmaX(:,1),SigmaY(:,1),SigmaZ(:,1)); axis equal;  title('Matlab starting magnetization')
                figure(figure3); subplot(2,2,4); quiver3(InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:),SigmaX(:,end),SigmaY(:,end),SigmaZ(:,end)); axis equal;  title('Matlab ending magnetization')
                xlabel('x'); ylabel('y'); zlabel('z');
            end
        end
    end
end
elapsedTime = toc

if (ShowTheResult)
    plot(fig10,L_loop,sum(E_arr(:,:,1),1),'.')
    plot(fig10,L_loop,sum(E_arr(:,:,2),1),'.')
    xlabel(fig10,'L [l_{ex}]')
    ylabel(fig10,'E [-]')
end

% disp(['Energy intersection: ' num2str(interp1(sum(E_arr(:,:,1),1)-sum(E_arr(:,:,2),1),L_loop,0))]);

end