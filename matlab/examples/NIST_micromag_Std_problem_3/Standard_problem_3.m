function [elapsedTime,problem,solution] = Standard_problem_3( resolution, use_CUDA, ShowTheResult, SaveTheResult )

clearvars -except resolution use_CUDA SaveTheResult ShowTheResult
% close all

if ~exist('resolution','var')
    resolution = [5,5,5];
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

if use_CUDA
    dir_m = 'GPU';
else
    dir_m = 'CPU';
end

% if (ShowTheResult)
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
    figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold on; grid on; box on
% end
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

% problem.alpha = 1e3;
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
% AlphaFct = @(t) problem.alpha * 10.^( 5 * min(t,t_end*0.1)/(t_end*0.1) );
% problem = problem.setAlpha( AlphaFct, linspace(0,t_end,100) );
% problem.alpha = 0;

solution = struct();
%convert the class obj to a struct so it can be loaded into fortran
prob_struct = struct(problem);

L_loop = linspace(8,9,2);

tic
for i = 1:length(L_loop)  
    for j = 1%1%2%1:2
        switch j
            %initial magnetization
            case 1
                %Flower state
                problem.m0(:) = 0;
                problem.m0(:,3) = 1 ;
                
                t_end = 200e-9;
                
                %time grid on which to solve the problem
                problem = problem.setTime( linspace(0,t_end,50) );
                
                %time-dependent applied field
                problem = problem.setHext( HextFct, linspace(0,t_end,2) );
            case 2
                %Vortex state
                [x,y,z]=ndgrid(linspace(-1,1,resolution(1)),linspace(-1,1,resolution(2)),linspace(-1,1,resolution(3)));
                xvec =  sin(atan2(z,x));
                yvec = -cos(atan2(z,x));
                xvec(round(resolution(1)/2),:,round(resolution(3)/2)) = 1/sqrt(2);
                yvec(round(resolution(1)/2),:,round(resolution(3)/2)) = -1/sqrt(2);
                problem.m0(:,1) = xvec(:) ;
                problem.m0(:,3) = yvec(:) ;
                % figure; quiver3(x(:),y(:),z(:),problem.m0(:,1),problem.m0(:,2),problem.m0(:,3)); axis equal;
                
                t_end = 200e-9;
                
                %time grid on which to solve the problem
                problem = problem.setTime( linspace(0,t_end,50) );
                
                %time-dependent applied field
                problem = problem.setHext( HextFct, linspace(0,2*t_end,2) );
        end
        
        problem.grid_L = [lex,lex,lex]*L_loop(i);%m

        solution = struct();
        prob_struct = struct(problem);  %convert the class obj to a struct so it can be loaded into fortran

        solution = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );

    %     indx = find(abs(solution.pts(:,3)) < 1e-15 );
    %     figure; M_end = squeeze(solution.M(200,indx,:)); quiver(solution.pts(indx,1),solution.pts(indx,2),M_end(:,1),M_end(:,2)); axis equal;
        figure; M_end = squeeze(solution.M(1,:,:)); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_end(:,1),M_end(:,2),M_end(:,3)); axis equal; title('Fortran starting')
        figure; M_end = squeeze(solution.M(end,:,:)); quiver3(solution.pts(:,1),solution.pts(:,2),solution.pts(:,3),M_end(:,1),M_end(:,2),M_end(:,3)); axis equal; title('Fortran ending')

        plot(fig1,problem.t,mean(solution.M(:,:,1),2),'rx'); 
        plot(fig1,problem.t,mean(solution.M(:,:,2),2),'gx'); 
        plot(fig1,problem.t,mean(solution.M(:,:,3),2),'bx'); 

        %--- Calculate the energy terms
        E_exc = sum((1/2)*(solution.M(:,:,1).*solution.H_exc(:,:,1) + solution.M(:,:,2).*solution.H_exc(:,:,2) + solution.M(:,:,3).*solution.H_exc(:,:,3)),2) ;        
        E_ext = sum(      (solution.M(:,:,1).*solution.H_ext(:,:,1) + solution.M(:,:,2).*solution.H_ext(:,:,2) + solution.M(:,:,3).*solution.H_ext(:,:,3)),2) ;
        E_dem = sum((1/2)*(solution.M(:,:,1).*solution.H_dem(:,:,1) + solution.M(:,:,2).*solution.H_dem(:,:,2) + solution.M(:,:,3).*solution.H_dem(:,:,3)),2) ;
        E_ani = sum((1/2)*(solution.M(:,:,1).*solution.H_ani(:,:,1) + solution.M(:,:,2).*solution.H_ani(:,:,2) + solution.M(:,:,3).*solution.H_ani(:,:,3)),2) ;
        disp(mu0*[E_exc(end) E_ext(end) E_dem(end) E_ani(end)])
%         plot(fig2,problem.t,mu0*E_exc,'.')
%         plot(fig2,problem.t,mu0*E_ext,'.')
%         plot(fig2,problem.t,mu0*E_dem,'.')
%         plot(fig2,problem.t,mu0*E_ani,'.')
        plot(fig2,problem.t,mu0*E_exc-mu0*E_exc(1),'.')
        plot(fig2,problem.t,mu0*E_ext-mu0*E_ext(1),'.')
        plot(fig2,problem.t,mu0*E_dem-mu0*E_dem(1),'.')
        plot(fig2,problem.t,mu0*E_ani-mu0*E_ani(1),'.')
        
        
        %----------- MATLAB -----------------------
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
        plot(fig1, problem.t,Mx,'ro')
        plot(fig1, problem.t,My,'go')
        plot(fig1, problem.t,Mz,'bo')
        figure; quiver3(InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:),SigmaX(:,1),SigmaY(:,1),SigmaZ(:,1)); axis equal;  title('Matlab starting')
        figure; quiver3(InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:),SigmaX(:,end),SigmaY(:,end),SigmaZ(:,end)); axis equal;  title('Matlab ending')
    end
end
elapsedTime = toc

end