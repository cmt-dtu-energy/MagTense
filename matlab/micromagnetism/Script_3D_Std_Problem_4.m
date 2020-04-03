function Script_3D_Std_Problem_4(fig1,MySim)
disp('Running Matlab model')

%% Calculates dynamic solution to mumag std prob 4
clearvars -except fig1 MySim
mu0 = 4*pi*1e-7;

%% addpaths
addpath('../MEX_files');
addpath('../util');

%% Create Preset
if (~exist('MySim','var')) 
    NIST_field = 1;
    
    resolution = [1*36,1*9,1];
    MySim = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

    MySim.dem_appr = getMicroMagDemagApproximation('none');
    MySim.alpha = 4.42e3;
    MySim.gamma = 0;

    %initial magnetization
    MySim.m0(:) = 1/sqrt(3);

    %time grid on which to solve the problem
    MySim = MySim.setTime( linspace(0,100e-9,200) );
    MySim.setTimeDis = int32(100);
    HystDir = 1/mu0*[1,1,1] ;

    %time-dependent applied field
    HextFct = @(t) (1e-9-t)' .* HystDir .* (t<1e-9)';

    MySim = MySim.setHext( HextFct, 10*numel(t) );
end

if ~(MySim.ShowTheResult == 0)
    if (~exist('fig1','var'))
        figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
    end
    if (~isempty(fig1))
        figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
    end
end

%% Run the simulation
MySim = MySim.setSolverType( 'UseExplicitSolver' );
[SigmaInit,~,~,InteractionMatrices] = ComputeTheSolution(MySim);

for k=1:size(SigmaInit,1) 
    Sigma = SigmaInit(k,:).' ;
    NN = round(numel(Sigma)/3) ;

    SigmaX = Sigma(0*NN+[1:NN]) ;
    SigmaY = Sigma(1*NN+[1:NN]) ;
    SigmaZ = Sigma(2*NN+[1:NN]) ;
    SigmaN = sqrt(SigmaX.^2+SigmaY.^2+SigmaZ.^2) ;
    Mx(k) = mean(SigmaX./SigmaN) ;
    My(k) = mean(SigmaY./SigmaN) ;
    Mz(k) = mean(SigmaZ./SigmaN) ;
end

if (MySim.ShowTheResult)
    figure; quiver(InteractionMatrices.X(:),InteractionMatrices.Y(:),SigmaX,SigmaY); axis equal;  title('Starting state - Matlab')
end

%% Perform dynamic calculation with first field
MySim = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
MySim.alpha = 4.42e3 ;
MySim.gamma = 2.21e5 ;
MySim = MySim.setUseCuda( 0 );
MySim.dem_appr = getMicroMagDemagApproximation('none');
MySim.dem_thres = 1e-4;
MySim = MySim.setTime( linspace(0,1e-9,200) );
MySim.setTimeDis = int32(10);

if (NIST_field == 1)
    %field 1
    HystDir = 1/mu0*[-24.6,4.3,0]/1000 ;
end
if (NIST_field == 2)
    %field 2
    HystDir = 1/mu0*[-35.5,-6.3,0]/1000 ;
end

HextFct = @(t) (t>-1)' .*HystDir;
MySim = MySim.setHext( HextFct );

MySim.m0(:) = SigmaInit(end,:);

MySim = MySim.setSolverType( 'UseDynamicSolver' );

tic
SigmaSol1 = ComputeTheSolution(MySim);
toc

for k=1:size(SigmaSol1,1) 
    Sigma = SigmaSol1(k,:).' ;
    NN = round(numel(Sigma)/3) ;

    SigmaX = Sigma(0*NN+[1:NN]) ;
    SigmaY = Sigma(1*NN+[1:NN]) ;
    SigmaZ = Sigma(2*NN+[1:NN]) ;
    SigmaN = sqrt(SigmaX.^2+SigmaY.^2+SigmaZ.^2) ;
    Mx(k) = mean(SigmaX./SigmaN) ;
    My(k) = mean(SigmaY./SigmaN) ;
    Mz(k) = mean(SigmaZ./SigmaN) ;
end

if (MySim.ShowTheResult)
    plot(fig1, MySim.t,Mx,'ro')
    plot(fig1, MySim.t,My,'go')
    plot(fig1, MySim.t,Mz,'bo')
end

end