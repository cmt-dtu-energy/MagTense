function Script_3D_Std_Problem_4(fig1,MySim_dym)
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
    MySim_ini = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

    MySim_ini.dem_appr = getMicroMagDemagApproximation('none');
    MySim_ini.alpha = 4.42e3;
    MySim_ini.gamma = 0;

    %initial magnetization
    MySim_ini.m0(:) = 1/sqrt(3);

    %time grid on which to solve the problem
    MySim_ini = MySim_ini.setTime( linspace(0,100e-9,200) );
    MySim_ini.setTimeDis = int32(100);
    HystDir = 1/mu0*[1,1,1] ;

    %time-dependent applied field
    HextFct = @(t) (1e-9-t)' .* HystDir .* (t<1e-9)';

    MySim_ini = MySim_ini.setHext( HextFct, linspace(0,100e-9,2000) );
end

if ~(MySim_ini.ShowTheResult == 0)
    if (~exist('fig1','var'))
        figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
    end
    if (~isempty(fig1))
        figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
    end
end

%% Run the simulation
MySim_ini = MySim_ini.setSolverType( 'UseExplicitSolver' );
[SigmaInit,~,InteractionMatrices] = ComputeTheSolution(MySim_ini);

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

if (MySim_ini.ShowTheResult)
    figure; quiver(InteractionMatrices.X(:),InteractionMatrices.Y(:),SigmaX,SigmaY); axis equal;  title('Starting state - Matlab')
end

%% Perform dynamic calculation with first field
MySim_dym = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
MySim_dym.alpha = 4.42e3 ;
MySim_dym.gamma = 2.21e5 ;
MySim_dym = MySim_dym.setUseCuda(MySim_ini.useCuda);
MySim_dym.dem_appr = MySim_ini.dem_appr;
MySim_dym.dem_thres = MySim_ini.dem_thres;
MySim_dym = MySim_dym.setTime( linspace(0,1e-9,200) );
MySim_dym.setTimeDis = int32(10);

if (NIST_field == 1)
    %field 1
    HystDir = 1/mu0*[-24.6,4.3,0]/1000 ;
end
if (NIST_field == 2)
    %field 2
    HystDir = 1/mu0*[-35.5,-6.3,0]/1000 ;
end

HextFct = @(t) (t>-1)' .*HystDir;
MySim_dym = MySim_dym.setHext( HextFct, linspace(0,1e-9,2000) );

MySim_dym.m0(:) = SigmaInit(end,:);

MySim_dym = MySim_dym.setSolverType( 'UseDynamicSolver' );

tic
SigmaSol1 = ComputeTheSolution(MySim_dym);
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

if (MySim_dym.ShowTheResult)
    plot(fig1, MySim_dym.t,Mx,'ro')
    plot(fig1, MySim_dym.t,My,'go')
    plot(fig1, MySim_dym.t,Mz,'bo')
end

end