function Script_3D_Std_Problem_2(fig1,MySim)
%Calculates "implicit" solution (circles) and compares with explicit solution (red line, dots)
disp('Running Matlab model')

%% addpaths
addpath('../MEX_files');
addpath('../util');


if (~exist('fig1','var'))
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
end

clearvars -except fig1 MySim  

%% Setup problem
mu0 = 4*pi*1e-7;

resolution = [20;4;1];
MySim = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

MySim.K0 = 0 ;
MySim.A0 = 1.74532925199e-10;
MySim.Ms = 1000e3 ;
MySim.grid_L = [5e-6,1e-6,1e-7];%m

MySim.alpha = @(t) 1e3*(10.^(5*min(t,2e-9)/2e-9));

MaxH = 0.1;
MySim = MySim.setTime( linspace(MaxH,-0.159*MaxH,101) );

%time-dependent applied field
HystDir = 1/mu0*[1,1,1]/sqrt(3) ;
HextFct = @(t) HystDir .* t';
MySim = MySim.setHext( HextFct, linspace(MaxH,-0.159*MaxH,101) );

ddHextFct = @(t) HystDir + 0.*t';
MySim = MySim.setddHext( ddHextFct, linspace(MaxH,-0.159*MaxH,101) );

t_explicit = 5e-9 ;
MySim = MySim.setTime( linspace(0,t_explicit,2) );

%initial magnetization
MySim.m0(:) = 1/sqrt(3);

%% The implicit solver
MySim = MySim.setSolverType( 'UseImplicitSolver' );
SigmaSol = ComputeTheSolution(MySim);

for k=1:size(SigmaSol,1) 
    Sigma = SigmaSol(k,:).' ;
    NN = round(numel(Sigma)/3) ;
    SigmaX = Sigma(0*NN+[1:NN]) ;
    SigmaY = Sigma(1*NN+[1:NN]) ;
    SigmaZ = Sigma(2*NN+[1:NN]) ;
    SigmaN = sqrt(SigmaX.^2+SigmaY.^2+SigmaZ.^2) ;
    Mx(k) = mean(SigmaX./SigmaN) ;
    My(k) = mean(SigmaY./SigmaN) ;
    Mz(k) = mean(SigmaZ./SigmaN) ;
    Mk(k) = Mx(k)*HystDir(1) + My(k)*HystDir(2) + Mz(k)*HystDir(3) ;
end
plot(fig1,MySim.Hext(:,1),mu0*Mk,'ro') %Minus signs added to correspond to regular hysteresis plots.


%% The explicit solver
MySim = MySim.setSolverType( 'UseExplicitSolver' );

MySim = MySim.setTime( linspace(MaxH,-MaxH,26) );
MySim = MySim.setHext( HextFct, linspace(MaxH,-MaxH,26) );
MySim = MySim.setTimeExplicit( linspace(0,1,numel(MySim.t)) );

MySim = MySim.setTime( linspace(0,t_explicit,2) );

SigmaSol2 = ComputeTheSolution(MySim);

for k=1:size(SigmaSol2,1) 
    Sigma = SigmaSol2(k,:).' ;
    NN = round(numel(Sigma)/3) ;
    SigmaX = Sigma(0*NN+[1:NN]) ;
    SigmaY = Sigma(1*NN+[1:NN]) ;
    SigmaZ = Sigma(2*NN+[1:NN]) ;
    Mx2(k) = mean(SigmaX) ;
    My2(k) = mean(SigmaY) ;
    Mz2(k) = mean(SigmaZ) ;
    Mk2(k) = Mx2(k)*HystDir(1) + My2(k)*HystDir(2) + Mz2(k)*HystDir(3) ;
end
plot(fig1,MySim.Hext(:,1),mu0*Mk2,'r.-')

end