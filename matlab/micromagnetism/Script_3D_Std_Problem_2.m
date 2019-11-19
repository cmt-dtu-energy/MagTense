function Script_3D_Std_Problem_2(fig1,MySim)
%Calculates "implicit" solution (circles) and compares with explicit solution (red line, dots)
disp('Running Matlab model')

if (~exist('fig1','var'))
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
end

clearvars -except fig1 MySim  

%% Setup problem
MySim.SaveTheResult = 0;

MySim.K0 = 0 ;
MySim.Kz = 1 ;
% MySim.A0 = 6.98131700798e-10 ;
MySim.A0 = 1.74532925199e-10;
MySim.Ms = 1000e3 ;
MySim.MaxH = 0.1 ;
MySim.MaxHx = MySim.MaxH*sqrt(3) ;  
MySim.nGrid = [20;4;1] ;

MySim.Lx = 5e-6 ;
MySim.Ly = 1e-6 ;
MySim.Lz = 1e-7 ; 
MySim.FreqH = pi/2 ;
 MySim.MaxT = 5 ;
MySim.MaxT0 = 2 ;
MySim.alpha = @(t) -10000e-10*(10.^(5*min(t,MySim.MaxT0)/MySim.MaxT0));

MySim.HystDir = [1,1,1]./sqrt(3) ;

MySim.ddHsX = MySim.HystDir(1) ;
MySim.ddHsY = MySim.HystDir(2) ;
MySim.ddHsZ = MySim.HystDir(3) ;

MySim.HsX = @(t) MySim.HystDir(1).*t ;
MySim.HsY = @(t) MySim.HystDir(2).*t ;
MySim.HsZ = @(t) MySim.HystDir(3).*t ;
    

MySim.t = linspace(-MySim.MaxH,0.159*MySim.MaxH,101);

MySim.SigmaXfun = @(X,Y,Z) 1/sqrt(3)+0.*X ;
MySim.SigmaYfun = @(X,Y,Z) 1/sqrt(3)+0.*Y ;
MySim.SigmaZfun = @(X,Y,Z) 1/sqrt(3)+0.*Z ;

MySim.InitialState = 'Input' ;

%% The implicit solver
MySim.UseImplicitSolver = 1;
MySim.UseExplicitSolver = 0;

MySim.HeffLimMagn = -1; % The magic number that kinda sorta works

MySim.SimulationName = 'PhysParams_TestDirectHaODEHRdlex60' ;

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
    Mk(k) = Mx(k)*MySim.HystDir(1) + My(k)*MySim.HystDir(2) + Mz(k)*MySim.HystDir(3) ;
end
plot(fig1,-MySim.t,Mk,'ok') %Minus signs added to correspond to regular hysteresis plots.


%% The explicit solver
MySim.UseImplicitSolver = 0;
MySim.UseExplicitSolver = 1;

MySim.MaxT = 2 ;
MySim.HsX = @(t) MySim.HystDir(1).*(-MySim.MaxH +t.*(2*MySim.MaxH./MySim.MaxT)) ;
MySim.HsY = @(t) MySim.HystDir(2).*(-MySim.MaxH +t.*(2*MySim.MaxH./MySim.MaxT)) ;
MySim.HsZ = @(t) MySim.HystDir(3).*(-MySim.MaxH +t.*(2*MySim.MaxH./MySim.MaxT)) ;
MySim.tt = linspace(0,1,26).*MySim.MaxT ;

MySim.CalcEigenvalue = 1 ;

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
    Mk2(k) = Mx2(k)*MySim.HystDir(1) + My2(k)*MySim.HystDir(2) + Mz2(k)*MySim.HystDir(3) ;
end
plot(fig1,-linspace(-MySim.MaxH,MySim.MaxH,numel(Mx2)),Mk2,'.-r')
legend('"Implicit method"', 'Explicit method','Location','SouthEast')
ylabel('M')
xlabel('H_{applied} [T]')
end