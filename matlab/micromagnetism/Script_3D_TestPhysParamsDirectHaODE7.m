%Calculates "implicit" solution (circles) and compares with explicit
%solution (red line, dots)

clear 
% %% Create Preset
a = memory
MySim.K0 = 0 ;
MySim.Kz = 1 ;
% MySim.A0 = 6.98131700798e-10 ;
MySim.A0 = 1.74532925199e-10;
MySim.Ms = 1000e3 ;
MySim.MaxH = 0.1 ;
MySim.MaxHx = MySim.MaxH*sqrt(3) ;  
%nHalfTimes = [1,2;1,2;0,1] ;
%nHalfTimes = [0,19;0,4;0,1] ;
% nHalfTimes = [19;4;1] ;
% nHalfTimes = [6;3;0] ;
% nHalfTimes = [19;4;1] ; %OLD GRID
MySim.nHalfTimes = [19;4;1] ;
% MySim.nHalfTimes = [51;13;2] ;
% MySim.nHalfTimes = [40;9;0] ;
% MySim.nHalfTimes = [30;6;0] ;

% nHalfTimes = [6;2;0] ;

MySim.Lx = 5e-6 ;
MySim.Ly = 1e-6 ;
MySim.Lz = 1e-7 ; 
% nHalfTimes = [12] ;
MySim.FreqH = pi/2 ;
% MySim.MaxT = 3.0276 ;
 MySim.MaxT = 5 ;
MySim.MaxT0 = 2 ;
% MySim.MaxT1 = 3 ; 
% MySim.alpha = @(t) -10000e9*(10.^(5*min(t,MySim.MaxT0)/MySim.MaxT0));
MySim.alpha = @(t) -10000e-10*(10.^(5*min(t,MySim.MaxT0)/MySim.MaxT0));

% HsX = @(t) -MaxH.*(1-2.*t./MaxT) ;

MySim.HystDir = [1,1,1]./sqrt(3) ;
% HystDir = [1,0,0] ;

MySim.ddHsX = MySim.HystDir(1) ;
MySim.ddHsY = MySim.HystDir(2) ;
MySim.ddHsZ = MySim.HystDir(3) ;

MySim.HsX = @(t) MySim.HystDir(1).*t ;
MySim.HsY = @(t) MySim.HystDir(2).*t ;
MySim.HsZ = @(t) MySim.HystDir(3).*t ;
    

MySim.t = linspace(-MySim.MaxH,0.159*MySim.MaxH,101)
% MySim.t = linspace(-MySim.MaxH,MySim.MaxH,101);    
% t = linspace(-MaxH,MaxH,101) ; 
% MySim.t = MySim.t(82:end);
%      nHalfTimes = [2,2,3] ;
% alpha = @(t) -1000e17 ;
MySim.GifFilename = 'ThisPhysTest.gif' ;

MySim.SigmaXfun = @(X,Y,Z) 1/sqrt(3)+0.*X ;
MySim.SigmaYfun = @(X,Y,Z) 1/sqrt(3)+0.*Y ;
MySim.SigmaZfun = @(X,Y,Z) 1/sqrt(3)+0.*Z ;

MySim.InitialState = 'Input' ;

MySim.UseImplicitSolver = 1;
MySim.UseExplicitSolver = 0;

MySim.HeffLimMagn = -1; % The magic number that kinda sorta works

MySim.SimulationName = 'PhysParams_TestDirectHaODEHRdlex60' ;

% [ProblemSetupStruct,ThisPreSetFileName] = SetupProblem(ThisPreSetFileName,t,HystDir,MaxT,K0,Kz,A0,Ms,MaxH,MaxHx,HsX,HsY,HsZ,ddHsX,ddHsY,ddHsZ,nHalfTimes,GifFilename,InitialState,SigmaXfun,SigmaYfun,SigmaZfun,Lx,Ly,Lz);
SigmaSol = ComputeTheSolution(MySim);
% [ProblemSetupStruct,ThisPreSetFileName] = SetupProblem(MySim);

% save(ThisPreSetFileName,'t','HystDir','MaxT','K0','Kz','A0','Ms','MaxH','MaxHx','HsX','HsY','HsZ','ddHsX','ddHsY','ddHsZ','nHalfTimes','GifFilename','InitialState','SigmaXfun','SigmaYfun','SigmaZfun','Lx','Ly','Lz') ;
% ThisSigma = LandauLifshitzEvolveDirectHaODEinit(ThisPreSetFileName) ;
% ThisSigma = LandauLifshitzEvolveDirectHaODEinit(ProblemSetupStruct);
% return 

% load([SaveFileName(1:end-4),'solution'])
for k=1:size(SigmaSol,1) 
Sigma = SigmaSol(k,:).' ;
NN = round(numel(Sigma)/3) ;
% N = round(NN^(1/3)) ;
SigmaX = Sigma(0*NN+[1:NN]) ;
SigmaY = Sigma(1*NN+[1:NN]) ;
SigmaZ = Sigma(2*NN+[1:NN]) ;
SigmaN = sqrt(SigmaX.^2+SigmaY.^2+SigmaZ.^2) ;
Mx(k) = mean(SigmaX./SigmaN) ;
My(k) = mean(SigmaY./SigmaN) ;
Mz(k) = mean(SigmaZ./SigmaN) ;
Mk(k) = Mx(k)*MySim.HystDir(1) + My(k)*MySim.HystDir(2) + Mz(k)*MySim.HystDir(3) ;
end
% 
figure
plot(MySim.t,Mk)
'' ;
%%
MySim.UseImplicitSolver = 0;
MySim.UseExplicitSolver = 1;

MySim.MaxT = 2 ;
MySim.HsX = @(t) MySim.HystDir(1).*(-MySim.MaxH +t.*(2*MySim.MaxH./MySim.MaxT)) ;
MySim.HsY = @(t) MySim.HystDir(2).*(-MySim.MaxH +t.*(2*MySim.MaxH./MySim.MaxT)) ;
MySim.HsZ = @(t) MySim.HystDir(3).*(-MySim.MaxH +t.*(2*MySim.MaxH./MySim.MaxT)) ;
MySim.tt = linspace(0,1,26).*MySim.MaxT ;
% ThisPreSetFileNameCOMP = 'PhysParams_TestDirectCompare05.mat' ;
MySim.SimulationName = 'PhysParams_TestDirectHRdlex60' ;
MySim.alpha = @(t) -10000e-10*(10.^(5*min(t,MySim.MaxT0)/MySim.MaxT0)); 
MySim.CalcEigenvalue = 1 ;
% save(ThisPreSetFileNameCOMP,'CalcEigs','alpha','HystDir','MaxT','K0','Kz','A0','Ms','MaxH','MaxHx','HsX','HsY','HsZ','ddHsX','ddHsY','ddHsZ','nHalfTimes','GifFilename','InitialState','SigmaXfun','SigmaYfun','SigmaZfun','Lx','Ly','Lz','tt') ;

%  LandauLifshitzEvolveStepByStep(ThisPreSetFileNameCOMP) ;
SigmaSol2 = ComputeTheSolution(MySim);
 %%
 %load([ThisPreSetFileNameCOMP(1:end-4),'solution'])
 clear Mx2
for k=1:size(SigmaSol2,1) 
Sigma = SigmaSol2(k,:).' ;
NN = round(numel(Sigma)/3) ;
% N = round(NN^(1/3)) ;
SigmaX = Sigma(0*NN+[1:NN]) ;
SigmaY = Sigma(1*NN+[1:NN]) ;
SigmaZ = Sigma(2*NN+[1:NN]) ;
Mx2(k) = mean(SigmaX) ;
My2(k) = mean(SigmaY) ;
Mz2(k) = mean(SigmaZ) ;
Mk2(k) = Mx2(k)*MySim.HystDir(1) + My2(k)*MySim.HystDir(2) + Mz2(k)*MySim.HystDir(3) ;
end
%% Plot the solutions. Minus signs added to correspond to regular hysteresis plots.
figure
plot(-MySim.t,Mk,'ok') 
hold on
plot(-linspace(-MySim.MaxH,MySim.MaxH,numel(Mx2)),Mk2,'.-r')
legend('"Implicit method"', 'Explicit method')
ylabel('M')
xlabel('H_{applied} [T]')
%% Eigenvalue analysis not yet implemented
% open('EigTestDirectEq01.fig')
% hold on
% plot(linspace(-MaxH,MaxH,numel(Mx2)),-AllVV(NN+1,:).','+r')
% grid on
% PlotTheSolutionWithLoop(ThisPreSetFileName)
% load()