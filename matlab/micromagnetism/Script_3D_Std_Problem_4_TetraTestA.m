function Script_3D_Std_Problem_4_TetraTestA
close all ;
disp('Running Matlab model')

%% Calculates dynamic solution to mumag std prob 4
clearvars -except fig1 MySim
mu0 = 4*pi*1e-7;

%% addpaths
addpath('../MEX_files');
addpath('../util');

%% TetraMesh
thisGridL = [125e-9,125e-9/2,125e-9/4] ;
GeomInfo = DefaultMicroMagProblem(2,2,2) ;

model = CreateTetraMesh(thisGridL,20e-9) ;

figure ; pdeplot3D(model,'FaceAlpha',0.1) ;
TetraMeshFileName = 'TestTetraMesh01.mat' ;
save(TetraMeshFileName,'model') ;

%% Create Preset

NIST_field = 1;

resolution = [size(model.Mesh.Elements,2),1,1];
MySim_ini = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
%         MySim_ini.A0 = 0  ;
%         MySim_ini.Ms = 1e-3 ;
MySim_ini.grid_L = thisGridL;% NOT STANDARD
MySim_ini.RecomputeInteractionMatrices = 0 ;
MySim_ini.ExternalMesh = 1 ;
MySim_ini.MeshType = 'Tetra' ;
MySim_ini.ExternalMeshFileName = TetraMeshFileName ;
MySim_ini.DemagTensorFileName = 'TetraInteractionMatrices.mat' ;
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

MySim_ini = MySim_ini.setHext( HextFct, linspace(0,100e-9,20) ); % was 2000 !
%     for n=1:3 ; MySim_ini.HextFct{n} = @(t) (1e-9-t)' .* HystDir(n) .* (t<1e-9)' ; end



%% Test Demag Field
load('CompareInteractionMatricesB.mat')
InteractionMatricesCartesian = InteractionMatrices ;
load(MySim_ini.DemagTensorFileName) ;
'' ;
DemagTensor =  InteractionMatrices.DemagTensor ; Mfact = 1 ;
DemagTensor0 =  InteractionMatricesCartesian.DemagTensor ; Mfact = 1 ;
AA.HmX = @(Sx,Sy,Sz,t) - Mfact.*(DemagTensor.KglobXX{1}*Sx+DemagTensor.KglobXY{1}*Sy+DemagTensor.KglobXZ{1}*Sz) ;
AA.HmY = @(Sx,Sy,Sz,t) - Mfact.*(DemagTensor.KglobXY{1}*Sx+DemagTensor.KglobYY{1}*Sy+DemagTensor.KglobYZ{1}*Sz) ;
AA.HmZ = @(Sx,Sy,Sz,t) - Mfact.*(DemagTensor.KglobXZ{1}*Sx+DemagTensor.KglobYZ{1}*Sy+DemagTensor.KglobZZ{1}*Sz) ;

AA2.HmX = @(Sx,Sy,Sz,t) - Mfact.*(DemagTensor0.KglobXX{1}*Sx+DemagTensor0.KglobXY{1}*Sy+DemagTensor0.KglobXZ{1}*Sz) ;
AA2.HmY = @(Sx,Sy,Sz,t) - Mfact.*(DemagTensor0.KglobXY{1}*Sx+DemagTensor0.KglobYY{1}*Sy+DemagTensor0.KglobYZ{1}*Sz) ;
AA2.HmZ = @(Sx,Sy,Sz,t) - Mfact.*(DemagTensor0.KglobXZ{1}*Sx+DemagTensor0.KglobYZ{1}*Sy+DemagTensor0.KglobZZ{1}*Sz) ;

IIin = (abs(InteractionMatrices.X)<(GeomInfo.grid_L(1)/4)) & (abs(InteractionMatrices.Y)<(GeomInfo.grid_L(2)/4)) ;

SigmaX = 0.*InteractionMatrices.X ; SigmaY = 0.*InteractionMatrices.X ; SigmaZ = 0.*InteractionMatrices.X ;
SigmaX(IIin) = 1 ;
ThisHmX = AA.HmX(SigmaX,SigmaY,SigmaZ,0) ; % fine
ThisHmY = AA.HmY(SigmaX,SigmaY,SigmaZ,0) ; % fine
ThisHmZ = AA.HmZ(SigmaX,SigmaY,SigmaZ,0) ; % fine

figure ; hold on
quiver3(InteractionMatrices.X,InteractionMatrices.Y,InteractionMatrices.Z,SigmaX,SigmaY,SigmaZ,'color','k') ;
quiver3(InteractionMatrices.X,InteractionMatrices.Y,InteractionMatrices.Z,ThisHmX,ThisHmY,ThisHmZ,'color','b') ;
axis equal
'' ;
A2 = InteractionMatrices.A2;
Jfact = 1 ;
AA.HjX = @(Sx,Sy,Sz,t) - (2*Jfact).*(A2*Sx) ;
AA.HjY = @(Sx,Sy,Sz,t) - (2*Jfact).*(A2*Sy) ;
AA.HjZ = @(Sx,Sy,Sz,t) - (2*Jfact).*(A2*Sz) ;

ThisHjX = AA.HjX(SigmaX,SigmaY,SigmaZ,0) ;
ThisHjY = AA.HjY(SigmaX,SigmaY,SigmaZ,0) ;
ThisHjZ = AA.HjZ(SigmaX,SigmaY,SigmaZ,0) ;

figure ; hold on
quiver3(InteractionMatrices.X,InteractionMatrices.Y,InteractionMatrices.Z,SigmaX,SigmaY,SigmaZ,'color','k') ;
quiver3(InteractionMatrices.X,InteractionMatrices.Y,InteractionMatrices.Z,ThisHjX,ThisHjY,ThisHjZ,'color','b') ;
axis equal
figure ; plot(ThisHjX(:),'.')
figure ; plot(ThisHmX(:),'.')
'' ;

return

%% Run the simulation
MySim_ini = MySim_ini.setSolverType( 'UseExplicitSolver' );
[SigmaInit,~,InteractionMatrices] = ComputeTheSolution(MySim_ini);
figure; quiver3(InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:),MySim_ini.m0(:,1),MySim_ini.m0(:,2),MySim_ini.m0(:,3)); axis equal;  title('Starting state - Matlab')
for k=1:size(SigmaInit,1)
    Sigma = SigmaInit(k,:).' ;
    NN = round(numel(Sigma)/3) ;
    
    SigmaX = Sigma(0*NN+[1:NN]) ;
    SigmaY = Sigma(1*NN+[1:NN]) ;
    SigmaZ = Sigma(2*NN+[1:NN]) ;
    if k==1 | k==size(SigmaInit,1)
        figure; quiver3(InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:),SigmaX(:),SigmaY(:),SigmaZ(:)); axis equal;  title('Starting state - Matlab')
    end
    SigmaN = sqrt(SigmaX.^2+SigmaY.^2+SigmaZ.^2) ;
    Mx(k) = mean(SigmaX./SigmaN) ;
    My(k) = mean(SigmaY./SigmaN) ;
    Mz(k) = mean(SigmaZ./SigmaN) ;
end

if (MySim_ini.ShowTheResult)
    figure ; plot(Mx,'o') ;
    %     figure; quiver3(InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:),SigmaX(:),SigmaY(:),SigmaZ(:)); axis equal;  title('Starting state - Matlab')
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
for n=1:3 ; MySim_dym.HextFct{n} = @(t) (t>-1)' .*HystDir(n) ; end ;
MySim_dym.m0(:) = SigmaInit(end,:);

MySim_dym = MySim_dym.setSolverType( 'UseDynamicSolver' );
MySim_dym.grid_L = MySim_ini.grid_L ;
tic
SigmaSol1 = ComputeTheSolution(MySim_dym);
toc
clear Mx My Mz
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
figure(figure1) ;
end