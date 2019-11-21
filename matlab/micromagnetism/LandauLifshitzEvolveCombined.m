function [SigmaSol,Applied_field,VV] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,Sigma)

%--- alpha and gamma are native functions in Matlab and must be initialized as variables
alpha = 0;
gamma = 0;

%% Evaluate all variables in the ProblemSetupStruct
names = fieldnames(ProblemSetupStruct);
for i=1:length(names)
    eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
end

%% Unfold the struct (NB: *after* unfolding ProblemSetupStruct)
DemagTensor = InteractionMatrices.DemagTensor;
N = InteractionMatrices.N;
A2 = InteractionMatrices.A2;
Kxx = InteractionMatrices.Kxx;
Kxy = InteractionMatrices.Kxy;
Kxz = InteractionMatrices.Kxz;
Kyx = InteractionMatrices.Kyx;
Kyy = InteractionMatrices.Kyy;
Kyz = InteractionMatrices.Kyz;
Kzx = InteractionMatrices.Kzx;
Kzy = InteractionMatrices.Kzy;
Kzz = InteractionMatrices.Kzz;
AvrgMatrix = InteractionMatrices.AvrgMatrix;
CopyMatrix = InteractionMatrices.CopyMatrix;

%% Effective field parameter prefactors

Jfact = A0/(MU0*Ms) ;   % "J" : exchange term
% Hfact = 1/MU0 ;         % "H" : external field term (b.c. user input is in Tesla)
Mfact = Ms ;            % "M" : demagnetization term
Kfact = K0/(MU0*Ms) ;   % "K" : anisotropy term

%% All Interaction Terms (function handles)

O = ones(prod(N),1) ;

% "J" : Effective field : exchange

AA.HjX = @(Sx,Sy,Sz,t) - (2*Jfact).*(A2*Sx) ;
AA.HjY = @(Sx,Sy,Sz,t) - (2*Jfact).*(A2*Sy) ;
AA.HjZ = @(Sx,Sy,Sz,t) - (2*Jfact).*(A2*Sz) ;

% "M" : Effective field : demagnetization (fine grid terms)

AA.HmX = @(Sx,Sy,Sz,t) - Mfact.*(DemagTensor.KglobXX{1}*Sx+DemagTensor.KglobXY{1}*Sy+DemagTensor.KglobXZ{1}*Sz) ;
AA.HmY = @(Sx,Sy,Sz,t) - Mfact.*(DemagTensor.KglobYX{1}*Sx+DemagTensor.KglobYY{1}*Sy+DemagTensor.KglobYZ{1}*Sz) ;
AA.HmZ = @(Sx,Sy,Sz,t) - Mfact.*(DemagTensor.KglobZX{1}*Sx+DemagTensor.KglobZY{1}*Sy+DemagTensor.KglobZZ{1}*Sz) ;

% "H" : Effective field : external field

if UseExplicitSolver | UseDynamicSolver
    AA.HhX = @(Sx,Sy,Sz,t) - HsX(t).*O ;
    AA.HhY = @(Sx,Sy,Sz,t) - HsY(t).*O ;
    AA.HhZ = @(Sx,Sy,Sz,t) - HsZ(t).*O ;
end

if UseImplicitSolver
    AA.HhX = @(Sx,Sy,Sz,t) - ddHsX*(t).*O ;
    AA.HhY = @(Sx,Sy,Sz,t) - ddHsY*(t).*O ;
    AA.HhZ = @(Sx,Sy,Sz,t) - ddHsZ*(t).*O ;

    %
    
    AA.ddHhX = @(Sx,Sy,Sz,t) - ddHsX.*O +0.*t ;
    AA.ddHhY = @(Sx,Sy,Sz,t) - ddHsY.*O +0.*t ;
    AA.ddHhZ = @(Sx,Sy,Sz,t) - ddHsZ.*O +0.*t ;
end

% "K" : Effective field : anisotropy

AA.HkX = @(Sx,Sy,Sz,t) - (2*Kfact).*((Kxx).*Sx + (Kxy).*Sy + (Kxz).*Sz) ;
AA.HkY = @(Sx,Sy,Sz,t) - (2*Kfact).*((Kyx).*Sx + (Kyy).*Sy + (Kyz).*Sz) ;
AA.HkZ = @(Sx,Sy,Sz,t) - (2*Kfact).*((Kzx).*Sx + (Kzy).*Sy + (Kzz).*Sz) ;

% Energies (these variables don't have any effects on the computation: they are just for plotting)

AA.Gj = @(Sx,Sy,Sz,Hx,Hy,Hz) (1/2)*(Sx.'*Hx + Sy.'*Hy + Sz.'*Hz) ; % -(J)*(Sx.'*A2*Sx + Sy.'*A2*Sy)  ;
AA.Gh = @(Sx,Sy,Sz,Hx,Hy,Hz)     + (Sx.'*Hx + Sy.'*Hy + Sz.'*Hz) ;
AA.Gm = @(Sx,Sy,Sz,Hx,Hy,Hz) (1/2)*(Sx.'*Hx + Sy.'*Hy + Sz.'*Hz) - 1 ; % MU.*(Sx.'*HHx(Sx,Sy) + Sy.'*HHy(Sx,Sy)) - 1  ;
AA.Gk = @(Sx,Sy,Sz,Hx,Hy,Hz) (1/2)*(Sx.'*Hx + Sy.'*Hy + Sz.'*Hz) ;

%% Hessian
if UseImplicitSolver | CalcEigenvalue
    NN = size(DemagTensor.KglobXX{1},1) ;
    HessGXX = - (2*Jfact).*(A2) - (Mfact).*DemagTensor.KglobXX{1} - (Kfact).*(Kxx)*eye(NN) ;
    HessGXY =                   - (Mfact).*DemagTensor.KglobXY{1} - (Kfact).*(Kxy)*eye(NN) ;
    HessGXZ =                   - (Mfact).*DemagTensor.KglobXZ{1} - (Kfact).*(Kxz)*eye(NN) ;
    
    HessGYX =                   - (Mfact).*DemagTensor.KglobYX{1} - (Kfact).*(Kyx)*eye(NN) ;
    HessGYY = - (2*Jfact).*(A2) - (Mfact).*DemagTensor.KglobYY{1} - (Kfact).*(Kyy)*eye(NN) ;
    HessGYZ =                   - (Mfact).*DemagTensor.KglobYZ{1} - (Kfact).*(Kyz)*eye(NN) ;
    
    HessGZX =                   - (Mfact).*DemagTensor.KglobZX{1} - (Kfact).*(Kzx)*eye(NN) ;
    HessGZY =                   - (Mfact).*DemagTensor.KglobZY{1} - (Kfact).*(Kzy)*eye(NN) ;
    HessGZZ = - (2*Jfact).*(A2) - (Mfact).*DemagTensor.KglobZZ{1} - (Kfact).*(Kzz)*eye(NN) ;
    
    for k=2:numel(AvrgMatrix)
        
        HessGXX = HessGXX - Mfact.*(CopyMatrix{k}*DemagTensor.KglobXX{k}*AvrgMatrix{k}) ;
        HessGXY = HessGXY - Mfact.*(CopyMatrix{k}*DemagTensor.KglobXY{k}*AvrgMatrix{k}) ;
        HessGXZ = HessGXZ - Mfact.*(CopyMatrix{k}*DemagTensor.KglobXZ{k}*AvrgMatrix{k}) ;
        
        HessGYX = HessGYX - Mfact.*(CopyMatrix{k}*DemagTensor.KglobYX{k}*AvrgMatrix{k}) ;
        HessGYY = HessGYY - Mfact.*(CopyMatrix{k}*DemagTensor.KglobYY{k}*AvrgMatrix{k}) ;
        HessGYZ = HessGYZ - Mfact.*(CopyMatrix{k}*DemagTensor.KglobYZ{k}*AvrgMatrix{k}) ;
        
        HessGZX = HessGZX - Mfact.*(CopyMatrix{k}*DemagTensor.KglobZX{k}*AvrgMatrix{k}) ;
        HessGZY = HessGZY - Mfact.*(CopyMatrix{k}*DemagTensor.KglobZY{k}*AvrgMatrix{k}) ;
        HessGZZ = HessGZZ - Mfact.*(CopyMatrix{k}*DemagTensor.KglobZZ{k}*AvrgMatrix{k}) ;
    end
    
    Heff = @(t,Sigma) CalculateEffectiveFieldManyRanges3D(t,Sigma,AvrgMatrix,CopyMatrix,Mfact,DemagTensor,AA) ;
    HHess = @(t,Sigma) CalculateMemoryEfficientHessianMatrixOfH(Heff(t,NormalizeSigmaCAT(Sigma)),Sigma,HessGXX,HessGXY,HessGXZ,HessGYX,HessGYY,HessGYZ,HessGZX,HessGZY,HessGZZ) ;
end

%% Global data
global TheData
InitialData.LastT = 0;
TheData = InitialData;

%% Define time-derivative function
if UseExplicitSolver | UseDynamicSolver
    TheData.dSigmaRMS = inf;    
    TheData.LastPlottedT = 0;
    dSigma2 = @(t,Sigma) CalculateDeltaSigmaManyRanges3D(t,Sigma,AvrgMatrix,CopyMatrix,Mfact,DemagTensor,alpha,gamma,AA) ;

    if ~exist('HeffLimMagn','var') 
        HeffLimMagn = -inf ;
    end
    ThatOutPutFunct = @(t,y,flag) MyOutputFunction(t,y,flag,UseImplicitSolver,HeffLimMagn) ;
end

if UseImplicitSolver
    InitialData.p = 0;
    TheData = InitialData;
    
    dSigma3 = @(t,Sigma) CalculateDeltaSigmaEquilibriumExplicitOut(t,Sigma,HHess,AA,MaxComputationalTimePerStep) ;
    
    ThatOutPutFunct = @(t,y,flag) MyOutputFunction(t,y,flag,UseImplicitSolver) ;
end


%% ODE
disp('Integrating Equation of Motion')
if UseDynamicSolver
    options = odeset('RelTol',1e-12) ;
    [t,SigmaSol] = ode45(dSigma2,t,Sigma,options);
%         [t,SigmaSol] = ode23(dSigma2,t,Sigma,options);
end
if UseExplicitSolver
    options = odeset('RelTol',1e-12,'OutputFcn',ThatOutPutFunct) ;
    [t,SigmaSol] = ode45(dSigma2,linspace(0,MaxT,nT),Sigma,options);
end

if UseImplicitSolver
    options = odeset('OutputFcn',ThatOutPutFunct,'RelTol',1e-3,'MaxStep',abs(t(1)-t(end))) ; % ,'MaxStep',0.01*abs(t(1)-t(end))) ; % ,'Events',ThatEventFunct) ;
    [t,SigmaSol] = ode45(dSigma3,t,Sigma,options);
end


Applied_field = [HsX(t), HsY(t), HsZ(t)+0.*t];

%% Calculate the Eigenvalue
if CalcEigenvalue
    ThisHHess5 = HHess(t(end),SigmaSol(end,:)') ;
    VV = CalculateHessianEigenvaluesFromHH(ThisHHess5) ;
    if UseImplicitStepsSolver
        if VV(NN+1) < 0
            error('TheEigenvalue has wrong sign') ;
        end
    end
else
   VV = [];
end

clear TheData
end
