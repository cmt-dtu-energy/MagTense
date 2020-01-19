% function dSigma = CalculateDeltaSigmaManyRanges3D(t,Sigma,N,AvrgMatrix,CopyMatrix,Mfact,DemagTensor,alpha,gamma,PlotStruct,AA,dx,dy,dz)
function dSigma = CalculateDeltaSigmaManyRanges3D(t,Sigma,AvrgMatrix,CopyMatrix,Mfact,DemagTensor,alpha,gamma,AA)

NN = round(numel(Sigma)/3) ;

SigmaX = Sigma(0*NN+[1:NN]) ;
SigmaY = Sigma(1*NN+[1:NN]) ;
SigmaZ = Sigma(2*NN+[1:NN]) ;

%% Demagnetization: Long range (coarse grids)

HmXcTOT = zeros(NN,1) ;
HmYcTOT = zeros(NN,1) ;
HmZcTOT = zeros(NN,1) ;
for k=2:numel(AvrgMatrix)
    
    % get sigma over coarse grid (average)
    SigmaXC = AvrgMatrix{k}*SigmaX ;
    SigmaYC = AvrgMatrix{k}*SigmaY ;
    SigmaZC = AvrgMatrix{k}*SigmaZ ;
    
    % calculate demag. field over coarse grid   
    HmXc = DemagTensor.KglobXX{k}*SigmaXC+DemagTensor.KglobXY{k}*SigmaYC+DemagTensor.KglobXZ{k}*SigmaZC ;
    HmYc = DemagTensor.KglobYX{k}*SigmaXC+DemagTensor.KglobYY{k}*SigmaYC+DemagTensor.KglobYZ{k}*SigmaZC ;
    HmZc = DemagTensor.KglobZX{k}*SigmaXC+DemagTensor.KglobZY{k}*SigmaYC+DemagTensor.KglobZZ{k}*SigmaZC ;

    % get demag. field over fine grid (copy)   
    HmXcTOT = HmXcTOT - Mfact*CopyMatrix{k}*HmXc ;  % Coarser
    HmYcTOT = HmYcTOT - Mfact*CopyMatrix{k}*HmYc ;  % Coarser
    HmZcTOT = HmZcTOT - Mfact*CopyMatrix{k}*HmZc ;  % Coarser
end

%%  Exchange 
ThisHjX = AA.HjX(SigmaX,SigmaY,SigmaZ,t) ;
ThisHjY = AA.HjY(SigmaX,SigmaY,SigmaZ,t) ;
ThisHjZ = AA.HjZ(SigmaX,SigmaY,SigmaZ,t) ;

%% External field
ThisHhX = AA.HhX(SigmaX,SigmaY,SigmaZ,t) ;
ThisHhY = AA.HhY(SigmaX,SigmaY,SigmaZ,t) ;
ThisHhZ = AA.HhZ(SigmaX,SigmaY,SigmaZ,t) ;

%% Demagnetization
if numel(CopyMatrix)>0
    ThisHmX = AA.HmX(SigmaX,SigmaY,SigmaZ,t) + HmXcTOT ; % fine + coarser
    ThisHmY = AA.HmY(SigmaX,SigmaY,SigmaZ,t) + HmYcTOT ; % fine + coarser
    ThisHmZ = AA.HmZ(SigmaX,SigmaY,SigmaZ,t) + HmZcTOT ; % fine + coarser
else
    ThisHmX = AA.HmX(SigmaX,SigmaY,SigmaZ,t) ; % fine 
    ThisHmY = AA.HmY(SigmaX,SigmaY,SigmaZ,t) ; % fine
    ThisHmZ = AA.HmZ(SigmaX,SigmaY,SigmaZ,t) ; % fine
end
%% Anisotropy
ThisHkX = AA.HkX(SigmaX,SigmaY,SigmaZ,t) ;
ThisHkY = AA.HkY(SigmaX,SigmaY,SigmaZ,t) ;
ThisHkZ = AA.HkZ(SigmaX,SigmaY,SigmaZ,t) ;

ThisHeffX = ThisHjX + ThisHhX + ThisHmX + ThisHkX ;
ThisHeffY = ThisHjY + ThisHhY + ThisHmY + ThisHkY ;
ThisHeffZ = ThisHjZ + ThisHhZ + ThisHmZ + ThisHkZ ;

%% Calculate Precession and Damping terms from Heff

% Compute m x heff (Precession term)
TheCrossX = -(SigmaY.*ThisHeffZ - SigmaZ.*ThisHeffY) ;
TheCrossY = -(SigmaZ.*ThisHeffX - SigmaX.*ThisHeffZ) ;
TheCrossZ = -(SigmaX.*ThisHeffY - SigmaY.*ThisHeffX) ;

% Compute m x m x heff (Damping term)
ThisHeffX2 = +SigmaY.*TheCrossZ-SigmaZ.*TheCrossY ;
ThisHeffY2 = +SigmaZ.*TheCrossX-SigmaX.*TheCrossZ ;
ThisHeffZ2 = +SigmaX.*TheCrossY-SigmaY.*TheCrossX ;

%% Calculate time-derivative of Sigma

if isequal(class(alpha),'double') % alpha is time-independent
    dSigmaX = -alpha.*ThisHeffX2 - gamma.*TheCrossX ;
    dSigmaY = -alpha.*ThisHeffY2 - gamma.*TheCrossY ;
    dSigmaZ = -alpha.*ThisHeffZ2 - gamma.*TheCrossZ ;
    dSigma = [dSigmaX;dSigmaY;dSigmaZ] ;
    dSigmaRMS = sqrt(sum((dSigma./-alpha).^2)/NN) ;
else % alpha is time-dependent
    dSigmaX = -alpha(t).*ThisHeffX2 - gamma.*TheCrossX ;
    dSigmaY = -alpha(t).*ThisHeffY2 - gamma.*TheCrossY ;
    dSigmaZ = -alpha(t).*ThisHeffZ2 - gamma.*TheCrossZ ;
    dSigma = [dSigmaX;dSigmaY;dSigmaZ] ;
    dSigmaRMS = sqrt(sum((dSigma./-alpha(t)).^2)/NN) ;
end

%--- If we are using gpuArray, bring the magnetization back for ode45
dSigma = gather(dSigma);

%% Pass data using a globa variable
global TheData
LastT = TheData.LastT ;
LastPlottedT = TheData.LastPlottedT ;
TheData.LastT = t ;
if isfield(TheData,'LastDSigma')
    LastDSigma = TheData.LastDSigma ;
    UsePrevDSigma = 1 ;
else
    UsePrevDSigma = 0 ;    
end
TheData.LastDSigma = dSigma ;
TheData.dSigmaRMS = dSigmaRMS ;

% disp(['Here ',num2str(t)]) ;

