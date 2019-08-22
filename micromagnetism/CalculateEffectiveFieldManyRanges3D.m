function ThisHeff = CalculateEffectiveFieldManyRanges3D(t,Sigma,AvrgMatrix,CopyMatrix,Mfact,DemagTensor,AA)

NN = round(numel(Sigma)/3) ;
% N = round(NN^(1/3)) ;
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
    '' ; 

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

%%
ThisHeffX = ThisHjX + ThisHhX + ThisHmX + ThisHkX ;
ThisHeffY = ThisHjY + ThisHhY + ThisHmY + ThisHkY ;
ThisHeffZ = ThisHjZ + ThisHhZ + ThisHmZ + ThisHkZ ;

ThisHeff = [ThisHeffX;ThisHeffY;ThisHeffZ] ;


