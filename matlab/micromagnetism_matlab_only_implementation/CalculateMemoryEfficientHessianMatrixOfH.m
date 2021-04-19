function HessH = CalculateMemoryEfficientHessianMatrixOfH(Heff,Sigma,HessGXX,HessGXY,HessGXZ,HessGYX,HessGYY,HessGYZ,HessGZX,HessGZY,HessGZZ) 
%% Disclaimer: May not actually be memory efficient. Result differs from CalculateHessianMatrixOfH by errors roughly 14 orders of magnitude smaller than input matrices.

NN = round(numel(Sigma)/3) ;
% N = round(NN^(1/3)) ;
SigmaX = Sigma(0*NN+[1:NN]) ;
SigmaY = Sigma(1*NN+[1:NN]) ;
SigmaZ = Sigma(2*NN+[1:NN]) ;

%%
[dfxdmx,dfydmx,dfzdmx,dfxdmy,dfydmy,dfzdmy,dfxdmz,dfydmz,dfzdmz,...
    dfxdmxdmx,dfydmxdmx,dfzdmxdmx,dfxdmxdmy,dfydmxdmy,dfzdmxdmy,dfxdmxdmz,dfydmxdmz,dfzdmxdmz,...
    dfxdmydmx,dfydmydmx,dfzdmydmx,dfxdmydmy,dfydmydmy,dfzdmydmy,dfxdmydmz,dfydmydmz,dfzdmydmz,...
    dfxdmzdmx,dfydmzdmx,dfzdmzdmx,dfxdmzdmy,dfydmzdmy,dfzdmzdmy,dfxdmzdmz,dfydmzdmz,dfzdmzdmz] = FirstAndSecondDerivativeOfFcomponents(SigmaX,SigmaY,SigmaZ) ;

%% (Hopefully) Exploits In-place operation syntax described at: https://blogs.mathworks.com/loren/2007/03/22/in-place-operations-on-data/
HessH = zeros(3*NN,3*NN);
HessH = AddHessDiag(HessH,HessGXX,dfxdmx,dfxdmy,dfxdmz,NN) ;
HessH = AddHessDiag(HessH,HessGYY,dfydmx,dfydmy,dfydmz,NN) ;
HessH = AddHessDiag(HessH,HessGZZ,dfzdmx,dfzdmy,dfzdmz,NN) ;
HessH = AddHessNonDiag(HessH,HessGXY,dfxdmx,dfxdmy,dfxdmz,dfydmx,dfydmy,dfydmz,NN) ;
HessH = AddHessNonDiag(HessH,HessGXZ,dfxdmx,dfxdmy,dfxdmz,dfzdmx,dfzdmy,dfzdmz,NN) ;
HessH = AddHessNonDiag(HessH,HessGYX,dfydmx,dfydmy,dfydmz,dfxdmx,dfxdmy,dfxdmz,NN) ;
HessH = AddHessNonDiag(HessH,HessGYZ,dfydmx,dfydmy,dfydmz,dfzdmx,dfzdmy,dfzdmz,NN) ;
HessH = AddHessNonDiag(HessH,HessGZX,dfzdmx,dfzdmy,dfzdmz,dfxdmx,dfxdmy,dfxdmz,NN) ;
HessH = AddHessNonDiag(HessH,HessGZY,dfzdmx,dfzdmy,dfzdmz,dfydmx,dfydmy,dfydmz,NN) ;
       
HessH = AddHessAIP(HessH,dfxdmxdmx,dfydmxdmx,dfzdmxdmx,dfxdmxdmy,dfydmxdmy,dfzdmxdmy,dfxdmxdmz,dfydmxdmz,dfzdmxdmz,...
    dfxdmydmx,dfydmydmx,dfzdmydmx,dfxdmydmy,dfydmydmy,dfzdmydmy,dfxdmydmz,dfydmydmz,dfzdmydmz,...
    dfxdmzdmx,dfydmzdmx,dfzdmzdmx,dfxdmzdmy,dfydmzdmy,dfzdmzdmy,dfxdmzdmz,dfydmzdmz,dfzdmzdmz,Heff,NN) ;