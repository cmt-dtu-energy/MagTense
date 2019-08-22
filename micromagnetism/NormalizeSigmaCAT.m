function [Sigma] = NormalizeSigmaCAT(Sigma)
NN = round(numel(Sigma)/3) ;
% N = round(NN^(1/3)) ;
SigmaX = Sigma(0*NN+[1:NN]) ;
SigmaY = Sigma(1*NN+[1:NN]) ;
SigmaZ = Sigma(2*NN+[1:NN]) ;

SigmaN = sqrt(SigmaX.^2+SigmaY.^2+SigmaZ.^2) ;
SigmaX = SigmaX./SigmaN ;
SigmaY = SigmaY./SigmaN ;
SigmaZ = SigmaZ./SigmaN ;

Sigma = [SigmaX;SigmaY;SigmaZ];