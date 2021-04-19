function TheCross = CrossProdCAT(Sigma,ThisHeff)
NN = round(numel(Sigma)/3) ;
% N = round(NN^(1/3)) ;
SigmaX = Sigma(0*NN+[1:NN]) ;
SigmaY = Sigma(1*NN+[1:NN]) ;
SigmaZ = Sigma(2*NN+[1:NN]) ;

ThisHeffX = ThisHeff(0*NN+[1:NN]) ;
ThisHeffY = ThisHeff(1*NN+[1:NN]) ;
ThisHeffZ = ThisHeff(2*NN+[1:NN]) ;


% Compute m x heff (Precession term)
TheCrossX = -(SigmaY.*ThisHeffZ - SigmaZ.*ThisHeffY) ;
TheCrossY = -(SigmaZ.*ThisHeffX - SigmaX.*ThisHeffZ) ;
TheCrossZ = -(SigmaX.*ThisHeffY - SigmaY.*ThisHeffX) ;

TheCross = [TheCrossX;TheCrossY;TheCrossZ] ;