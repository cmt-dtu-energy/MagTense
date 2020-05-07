function [Mx,My,Mz] = ComputeMagneticMomentGeneralMesh(SigmaAll,Vols)
% Vols = prod(dims,2) ;
sumVol = sum(Vols) ;
for k=1:size(SigmaAll,1) 
    Sigma = SigmaAll(k,:).' ;
    NN = round(numel(Sigma)/3) ;

    SigmaX = Sigma(0*NN+[1:NN]) ;
    SigmaY = Sigma(1*NN+[1:NN]) ;
    SigmaZ = Sigma(2*NN+[1:NN]) ;
%     SigmaN = sqrt(SigmaX.^2+SigmaY.^2+SigmaZ.^2) ;
    Mx(k) = sum(Vols.*SigmaX)./sumVol ;
    My(k) = sum(Vols.*SigmaY)./sumVol ;
    Mz(k) = sum(Vols.*SigmaZ)./sumVol ;
end