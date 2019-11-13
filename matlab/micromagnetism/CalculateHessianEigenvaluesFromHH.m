function VVsort = CalculateHessianEigenvaluesFromHH(ThisHHess5)


% ThisHeffN2 = Heff(NormalizeSigmaCAT(Sigma)) ;
% ThisHHess5 = CalculateHessianMatrixOfH(ThisHeffN2,Sigma,HessGXX,HessGXY,HessGXZ,HessGYX,HessGYY,HessGYZ,HessGZX,HessGZY,HessGZZ) ;

if 0
    [~,VV] = eig(ThisHHess5) ;
    VV = diag(VV) ;
else
    VV = eig(ThisHHess5) ;
end
VV = real(VV) ; % ????????
VVsort = sort(VV) ;

