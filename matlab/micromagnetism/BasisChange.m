function A = BasisChange(A,transMat)
    A = transMat.'*A*transMat;
end