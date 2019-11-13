function [Kxx,Kxy,Kxz,Kyx,Kyy,Kyz,Kzx,Kzy,Kzz] = CreateAnisotropyMatrixBase(Kx,Ky,Kz)
% the function creates the matrices for the anisotropy term
% when only one "grain" is present
% Kx,Ky, and Kz are the three components of the easy axis vector

Kxx = Kx.*Kx ;
Kxy = Kx.*Ky ;
Kxz = Kx.*Kz ;

Kyx = Ky.*Kx ;
Kyy = Ky.*Ky ;
Kyz = Ky.*Kz ;

Kzx = Kz.*Kx ;
Kzy = Kz.*Ky ;
Kzz = Kz.*Kz ;