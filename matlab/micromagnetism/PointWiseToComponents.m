function [Axx,Axy,Axz,Ayx,Ayy,Ayz,Azx,Azy,Azz] = PointWiseToComponents(A)
% A takes a vector where the three vector components of the same space point
% are consecutive: 
% m = [mx(1);my(1);mz(1);mx(2);my(2);mz(2);...;mx(N);my(N);mz(N)]
% this function separates the blocks for the different components.
% The inputs are:
% mx = [mx(1);mx(2);...;mx(N)] ;
% my = [my(1);my(2);...;my(N)] ;
% mz = [mz(1);mz(2);...;mz(N)] ;
% The outputs are organized in the same way

Axx = A(1:3:end,1:3:end) ;
Axy = A(1:3:end,2:3:end) ;
Axz = A(1:3:end,3:3:end) ;

Ayx = A(2:3:end,1:3:end) ;
Ayy = A(2:3:end,2:3:end) ;
Ayz = A(2:3:end,3:3:end) ;

Azx = A(3:3:end,1:3:end) ;
Azy = A(3:3:end,2:3:end) ;
Azz = A(3:3:end,3:3:end) ;