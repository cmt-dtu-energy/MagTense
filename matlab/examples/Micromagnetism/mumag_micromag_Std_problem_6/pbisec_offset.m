% The function finds perpendicular bisector between two points in 2D/3D
% Hyongju Park / hyongju@gmail.com
% input: two points in 2D/3D
% output: inequality Ax <= b

function [A,b] = pbisec_offset(x1, x2,d)

% middle_pnt = mean([x1;x2],1);
% middle_pnt = .3.*x1 + .7.*x2 ;
% middle_pnt = .7.*x1 + .3.*x2 ;
n_vec = (x2 - x1) / norm(x2 - x1);
middle_pnt = mean([x1;x2],1) - d.*n_vec ;

Ad = n_vec;
bd = dot(n_vec,middle_pnt);

if Ad * x1' <= bd
    A = Ad;
    b = bd;
else
    A = -Ad;
    b = -bd;
end
