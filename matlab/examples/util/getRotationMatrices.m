%%ang is an array with three elements (rotation angles about x-, y- and
%%z-axes, respectively). Angles are assumed in radians
%%Rerturns a cell array with three elements: the rotation matrices for
%%rotation about the x-, y- and z-axes, respectively.
function [RotOut] = getRotationMatrices( ang )

RotOut = cell(3,1);

%rotation about x-axis
RotOut{1} = [ 1, 0, 0; 0, cos(ang(1)), -sin(ang(1)); 0, sin(ang(1)), cos(ang(1))];

%rotation about y-axis
RotOut{2} = [cos(ang(2)), 0, sin(ang(2)); 0, 1, 0; -sin(ang(2)), 0, cos(ang(2))];

%rotation about z-axis
RotOut{3} = [ cos(ang(3)), -sin(ang(3)), 0; sin(ang(3)), cos(ang(3)), 0; 0, 0, 1];

end