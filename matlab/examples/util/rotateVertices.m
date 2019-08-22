%%rot is an array with n elements, each representing an angle in radians
%%ver is an array with vertices assumed to be [m,3] in dimensions
%%order is an array with n elements each representing a rotation about an
%%axis, i.e. the elements of order should be 1, 2 or 3. All other values
%%are not allowed
function [ver] = rotateVertices( rot, ver, order )


RotOut = getRotationMatrices( rot );

for i=1:length(order)
    %rotation about the axis defined in order    
    ver = ( RotOut{order(i)} * ver')';
end


end