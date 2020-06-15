%returns the rotation matrix about the z-axis given ang in radians
function rot = rotz( ang )

rot = zeros(3,3);

rot(1,1) = cos(ang);
rot(1,2) = -sin(ang);
rot(2,1) = sin(ang);
rot(2,2) = cos(ang);
rot(3,3) = 1;
end