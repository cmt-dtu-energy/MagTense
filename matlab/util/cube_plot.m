%function cube_plot(origin,dims,rot,order,color,grphRotAng)
function cube_plot( tile, order )
origin = tile.offset;
dims = tile.abc;
rot = tile.rotAngles;
color = tile.color;
if isfield(tile,'graphRotxAng')
    grphRotAng = tile.graphRotxAng;
else
    grphRotAng = 0;
end
%addpath('../utils/');
% CUBE_PLOT plots a cube with dimension of X, Y, Z.
%
% INPUTS:
% origin = set origin point for the cube in the form of [x,y,z].
% X      = cube length along x direction.
% Y      = cube length along y direction.
% Z      = cube length along z direction.
% color  = STRING, the color patched for the cube.
%         List of colors
%         b blue
%         g green
%         r red
%         c cyan
%         m magenta
%         y yellow
%         k black
%         w white
% OUPUTS:
% Plot a figure in the form of cubics.
%
% EXAMPLES
% cube_plot(2,3,4,'red')
%
% ------------------------------Code Starts Here------------------------------ %
% Define the vertexes of the unit cubic
X = dims(1);
Y = dims(2);
Z = dims(3);

ver = [1 1 0;
    0 1 0;
    0 1 1;
    1 1 1;
    0 0 1;
    1 0 1;
    1 0 0;
    0 0 0];
%move the vertices so the cube is centered on origo
for i=1:3
    ver(:,i) = ver(:,i) - 0.5;
end
%  Define the faces of the unit cubic
fac = [1 2 3 4;
    4 3 5 6;
    6 7 8 5;
    1 2 8 7;
    6 7 1 4;
    2 3 5 8];
%rotate the cube according to the angles specified in rot and the order of
%rotation specified in order
ver = rotateVertices( rot, ver, order );

%move the vertices back to have lower left corner in (0,0,0)
for i=1:3
    ver(:,i) = ver(:,i) + 0.5;
end

cube = [ver(:,1)*X-X/2+origin(1),ver(:,2)*Y-Y/2+origin(2),ver(:,3)*Z-Z/2+origin(3)];
t = hgtransform;

patch('Faces',fac,'Vertices',cube,'FaceColor',color,'parent',t);

Rx=makehgtform('xrotate',grphRotAng);
t.Matrix=Rx;
end
% ------------------------------Code Ends Here-------------------------------- %
