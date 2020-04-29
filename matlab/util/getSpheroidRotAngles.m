function [RotAngles, ax] = getSpheroidRotAngles( tile )

if ~(strcmp(tile.rot_axis,'c') || strcmp(tile.rot_axis,'symm'))
    error('rot_axis has to be "c" or "symm"')
end

a = tile.abc(1);
b = tile.abc(2);
c = tile.abc(3);

%--- If the user specifies the c-axis as the axis of rotation
if strcmp(tile.rot_axis, 'c')

    if (a == b)
        ComputedAxisOfRevIndex = 3;
    elseif (b == c)
        ComputedAxisOfRevIndex = 1;
    elseif (c == a)
        ComputedAxisOfRevIndex = 2;
    end

    TheEye = eye(3) ;
    AxOriginal = TheEye(:,ComputedAxisOfRevIndex) ;
    TheMatrix = RotMatrixZ(tile.ax(2))*RotMatrixY(tile.ax(1))*RotMatrixZ(tile.ax(3)) ;
    ax = (TheMatrix*AxOriginal).' ;
else
    ax = tile.ax;
end

%% MagTenseRotationPart 

Rot_mat = [[cos(pi/2) 0 sin(pi/2)];[0 1 0];[-sin(pi/2) 0 cos(pi/2)]];
for i = 1:length(ax(:,1))
    ax(i,:) = Rot_mat*ax(i,:)';
end

%--- Calculate the spherical coordinates (which Kaspar calls yaw and pitch)
%--- The azimuthal angle is offset by pi/2.
rot = zeros(1,3);
rot(:,1) = -atan2(ax(:,2),ax(:,1));
rot(:,2) = acos(ax(:,3)./sqrt(ax(:,1).^2+ax(:,2).^2+ax(:,3).^2))-pi/2;
rot(:,3) =  0;

RotAngles = rot;


end

function TheMatrix = RotMatrixX(t)
c = @(t) cos(t) ;
s = @(t) sin(t) ;
TheMatrix = [1,0,0;0,c(t),-s(t);0,s(t),c(t)] ;
end

function TheMatrix = RotMatrixY(t)
c = @(t) cos(t) ;
s = @(t) sin(t) ;
TheMatrix = [c(t),0,s(t);0,1,0;-s(t),0,c(t)] ;
end

function TheMatrix = RotMatrixZ(t)
    c = @(t) cos(t) ;
    s = @(t) sin(t) ;
    TheMatrix = [c(t),-s(t),0;s(t),c(t),0;0,0,1] ;
end