function model = getSpheroidComsolCompare( tile )

%--- Start the Comsol model
model = getGeneralComsolCompareBegin();

%% The rotation axis is the c-axis of the spheroid  
if strcmp(tile.rot_axis, 'c')
    % The rotation about the C axis can be performed in any moment,
    % since the C axis rotates with the object 
    % (i.e. it's not given in the global coord. syst.).
    % However, since we know that at the beginning the c axis 
    % is oriented in the z direction, it is more convenient to do it first 
    % so we can use the z-Rotation matrix.        
    model.param.set('theTheta', num2str(tile.ax(1)) , 'Inclination ( rotation about y axis - performed before azimuth )');
    model.param.set('thePhi'  , num2str(tile.ax(2)) , 'Azimuth (  rotation about z axis - performed after inclination )');
    model.param.set('thePsi'  , num2str(tile.ax(3)) , 'Rotation ( rotation about C axis - order does not matter )');

    a = tile.abc(1);
    b = tile.abc(1);
    c = tile.abc(1);
end

%% The rotation axis is the symmetry-axis of the spheroid
if strcmp(tile.rot_axis, 'symm')
    model.param.set('theTheta', num2str(acos(tile.ax(3)/sqrt(tile.ax(1)^2+tile.ax(2)^2+tile.ax(3)^2))) , 'Inclination ( rotation about y axis - performed before azimuth )');
    model.param.set('thePhi'  , num2str(atan2(tile.ax(2),tile.ax(1))) , 'Azimuth (  rotation about z axis - performed after inclination )');
    model.param.set('thePsi'  , num2str(0) , 'Rotation ( rotation about C axis - order does not matter )');
    
    %--- Rearrange dimensions a,b,c so that the different axis is always the c-axis
    %--- i.e. a & b axes: equal to each other, while c axis: The one that is different
    if tile.abc(1)==tile.abc (2)
        %--- Do nothing
    elseif tile.abc(1) == tile.abc (3)
        a = tile.abc(1) ; 
        b = tile.abc(3) ; 
        c = tile.abc(2) ;
    elseif tile.abc(2) == tile.abc (3)
        a = tile.abc(2) ; 
        b = tile.abc(3) ; 
        c = tile.abc(1) ;
    end
end

%% Set the common ellipsoid parameters
model.param.set('A_axis', num2str(a), 'length along a axis: i.e. the one that is originally in the x direction');
model.param.set('B_axis', num2str(b), 'length along b axis: i.e. the one that is originally in the y direction');
model.param.set('C_axis', num2str(c), 'length along c axis: i.e. the one that is originally in the z direction');

model.param.set('x0',  num2str(tile.offset(1))) ;
model.param.set('y0',  num2str(tile.offset(2))) ;
model.param.set('z0',  num2str(tile.offset(3))) ;

model.component('comp1').geom('geom1').create('elp1', 'Ellipsoid');
model.component('comp1').geom('geom1').feature('elp1').set('pos', {'x0' 'y0' 'z0'});
model.component('comp1').geom('geom1').feature('elp1').set('rot', 'thePsi*(180/pi)');
model.component('comp1').geom('geom1').feature('elp1').set('axis', {'theTheta*(180/pi)' 'thePhi*(180/pi)'});
model.component('comp1').geom('geom1').feature('elp1').set('semiaxes', {'A_axis' 'B_axis' 'C_axis'});

%% Call the general Comsol model
model = getGeneralComsolCompare( model, tile );

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