%Creates and returnes a default tile that is compatible with the
%Fortran/mex interface
function tile = getDefaultMagTile()

    tile = struct();
    
    tile.r0 = 0;
    tile.theta0 = 0;
    tile.z0 = 0;

    tile.dr = 0;
    tile.dtheta = 0;
    tile.dz = 0;

    %unit vector in the easy-axis direction    
    tile.u_ea = [1,0,0];
    %unit vector in one of the directions orthogonal to the easy axis
    tile.u_oa1 = [0,1,0];
    %unit vector in the other direction orthogonal to the easy axis
    tile.u_oa2 = [0,0,1];
    %the remanence in the direction of the easy axis given no internal
    %field
    tile.Mrem = 0;
    %The permeability in the easy-axis direction. Value from Oxford
    %Instruments (private communcation via email with Hetal Patel 20-11-2017)
    tile.mu_r_ea = 1.0;
    %The permeability in the directions orthogonal to the easy-axis direction    
    tile.mu_r_oa = 1.0;
    %add the field containing the magnetization vector
    tile.M = [0,0,0];

    %1 = cylindrical piece, 2 = prism, 3 = ellipsoid
    tile.tileType = int32(1);

    %center position of the tile (r0,theta0,z0) are defined with
    %respect to this, although it is given in Cartesian coords
    tile.offset = [0,0,0];
    tile.rotAngles = [0,0,0];
    tile.abc = [0,0,0];
    %1 = hard magnet, 2 = soft magnet
    tile.magnetType = int32(1);
    %default index into the state function array if such is used
    tile.stfcnIndex = int32(1);
    %if equal to zero the tile is not included in the iteration. Else it is
    tile.inclIter = int32(1);
    
    tile.color = [1,0,0];
    %rotation of the magnet
    %tile.theta_off = 0;
    
    %tile.graphRotxAng = 0;
    
    %whether to exploit symmetry or not (0 not, 1 do) 
    tile.useSymm = int32(0);
    %symmetry operation. 1 for symmetry and -1 for anti-symmetry (about the
    %xy, xz and yz planes, respectively)
    tile.symmOps = [1,1,1];
    
    %latest relative error of the iteration (set to dummy value initially)
    tile.Mrel = 0.;
    
    %vertices used when the tile is a tetrahedron (and later, when it can
    %be a general surface made of triangular elements)
    tile.vertices = zeros(3,4);
end