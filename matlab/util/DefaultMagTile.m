classdef DefaultMagTile
    
    %% Defines a class with the default settings for running a magnetostatic
    %% problem using the Fortran implementation in MagTense.
    %% Note that the naming convention is dictated in the Fortran subroutine 
    %% getProblemFieldNames 

properties
    r0
    theta0
    z0

    dr
    dtheta
    dz

    %unit vector in the easy-axis direction    
    u_ea
    %unit vector in one of the directions orthogonal to the easy axis
    u_oa1
    %unit vector in the other direction orthogonal to the easy axis
    u_oa2
    %the remanence in the direction of the easy axis given no internal field
    Mrem
    %The permeability in the easy-axis direction.
    mu_r_ea
    %The permeability in the directions orthogonal to the easy-axis direction    
    mu_r_oa
    %The magnetization vector
    M  

    %center position of the tile (r0,theta0,z0) are defined with
    %respect to this, although it is given in Cartesian coords
    offset
    rotAngles
    abc
    
    %default index into the state function array if such is used
    stfcnIndex
    %if equal to zero the tile is not included in the iteration. Else it is
    inclIter
    
    %The rotation axis of a spheroid.
    %It can be either the axis of symmetry, 'symm', or the c-axis, 'c'.
    rot_axis
    
    %The rotation axis of a spheroid in vector coordinates.
    %The values given are the theta, phi and psi angles.
    %The rotAngles are computed from this property in the code.
    ax
    
    color
    %rotation of the magnet
    %tile.theta_off = 0;
    
    %tile.graphRotxAng = 0;
    
    %whether to exploit symmetry or not (0 not, 1 do) 
    useSymm
    %symmetry operation. 1 for symmetry and -1 for anti-symmetry (about the
    %xy, xz and yz planes, respectively)
    symmOps
    
    %latest relative error of the iteration (set to dummy value initially)
    Mrel
    
    %vertices used when the tile is a tetrahedron (and later, when it can
    %be a general surface made of triangular elements)
    vertices
    
    %boundary conditions represented as a list of properties (a struct
    %array)
    %refer to MagTenseTransientGeometry for the definitions of the various
    %fields in this struct array
    bdryCdts
end

properties (SetAccess=private,GetAccess=public)
    %1 = cylindrical piece, 2 = prism, 3 = ellipsoid
    tileType
     
    %1 = hard magnet, 2 = soft magnet
    magnetType    
end

methods
    function obj = DefaultMagTile()
        obj.r0 = 0;
        obj.theta0 = 0;
        obj.z0 = 0;
    
        obj.dr = 0;
        obj.dtheta = 0;
        obj.dz = 0;
    
        %unit vector in the easy-axis direction    
        obj.u_ea = [1,0,0];
        %unit vector in one of the directions orthogonal to the easy axis
        obj.u_oa1 = [0,1,0];
        %unit vector in the other direction orthogonal to the easy axis
        obj.u_oa2 = [0,0,1];
        %the remanence in the direction of the easy axis given no internal field
        obj.Mrem = 0;
        %The permeability in the easy-axis direction.
        obj.mu_r_ea = 1.0;
        %The permeability in the directions orthogonal to the easy-axis direction    
        obj.mu_r_oa = 1.0;
        %The magnetization vector
        obj.M = [0,0,0];
    
        %1 = cylindrical piece, 2 = prism, 3 = ellipsoid
        obj = obj.setMagTileType('Cylinder');
    
        %center position of the tile (r0,theta0,z0) are defined with
        %respect to this, although it is given in Cartesian coords
        obj.offset = [0,0,0];
        obj.rotAngles = [0,0,0];
        obj.abc = [0,0,0];
        %1 = hard magnet, 2 = soft magnet
        obj = obj.setMagnetType('Hard');
        %default index into the state function array if such is used
        obj.stfcnIndex = int32(1);
        %if equal to zero the tile is not included in the iteration. Else it is
        obj.inclIter = int32(1);
        
        %The rotation axis of a spheroid.
        %It can be either the axis of symmetry, 'symm', or the c-axis, 'c'.
        obj.rot_axis = 'symm';
        
        %The rotation axis of a spheroid in vector coordinates.
        %The values given are the theta, phi and psi angles.
        %The rotAngles are computed from this property in the code.
        obj.ax = [0 0 1];
        
        obj.color = [1,0,0];
        %rotation of the magnet
        %tile.theta_off = 0;
        
        %tile.graphRotxAng = 0;
        
        %whether to exploit symmetry or not (0 not, 1 do) 
        obj.useSymm = int32(0);
        %symmetry operation. 1 for symmetry and -1 for anti-symmetry (about the
        %xy, xz and yz planes, respectively)
        obj.symmOps = [1,1,1];
        
        %latest relative error of the iteration (set to dummy value initially)
        obj.Mrel = 0.;
        
        %vertices used when the tile is a tetrahedron (and later, when it can
        %be a general surface made of triangular elements)
        obj.vertices = zeros(3,4);
        
        %boundary conditions represented as a list of properties (a struct
        %array)
        %refer to MagTenseTransientGeometry for the definitions of the various
        %fields in this struct array
        obj.bdryCdts = struct( 'Type', -1, 'l', 0, 'A', 0, 'n_ind', 0, 'FaceID', -1, 'bdryFun', @dummyfun );
    end
    
    function obj = setMagTileType( obj, type_var )
        switch ( type_var )
            case {'Cylinder','cylinder'}
                obj.tileType = MagTenseTilesUtil.getMagTileType( 'Cylinder' );
            case {'Prism','prism'}
                obj.tileType = MagTenseTilesUtil.getMagTileType( 'Prism' );
            case {'Circpiece','circpiece'}
                obj.tileType = MagTenseTilesUtil.getMagTileType( 'Circpiece' );
            case {'Circpieceinv','circpieceinv'}
                obj.tileType = MagTenseTilesUtil.getMagTileType( 'Circpieceinv' );
            case {'Tetrahedron','tetrahedron'}
                obj.tileType = MagTenseTilesUtil.getMagTileType( 'Tetrahedron' );
            case {'Sphere','sphere'}
                obj.tileType = MagTenseTilesUtil.getMagTileType( 'Sphere' );
             case {'Spheroid','spheroid'}
                obj.tileType = MagTenseTilesUtil.getMagTileType( 'Spheroid' );
            case {'Planarcoil','planarcoil'}
                obj.tileType = MagTenseTilesUtil.getMagTileType( 'Planarcoil' );
        end
    end

    function res = getMagTileType( obj )
        switch obj.tileType 
            case MagTenseTilesUtil.getMagTileType( 'Cylinder' )
                res = 'Cylinder';
            case MagTenseTilesUtil.getMagTileType( 'Prism' )
                res = 'Prism';
            case MagTenseTilesUtil.getMagTileType( 'Circpiece' )
                res = 'Circpiece';
            case MagTenseTilesUtil.getMagTileType( 'Circpieceinv' )
                res = 'Circpieceinv';
            case MagTenseTilesUtil.getMagTileType( 'Tetrahedron' )
                res = 'Tetrahedron';
            case MagTenseTilesUtil.getMagTileType( 'Sphere' )
                res = 'Sphere';
            case MagTenseTilesUtil.getMagTileType( 'Spheroid' )
                res = 'Spheroid';
            case MagTenseTilesUtil.getMagTileType( 'Planarcoil' )
                res = 'Planarcoil';
        end
    end

    function obj = setMagnetType( obj, type_var )
        switch type_var 
            case {'Hard','hard'}
                obj.magnetType = MagTenseTilesUtil.getMagnetType( 'Hard' );
            case {'Soft','soft'}
                obj.magnetType = MagTenseTilesUtil.getMagnetType( 'Soft' );
            case {'Soft_const_mur','soft_const_mur'}
                obj.magnetType = MagTenseTilesUtil.getMagnetType( 'Soft_const_mur' );
        end
    end

    %Override struct function for a final check before handing to Fortran
    function obj2 = struct(obj)
        warning('off','MATLAB:structOnObject')
        for i = 1:length(obj)
            obj2(i)=builtin('struct',obj(i)); % Actual struct conversion
        end
        warning('on','MATLAB:structOnObject')
    end
    
end
end