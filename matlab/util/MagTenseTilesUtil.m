
%defines a class for doing static evaluations of tiles
classdef MagTenseTilesUtil
   properties
       
   end
   
   methods (Static)
       function type = getMagTileType( str_type )
            switch ( str_type )
                case 'Cylinder'
                    type = int32( 1 );
                case 'Prism'
                    type = int32( 2 );
                case 'Circpiece'
                    type = int32( 3 );
                case 'Circpieceinv'
                    type = int32( 4 );
                case 'Tetrahedron'
                    type = int32( 5 );
                case 'Sphere'
                    type = int32( 6 );
                case 'Spheroid'
                    type = int32( 7 );
                case 'Planarcoil'
                    type = int32( 101 );
            end
        end
       
       function type = getMagnetType( str_type )
            switch ( str_type )
                case 'Hard'
                    type = int32( 1 );
                case 'Soft'
                    type = int32( 2 );
                case 'Soft_const_mur'
                    type = int32( 3 );
            end
        end
       

        function stateFunction = getDefaultMagStatStateFunction()
            stateFunction = struct();
            stateFunction.nT = int32(1);
            stateFunction.nH = int32(1);
            stateFunction.M = zeros(1,1);
            stateFunction.T = zeros(1,1);
            stateFunction.H = zeros(1,1);
        end

        function stateFcn = getStateFunctionSoftIron(name)
            %%Name is not used currently, but enables the possibility for having
            %%different state functions

            stateFcn = MagTenseTilesUtil.getDefaultMagStatStateFunction();
            %%setup the state function
            
            % fullfile keeps this working on both Windows and Linux; mfilename gives the
            % path of this file without its extension, so its folder is matlab/util
            util_dir = fileparts(mfilename('fullpath'));
            k = load(fullfile(util_dir, '..', '..', 'documentation', 'MH_Comsol_low_carb_annealed_extrap.txt'));
            
            stateFcn.nT = int32(3);
            stateFcn.nH = int32(length(k(:,1)));
            stateFcn.T = [200,300,400];
            stateFcn.H = k(:,1);
            stateFcn.M = zeros( stateFcn.nT, stateFcn.nH );
            stateFcn.M(1,:) = k(:,2);
            stateFcn.M(2,:) = k(:,2);
            stateFcn.M(3,:) = k(:,2);
        end

        %mur is the set relative permeability
        %Ms is the saturation magnetization in T
        function stateFcn = MakeMH_Fe_const_mur( mur, Ms )
        
            mu0 = 4 * pi * 1e-7;
            
            nH = 100;
            
            H = 10.^(linspace( 0, log10(Ms*2/mu0), nH ) );
            H(1) = 0;
            %H = linspace( 0, Ms*2/mu0, nH );
        
            M = (mur-1) * H;
            
            M( M>Ms/mu0 ) = Ms/mu0;
            
            plot(H*mu0,M*mu0,'-o');
            
            out = zeros( nH, 2 );
            out(:,1) = H;
            out(:,2) = M;
            
            dlmwrite( ['Fe_mur_' num2str(mur) '_Ms_' num2str(Ms) '.txt'], out, 'delimiter','\t','precision','%15.7f');
            
            
            stateFcn.nT = int32(3);
            stateFcn.nH = int32(nH);
            stateFcn.T = [200,300,400];
            stateFcn.H = out(:,1);
            stateFcn.M = zeros( stateFcn.nT, stateFcn.nH );
            stateFcn.M(1,:) = out(:,2);
            stateFcn.M(2,:) = out(:,2);
            stateFcn.M(3,:) = out(:,2);
        end


       %returns the volume of each tile in an array of dimensions n x 1
       %with n = length(tiles)
       function dV = getVolume( tiles )
           n = length(tiles);
           
           dV = zeros(n,1);
           
           for i=1:n
               switch tiles(i).tileType
                   case getMagTileType('prism')
                       dV(i) = prod( tiles(i).abc );
                   case getMagTileType('tetrahedron')
                       v = tiles(i).vertices;
                       dV(i) = 1./6. * abs( dot( v(:,1)-v(:,4) , cross( (v(:,2)-v(:,4)), (v(:,3)-v(:,4)) ) ) );
                       %the other tile types also need to be implemented
                       %here
               end
           end
       end
       
       %finds and returns the center position of each tile given the tile
       %type
       function C = getCenter( tiles )
           n = length(tiles);
           C = zeros(n,3);
           
           for i=1:n
               switch tiles(i).tileType
                   case getMagTileType('prism')
                       C(i,:) = tiles(i).offset;
                   case getMagTileType('tetrahedron')
                       C(i,:) = mean( tiles(i).vertices, 2 );
               end
           end
       end
       
       
   end
    
    
end
