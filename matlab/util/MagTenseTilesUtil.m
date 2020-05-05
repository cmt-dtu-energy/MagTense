
%defines a class for doing static evaluations of tiles
classdef MagTenseTilesUtil
   properties
       
   end
   
   methods (Static)
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
       
       
   end
    
    
end