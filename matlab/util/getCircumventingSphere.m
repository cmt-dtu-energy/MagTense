
%based on the type of the tile a center and a radius are found and attached
%to the tile such as to define a circumventing sphere for easy checking
%whether two tiles can possibly be sharing faces
function circSph = getCircumventingSphere( tile )


%many pairs need to be implemented here eventually!
switch tile.tileType
   case getMagTileType('prism')
        C = tile.offset;
        R = sqrt( sum((tile.offset-tile.abc./2).^2 ) );
   case getMagTileType('cylinder')
   case getMagTileType('circpiece')
   case getMagTileType('tetrahedron')
       C = mean( tile.vertices, 2 )';
       R = max( sqrt( sum( ( repmat(C',1,4) - tile.vertices  ).^2, 1 ) ) );
end

circSph = struct('R',R,'C',C);
           

end