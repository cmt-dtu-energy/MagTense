
%By Kaspar K. Nielsen, kasparkn@gmail.com
%First version on 12 May 2020
%Goes through all tiles and finds their respective neighbors and
%subsequently sets up their internal boundary conditions
%takes tiles as an input array and adds the struct-array
%bdryCdts with one element per boundary condition on the followin form:
%struct( 'Type', MagTenseTransientGeometry.FC_INTERNAL, 'l', d/2, 'A', d*d, 'n_ind', 2 );
%thus taking the type (as defined in MagTenseTransientGeometry)
%Distance l
%from center of the current tile to the center of the shared face (or
%external face in case of an external boundary condition)
%Surface area A of the shared face or the surface area towards the external
%boundary condition
%The index into the tile array, n_ind, of the neighbor with whom this
%(internal) boundary condition is shared. 
function tiles = MagtenseApplyBoundaryConditions( tiles )

n = length(tiles);
%loop over each tile
for i=1:n
    %loop over all the remaining tiles
    for j=i+1:n
        %get the shared faces between the two tiles if any exist
        face = MagTenseGetSharedFaces(tiles(i),tiles(j));
        if ~isempty(face)
           
           tiles(i).bdryCdts(length(tiles(i).bdryCdts)+1) = struct( 'Type', MagTenseTransientGeometry.FC_INTERNAL, 'l', face.lA, 'A', face.A, 'n_ind', j, 'FaceID', face.ID_A, 'bdryFun', @dummyfun );
           
           tiles(j).bdryCdts(length(tiles(j).bdryCdts)+1) = struct( 'Type', MagTenseTransientGeometry.FC_INTERNAL, 'l', face.lB, 'A', face.A, 'n_ind', i, 'FaceID', face.ID_B, 'bdryFun', @dummyfun );
           
        end
    end
end

end