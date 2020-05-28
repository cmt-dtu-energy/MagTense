
function plotTileWithNeighbors( tiles, id )

tiles(id).color=[1,0,0];

tl = tiles(id);
for i=1:length(tiles(id).bdryCdts)
    ind = tiles(id).bdryCdts(i).n_ind;
    if ind>0
        tiles(ind).color = [0,0,1];
      tl = [tl tiles(ind)];
   end
end
plotTiles(tl);
end