

function [] = plotTiles( tiles, arrows )

if ~exist( 'arrows', 'var' )
    arrows = false;
end
for i=1:length(tiles)
    switch tiles(i).tileType
        case getMagTileType( 'cylinder' )
            %cylinder piece
            plotCylPiece( tiles(i));
        case getMagTileType( 'prism' )
            %a cube
            %cube_plot(tiles(i).offset,tiles(i).abc,tiles(i).rotAngles,3,tiles(i).color, tiles(i).graphRotxAng);
            cube_plot( tiles(i), [1,2,3] );
        case getMagTileType( 'circpiece' )
            plotCircPiece( tiles(i) );
        case getMagTileType( 'circpieceinv' )
            plotCircPieceInv( tiles(i) );
    end
end

if arrows
    plotFieldArrows( tiles );
end