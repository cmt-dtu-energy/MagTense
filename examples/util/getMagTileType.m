function type = getMagTileType( tileName )

switch ( tileName )
    case 'cylinder'
        type = int32( 1 );
    case 'prism'
        type = int32( 2 );
    case 'circpiece'
        type = int32( 3 );
    case 'circpieceinv'
        type = int32( 4 );
    case 'planarcoil'
        type = int32( 101 );
end

end