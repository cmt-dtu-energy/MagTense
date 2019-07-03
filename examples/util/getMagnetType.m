function type = getMagnetType( magnetName )

switch ( magnetName )
    case 'hard'
        type = int32( 1 );
    case 'soft'
        type = int32( 2 );
        
end

end