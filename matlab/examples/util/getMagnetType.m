function type = getMagnetType( magnetName )

switch ( magnetName )
    case 'hard'
        type = int32( 1 );
    case 'soft'
        type = int32( 2 );
    case 'soft_const_mur'
        type = int32( 3 );
        
end

end