module MagneticForce
    use IntegrationDataTypes
    use TileNComponents
    implicit none

    contains
    
    
    subroutine getForce( tiles, n_tiles,  surf )
    type( MagTile ), intent(in), dimension(n_tiles) :: tiles
    integer,intent(in) :: n_tiles
    type( surf_carth ), intent(in) :: surf
    
    end subroutine getForce
    
end module MagneticForce
    