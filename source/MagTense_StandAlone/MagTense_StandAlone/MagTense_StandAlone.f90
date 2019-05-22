!  MagTense_StandAlone.f90 
!
!  FUNCTIONS:
!  MagTense_StandAlone - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: MagTense_StandAlone
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program MagTense_StandAlone
    use TileNComponents
    use IterateMagnetSolution
    use IO_CALL_Mag
    implicit none

    character(len=1000) :: file_tiles_in,file_tiles_out,file_setup,file_field_pts,file_field_out
    type(MagTile),dimension(:),allocatable :: tiles
    type(MagTileIOSetts) :: setts
    integer :: n_tiles,n_pts
    real,dimension(:,:),allocatable :: pts,H
    real :: start,finish,resumeIteration
    procedure(displayIteration_fct),pointer :: disp_fct => null()
    call cpu_time(start)
    
    disp_fct => dispIte_fct
    
    !!no support for resuming iterations at the moment
    resumeIteration = 0
    
    !::Read the file named io.txt which should be located in the same folder as the executable
    !pause
    open(11,file='io.txt',status='old',access='sequential',form='formatted',action='read')
    read(11,*) file_tiles_in
    read(11,*) file_tiles_out
    read(11,*) file_setup
    read(11,*) file_field_pts
    read(11,*) file_field_out
    close(11)
    
    !::Load the input file containing all the tile information
    call loadTiles( tiles, file_tiles_in, 0 )
    n_tiles = size(tiles)
    !::Load the model settings file
    call loadSettings( setts, file_setup )
        
    !::If requested, save the iterated tiles
    if ( setts%iterateSolution ) then
        write(*,*) 'Doing iteration'
        call iterateMagnetization( tiles, n_tiles, setts%stateFcn, 1, setts%T, setts%maxErr, setts%nIteMax, disp_fct, resumeIteration )
        
        !< save the iterated tiles
        call loadTiles( tiles, file_tiles_out, 1 )
    endif
    
        
    !::If requested, save the found magnetic field
    !::Load the points where the solution is wanted
    if ( setts%returnSolution .eq. .true. ) then
        write(*,*) 'Finding solution at requested points'
        call loadSolutionPoints( pts, file_field_pts, n_pts )
        allocate(H(n_pts,3))
        call getFieldFromTiles( tiles, H, pts, n_tiles, n_pts )
        !< save the found field
        call writeField( H, file_field_out, n_pts )
    endif    
    call cpu_time(finish)
    
    write(*,*) 'Elapsed time', finish-start

    end program MagTense_StandAlone

