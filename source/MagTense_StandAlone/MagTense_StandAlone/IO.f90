MODULE IO_CALL_Mag
    use TileNComponents
    use MagParameters
    !use MagTileIO
    implicit none
    
    type MagTileIOSetts
        real :: maxErr, T
        integer :: nIteMax
        logical :: iterateSolution
        type(MagStateFunction),dimension(1) :: stateFcn
        logical :: returnSolution
    end type MagTileIOSetts

    
    contains
    
    !> function for displaying output (to the terminal)
    !! @param err the current relative error
    !! @param err_max the current maximum allowed error
    function dispIte_fct( err, err_max )
    real,intent(in) :: err,err_max
    integer :: dispIte_fct
    
    write(*,*) err,err_max
    
    dispIte_fct = 0
    
    end function dispIte_fct
    
    !> Wrapper function for loading tiles from a file of any supported type.
    !! However, we only support a simple ASCII file format at the moment
    !! @param tiles array of tiles to be returned
    !! @param file Absolute file with the tile information
    !! @param mode 0 for read, 1 for write to file
    subroutine loadTiles( tiles, file, mode )
    type(MagTile),intent(inout),dimension(:),allocatable :: tiles
    character(len=1000),intent(in) :: file
    integer, intent(in) :: mode
    
    integer :: n,i
    
    write(*,*) 'Loading tiles',mode
    !! open file as read and allocate array of tiles
    if ( mode .eq. 0 ) then
        open(11,file=file,status='old',access='sequential',form='formatted',action='read')    
        !! read in the number of tiles (should be first element in the file)
        read(11,*) n
        !! allocate tiles
        allocate(tiles(n))
    else
        open(12,file=file,status='unknown',access='sequential',form='formatted',action='write')    
        !! write the number of tiles
        n = size(tiles)
        write(12,*) n
    endif
    
    !! loop over all tiles
    do i=1,n
        
        !! Now go through each parameter of the tile struct and read/write it
        !! r0, theta0, z0, dr, dtheta, dz
        if ( mode .eq. 0 ) then
            read(11,*) tiles(i)%r0,tiles(i)%theta0,tiles(i)%z0,tiles(i)%dr,tiles(i)%dtheta,tiles(i)%dz
        else
            write(12,*) tiles(i)%r0,tiles(i)%theta0,tiles(i)%z0,tiles(i)%dr,tiles(i)%dtheta,tiles(i)%dz
        endif
        
        !! abc
        if ( mode .eq. 0 ) then
            read(11,*) tiles(i)%a,tiles(i)%b,tiles(i)%c
        else
            write(12,*) tiles(i)%a,tiles(i)%b,tiles(i)%c
        endif
        
        !! M, u_ea, u_oa1, u_oa2
        if ( mode .eq. 0 ) then
            read(11,*) tiles(i)%M(1),tiles(i)%M(2),tiles(i)%M(3)
            read(11,*) tiles(i)%u_ea(1),tiles(i)%u_ea(2),tiles(i)%u_ea(3)
            read(11,*) tiles(i)%u_oa1(1),tiles(i)%u_oa1(2),tiles(i)%u_oa1(3)
            read(11,*) tiles(i)%u_oa2(1),tiles(i)%u_oa2(2),tiles(i)%u_oa2(3)
        else
            write(12,*) tiles(i)%M(1),tiles(i)%M(2),tiles(i)%M(3)
            write(12,*) tiles(i)%u_ea(1),tiles(i)%u_ea(2),tiles(i)%u_ea(3)
            write(12,*) tiles(i)%u_oa1(1),tiles(i)%u_oa1(2),tiles(i)%u_oa1(3)
            write(12,*) tiles(i)%u_oa2(1),tiles(i)%u_oa2(2),tiles(i)%u_oa2(3)
        endif
        
        !! Mrem, mu_r_ea, mu_r_oa
        if ( mode .eq. 0 ) then
            read(11,*) tiles(i)%Mrem,tiles(i)%mu_r_ea,tiles(i)%mu_r_oa
        else
            write(12,*) tiles(i)%Mrem, tiles(i)%mu_r_ea,tiles(i)%mu_r_oa
        endif
        
        !! tile type, magnet type, statefunction index, includeInIteration
        if ( mode .eq. 0 ) then
            read(11,*) tiles(i)%tileType,tiles(i)%magnetType,tiles(i)%stateFunctionIndex,tiles(i)%includeInIteration
        else
            write(12,*) tiles(i)%tileType,tiles(i)%magnetType,tiles(i)%stateFunctionIndex,tiles(i)%includeInIteration
        endif
        
        !!offset and rot angles
        if ( mode .eq. 0 ) then
            read(11,*) tiles(i)%offset(1),tiles(i)%offset(2),tiles(i)%offset(3)
            read(11,*) tiles(i)%rotAngles(1),tiles(i)%rotAngles(2),tiles(i)%rotAngles(3)
        else
            write(12,*) tiles(i)%offset(1),tiles(i)%offset(2),tiles(i)%offset(3)
            write(12,*) tiles(i)%rotAngles(1),tiles(i)%rotAngles(2),tiles(i)%rotAngles(3)
        endif
        
        
        !!use symmetry and symmetry operations and color
        if ( mode .eq. 0 ) then
            read(11,*) tiles(i)%exploitSymmetry,tiles(i)%symmetryOps(1),tiles(i)%symmetryOps(2),tiles(i)%symmetryOps(3)
            read(11,*) tiles(i)%color(1),tiles(i)%color(2),tiles(i)%color(3)
        else
            write(12,*)  tiles(i)%exploitSymmetry,tiles(i)%symmetryOps(1),tiles(i)%symmetryOps(2),tiles(i)%symmetryOps(3)
            write(12,*) tiles(i)%color(1),tiles(i)%color(2),tiles(i)%color(3)
        endif
        !!The relative change in M at the last iteration
        if ( mode .eq. 0 ) then
            read(11,*) tiles(i)%Mrel
        else
            write(12,*) tiles(i)%Mrel
        endif
        
        
        !<finally, setup the evaluation points
        if ( mode .eq. 0 ) then
            call setupEvaluationPoints( tiles(i) )
        endif
        
        
    enddo
    
    if ( mode .eq. 0 ) then
        close(11)
    else
        close (12)
    endif
            
    end subroutine loadTiles
    
    !> Saves the magnetic field values to a file (Hx,Hy,Hz)
    !! @param H the field to be written (n,3)
    !! @param file the absolute file to write to
    !! @param n the number of points 
    subroutine writeField( H, file, n )
    real,dimension(n,3),intent(in) :: H
    character(len=1000),intent(in) :: file
    integer,intent(in) :: n
    integer :: i
    open(12,file=file,status='unknown',access='sequential',form='formatted',action='write')    
    do i=1,n
        write(12,*) H(i,1), H(i,2), H(i,3)
    enddo
    
    close(12)
    end subroutine writeField
    
    !> Loads the few needed settings and returns a struct with these
    subroutine loadSettings( setts, file )
    type(MagTileIOSetts),intent(inout) :: setts
    character(len=1000) :: file
    integer :: i,tmp
        write(*,*) 'Loading settings'
        open(11,file=file,status='old',access='sequential',form='formatted',action='read')    
        
        read(11,*) setts%maxErr        
        read(11,*) setts%nIteMax        
        read(11,*) tmp
        if ( tmp .eq. 0 ) then
            setts%iterateSolution = .false.
        else
            setts%iterateSolution = .true.
        endif
        
        read(11,*) tmp
        if ( tmp .eq. 0 ) then
            setts%returnSolution = .false.
        else
            setts%returnSolution = .true.
        endif
        
        read(11,*) setts%T        
        read(11,*) setts%stateFcn(1)%nT,setts%stateFcn(1)%nH
        
        allocate( setts%stateFcn(1)%M(setts%stateFcn(1)%nT,setts%stateFcn(1)%nH) )
        allocate( setts%stateFcn(1)%T(setts%stateFcn(1)%nT) )
        allocate( setts%stateFcn(1)%H(setts%stateFcn(1)%nH) )
        
        read(11,*) setts%stateFcn(1)%T(1),setts%stateFcn(1)%T(1),setts%stateFcn(1)%T(2),setts%stateFcn(1)%T(3)
        
        do i=1,setts%stateFcn(1)%nH
             read(11,*) setts%stateFcn(1)%H(i),setts%stateFcn(1)%M(1,i),setts%stateFcn(1)%M(2,i),setts%stateFcn(1)%M(3,i)        
        enddo
        
        
        close (11)
        
        
    end subroutine loadSettings
    
    !> Loads the points at which the field is required from a file
    !! @param pts allocatable array on the form (n,3)
    !! @param file absolute file from which to get the points
    !! @param n_pts the number of points returned as an integer
    subroutine loadSolutionPoints( pts, file, n_pts )
    real,dimension(:,:),allocatable,intent(inout) :: pts
    character(len=1000),intent(in) :: file
    integer,intent(inout) :: n_pts
    integer :: i
    
    open(11,file=file,status='old',access='sequential',form='formatted',action='read')    
    read(11,*) n_pts
    allocate( pts(n_pts,3) )
    
    do i=1,n_pts
        read(11,*) pts(i,1), pts(i,2), pts(i,3)    
    enddo
    
    
    close(11)
    
    end subroutine loadSolutionPoints
    
END MODULE IO_CALL_Mag
    
