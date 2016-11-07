MODULE parameters_CALL

    integer,parameter :: magTypeSoft=0,magTypeHard=1
    
    real,parameter :: pi=3.141592654
    
    type stateFunction
        real,dimension(:),allocatable :: T,H
        real,dimension(:,:), allocatable :: M
        real,dimension(:),allocatable :: Tcrit_cool
        real,dimension(:),allocatable :: Tcrit_heat
        real,dimension(:),allocatable :: Hcrit_cool
        real,dimension(:),allocatable :: Hcrit_heat
        integer :: nT,nH,nCritCool,nCritHeat
    endtype stateFunction
    
contains     
  subroutine writeDebugString( txt, valInt, valDbl )
    character(len=6) :: txt
    integer,optional :: valInt
    real, optional :: valDbl
    !open(11,file='debug.txt',status='old',access='sequential',form='formatted',position='append',action='write')
    !
    !
    !write(11,*) 'sjask'
    !if ( present( valint ) ) then
    !    write(11,*) txt,valint
    !else
    !    write(11,*) txt
    !endif
    !if ( present( valdbl ) ) then
    !    write(11,*) txt,valdbl
    !endif
    !
    !close(11)
    
  end subroutine writeDebugString
  
  subroutine writeDebugStringArr2D( txt, dblArr )
    character(len=6) :: txt    
    real, dimension(:,:),optional :: dblArr
    !open(11,file='debug.txt',status='old',access='sequential',form='formatted',position='append',action='write')
    !
    !
    !if ( present( dblArr ) ) then
    !    write(11,*) txt,dblArr
    !else
    !    write(11,*) txt
    !endif
    !
    !close(11)
    
  end subroutine writeDebugStringArr2D
  
  subroutine writeDebugStringArr1D( txt, dblArr )
    character(len=6) :: txt    
    real, dimension(:),optional :: dblArr
    !open(11,file='debug.txt',status='old',access='sequential',form='formatted',position='append',action='write')
    !
    !
    !if ( present( dblArr ) ) then
    !    write(11,*) txt,dblArr
    !else
    !    write(11,*) txt
    !endif
    !
    !close(11)
    
  end subroutine writeDebugStringArr1D
  
  
  subroutine writeDebugStringArr1DInt( txt, intArr )
    character(len=6) :: txt    
    integer*4, dimension(:),optional :: intArr
    !open(11,file='debug.txt',status='old',access='sequential',form='formatted',position='append',action='write')
    !
    !
    !if ( present( intArr ) ) then
    !    write(11,*) txt,intArr
    !else
    !    write(11,*) txt
    !endif
    !
    !close(11)
    !
    end subroutine writeDebugStringArr1DInt
  
    
END MODULE parameters_call
