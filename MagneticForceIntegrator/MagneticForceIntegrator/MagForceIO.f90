#include "fintrf.h"  
module MagForceIO
    use IntegrationDataTypes
    use MagStat2GetSolution
    implicit none

    !::Data model for the integration of the Maxwell stress tensor
    type, extends( dataCollectionModelBase) :: magStatModel
        type(MagTile),dimension(:),allocatable :: tiles
        integer :: n_tiles
        
        !::Variables below are deprecated as of 8 January 2018 (will not be further supported in the future and may eventually be deleted)
        integer :: n,spacedim
        real,dimension(:,:),allocatable :: dims,pos,M,rotAngles
        real,dimension(:,:,:),allocatable :: rotMat,rotMatInv         
    end type magStatModel
            
    real,parameter :: mu0=4*pi*1e-7
    
    contains
    
    !::
    !::Loads the data from Matlab into a fortran struct (the surface over which to integrate
    !::
    subroutine loadMagForceSurf( prhs, surf )
    mwPointer, intent(in) :: prhs
    type(surf_carth),intent(inout) :: surf
    
    character(len=10),dimension(:),allocatable :: fieldnames            
    mwSize :: sx
    integer :: nfields
    real*8,dimension(3) :: rectDims
    mwPointer :: xPtr,yPtr,zPtr,coneAnglePtr,z0Ptr,coordPtr
    mwPointer :: mxGetField, mxGetPr
    
        call getMagForceSurfFieldnames( fieldnames, nfields )
        sx = 1
        xPtr = mxGetField(prhs,sx,fieldnames(1))
        yPtr = mxGetField(prhs,sx,fieldnames(2))
        zPtr = mxGetField(prhs,sx,fieldnames(3))
        coneAnglePtr = mxGetField(prhs,sx,fieldnames(4))
        z0Ptr = mxGetField(prhs,sx,fieldnames(5))
        coordPtr = mxGetField(prhs,sx,fieldnames(6))        
        
        !::Copy from Matlab
        sx = 2
        call mxCopyPtrToReal8(mxGetPr(xPtr), surf%x, sx )
        call mxCopyPtrToReal8(mxGetPr(yPtr), surf%y, sx )
        call mxCopyPtrToReal8(mxGetPr(zPtr), surf%z, sx )
        
        sx = 1
        call mxCopyPtrToReal8(mxGetPr(coneAnglePtr), surf%cone_angle, sx )
        call mxCopyPtrToReal8(mxGetPr(z0Ptr), surf%z0, sx )        
        call mxCopyPtrToInteger4(mxGetPr(coordPtr), surf%coord, sx )        
        
        !::Setup the parameters
        select case ( surf%coord )
            case ( coord_sys_cone )
                 !::integrate over conical surface
                surf%r(1) = 0
                !::top circle, i.e. the biggest of the two conical circles
                surf%r(2) = surf%x(2)
                !::bottom circle
                surf%r(3) = surf%x(1) 
                surf%theta(1) = 0
                surf%theta(2) = 2*pi        
                surf%n_surfaces = 5
            case ( coord_sys_carth )
                surf%n_surfaces = 6
            case (coord_sys_cyl )
                !::integrate over cylindrical surface
                surf%r(1) = 0
                surf%r(2) = sqrt( surf%x(2)**2 + surf%y(2)**2 ) 
      
                surf%theta(1) = 0
                surf%theta(2) = 2*pi
                surf%n_surfaces = 4
            end select
            
        
        deallocate(fieldnames)
    end subroutine loadMagForceSurf
    
    !::
    !::Returns an array with the names of the fields expected in the tile struct
    subroutine getMagForceSurfFieldnames( fieldnames, nfields)
    integer,intent(out) :: nfields
    integer,parameter :: nf=9
    character(len=10),dimension(:),intent(out),allocatable :: fieldnames
            
        nfields = nf
        allocate(fieldnames(nfields))
        !::Setup the names of the members of the input struct
        fieldnames(1) = 'x'
        fieldnames(2) = 'y'
        fieldnames(3) = 'z'
        fieldnames(4) = 'cone_angle'
        fieldnames(5) = 'z0'
        fieldnames(6) = 'coord'
        fieldnames(7) = 'n_surfaces'        
        
        
        
    end subroutine getMagForceSurfFieldnames
    
end module MagForceIO
    