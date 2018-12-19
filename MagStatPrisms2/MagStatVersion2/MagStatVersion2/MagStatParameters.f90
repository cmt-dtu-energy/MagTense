module MagStatParameters
    
    implicit none 
    
    type MagStatStateFunction
        real,dimension(:),allocatable :: T,H
        real,dimension(:,:), allocatable :: M
        real,dimension(:),allocatable :: y2a !! y2a is the spline derivates (returned from splin2)
        integer :: nT,nH
    endtype MagStatStateFunction
    
    type NStoreArr
        real,dimension(:,:,:,:),allocatable :: N    
    endtype NStoreArr
    
    
    contains
end module MagStatParameters
    