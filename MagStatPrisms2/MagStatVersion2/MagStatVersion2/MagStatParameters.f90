module MagStatParameters
    
    implicit none 
    
    type MagStatStateFunction
        real,dimension(:),allocatable :: T,H
        real,dimension(:,:), allocatable :: M        
        integer :: nT,nH
    endtype MagStatStateFunction
    
    type NStoreArr
        real,dimension(:,:,:,:),allocatable :: N    
    endtype NStoreArr
    
    
    contains
end module MagStatParameters
    