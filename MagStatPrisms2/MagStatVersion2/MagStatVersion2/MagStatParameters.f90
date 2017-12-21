module MagStatParameters
    
    implicit none 
    
    type MagStatStateFunction
        real,dimension(:),allocatable :: T,H
        real,dimension(:,:), allocatable :: M        
        integer :: nT,nH
    endtype MagStatStateFunction
    
    contains
end module MagStatParameters
    