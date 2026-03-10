module MagParameters
    
    implicit none 
    
    type MagStateFunction
        real,dimension(:),allocatable :: T,H
        real,dimension(:,:), allocatable :: M
        real,dimension(:),allocatable :: y2a !! y2a is the spline derivates (returned from splin2)
        integer ( kind = 4 ) :: nT,nH
    endtype MagStateFunction
    
    type NStoreArr
        real,dimension(:,:,:,:),allocatable :: N    
    endtype NStoreArr
    
    
    contains
end module MagParameters
    