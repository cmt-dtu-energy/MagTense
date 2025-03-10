module DemagAuxFunctions

    use MicroMagParameters
    use DemagFieldGetSolution
    use iso_fortran_env
    use IO_GENERAL
    use UTIL_MICROMAG
    
    implicit none
    
    contains

    !>-----------------------------------------
    !> @author Rasmus Bjoerk, rabj@dtu.dk, DTU, 2022
    !> @brief
    !> Apply a threshold value to the dense demagnetization matrix that removes all absolute values smaller than epsilon
    !> It then converts the matrix to sparse matrix
    !> @params[in] problem the data structure containing the matrices
    !>-----------------------------------------
    subroutine ApplyThresholdDense( problem )
    type(MicroMagProblem),intent(inout) :: problem         !> Problem data structure  
    logical,dimension(:,:),allocatable :: mask_xx, mask_xy, mask_xz     !> mask used for finding non-zero values
    logical,dimension(:,:),allocatable :: mask_yy, mask_yz, mask_zz     !> mask used for finding non-zero values
    real(SP) :: threshold_var
    real(SP),dimension(:),allocatable  :: Kxx_abs, Kxy_abs, Kxz_abs, Kyy_abs, Kyz_abs, Kzz_abs  !> Temporary matrices with absolute values of the demag tensor, for threshold calculations
    integer ::  nx_K, ny_K, k_xx, k_xy, k_xz, k_yy, k_yz, k_zz 
    integer :: i,j,k,nx,ny,nz,ntot,ind                            !> Internal counters and index variables

    threshold_var =  problem%demag_threshold

    !If a fraction of entries is specified that are to be removed, find the function value for this
    !The demag tensor is considered as a whole, and the fraction specified concern the number of elements greater than epsilon
    if ( problem%demag_approximation .eq. DemagApproximationThresholdFraction ) then
            
        !Make a mask for each demag tensor with only the elements larger than zero
        !Make the masks only once - this is memory intensive, but computationally efficient
        nx_K = size(problem%Kxx(:,1))
        ny_K = size(problem%Kxx(1,:))
        allocate(mask_xx(nx_K,ny_K))
        allocate(mask_xy(nx_K,ny_K))
        allocate(mask_xz(nx_K,ny_K))
        allocate(mask_yy(nx_K,ny_K))
        allocate(mask_yz(nx_K,ny_K))
        allocate(mask_zz(nx_K,ny_K))
            
        mask_xx = problem%Kxx .gt. 0
        mask_xy = problem%Kxy .gt. 0
        mask_xz = problem%Kxz .gt. 0
        mask_yy = problem%Kyy .gt. 0
        mask_yz = problem%Kyz .gt. 0
        mask_zz = problem%Kzz .gt. 0
            
        !Make a copy (pack) of the tensors for speed so that the mask only has to be applied once
        allocate(Kxx_abs(count( mask_xx )))
        allocate(Kxy_abs(count( mask_xy )))
        allocate(Kxz_abs(count( mask_xz )))
        allocate(Kyy_abs(count( mask_yy )))
        allocate(Kyz_abs(count( mask_yz )))
        allocate(Kzz_abs(count( mask_zz )))
    
        !Pack the demag tensors into the temporary arrays. The intrinsic PACK routine overflows the stack, so use do loops
        k_xx = 1
        k_xy = 1
        k_xz = 1
        k_yy = 1
        k_yz = 1
        k_zz = 1
        do i=1,nx_K
            do j=1,ny_K
                if (mask_xx(i,j) .eqv. .true.) then
                    Kxx_abs(k_xx) = problem%Kxx(i,j) 
                    k_xx = k_xx + 1
                endif
                if (mask_xy(i,j) .eqv. .true.) then
                    Kxy_abs(k_xy) = problem%Kxy(i,j) 
                    k_xy = k_xy + 1
                endif
                if (mask_xz(i,j) .eqv. .true.) then
                    Kxz_abs(k_xz) = problem%Kxz(i,j) 
                    k_xz = k_xz + 1
                endif
                if (mask_yy(i,j) .eqv. .true.) then
                    Kyy_abs(k_yy) = problem%Kyy(i,j) 
                    k_yy = k_yy + 1
                endif
                if (mask_yz(i,j) .eqv. .true.) then
                    Kyz_abs(k_yz) = problem%Kyz(i,j) 
                    k_yz = k_yz + 1
                endif
                if (mask_zz(i,j) .eqv. .true.) then
                    Kzz_abs(k_zz) = problem%Kzz(i,j) 
                    k_zz = k_zz + 1
                endif
            enddo
        enddo
            
        call FindThresholdFraction(Kxx_abs, Kxy_abs, Kxz_abs, Kyy_abs, Kyz_abs, Kzz_abs, threshold_var)
            
        !cleanup
        deallocate(mask_xx)
        deallocate(mask_xy)
        deallocate(mask_xz)
        deallocate(mask_yy)
        deallocate(mask_yz)
        deallocate(mask_zz)
        deallocate(Kxx_abs)
        deallocate(Kxy_abs)
        deallocate(Kxz_abs)
        deallocate(Kyy_abs)
        deallocate(Kyz_abs)
        deallocate(Kzz_abs)
            
    endif
        
    call ConvertDenseToSparse_s( problem%Kxx, problem%K_s(1), threshold_var)
    call ConvertDenseToSparse_s( problem%Kxy, problem%K_s(2), threshold_var)
    call ConvertDenseToSparse_s( problem%Kxz, problem%K_s(3), threshold_var)
    call ConvertDenseToSparse_s( problem%Kyy, problem%K_s(4), threshold_var)
    call ConvertDenseToSparse_s( problem%Kyz, problem%K_s(5), threshold_var)
    call ConvertDenseToSparse_s( problem%Kzz, problem%K_s(6), threshold_var)
    
    end subroutine ApplyThresholdDense
    
    
    !>-----------------------------------------
    !> @author Rasmus Bjoerk, rabj@dtu.dk, DTU, 2022
    !> @brief
    !> Apply the fast fourier transform, then remove values below a certain threshold value to the FFT demagnetization 
    !> matrix that removes all absolute values smaller than epsilon. Then converts the matrix to a sparse matrix
    !> @params[in] problem the data structure containing the matrices
    !>-----------------------------------------
    subroutine ApplyThresholdFFT( problem )
    type(MicroMagProblem),intent(inout) :: problem         !> Problem data structure    
    logical,dimension(:,:),allocatable :: mask_xx, mask_xy, mask_xz     !> mask used for finding non-zero values
    logical,dimension(:,:),allocatable :: mask_yy, mask_yz, mask_zz     !> mask used for finding non-zero values
    integer ::  nx_K, ny_K, k_xx, k_xy, k_xz, k_yy, k_yz, k_zz 
    complex(kind=4) :: alpha_c, beta_c
    integer,dimension(3) :: L                                     !> Array specifying the dimensions of the fft
    complex(kind=4) :: thres
    real(SP),dimension(:),allocatable  :: Kxx_abs, Kxy_abs, Kxz_abs, Kyy_abs, Kyz_abs, Kzz_abs  !> Temporary matrices with absolute values of the demag tensor, for threshold calculations
    complex(kind=4),dimension(:,:),allocatable :: Kxx_c, Kxy_c, Kxz_c, Kyy_c, Kyz_c, Kzz_c      !> Temporary matrices for storing the complex version of the demag matrices
    integer :: status
    type(DFTI_DESCRIPTOR), POINTER :: desc_handle                       !> Handle for the FFT MKL stuff
    complex(kind=4),dimension(:,:),allocatable :: eye,FT,IFT,temp,temp2 !> Identity matrix, the indentity matrix' fourier transform and its inverse fourier transform
    integer :: i,j,k,nx,ny,nz,ntot,ind                            !> Internal counters and index variables
    real(SP) :: threshold_var
    
    !Apply the fast fourier transform concept, remove values below a certain threshold and then convert to a sparse matrix
    !for use in the field calculation later on
        
    !We have the following:
    ! H = N * M => 
    ! H = IFT * FT * N * IFT * FT * M
    ! with FT = FFT( eye(n,n) ) and IFT = IFFT( eye(n,n) )
    ! The matrix N' = FT * N * IFT is then reduced to a sparse matrix
    ! through N'(N'<demag_threshold) = 0 and then N' = sparse(N')
        
    !create the FT and IFT matrices
    allocate(eye(ntot,ntot),FT(ntot,ntot),IFT(ntot,ntot))
        
    !make the identity matrix
    eye(:,:) = 0
    do i=1,ntot
        eye(i,i) = cmplx(1.,0)
    enddo
        
    !get the FFT and the IFFT of the identity matrix
    L(1) = nx
    L(2) = ny
    L(3) = nz
    status = DftiCreateDescriptor( desc_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, L )
    !Set the descriptor to not override the input array (change of a default setting)
    status = DftiSetValue( desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
                
    status = DftiCommitDescriptor( desc_handle )
        
    do i=1,ntot
        status = DftiComputeForward( desc_handle, eye(:,i), FT(:,i) )            
    enddo
        
    !Find the inverse of the fourier transform
    allocate(temp(ntot,ntot) )
    temp = conjg(FT)
    IFT = transpose(temp) / ntot
        
    !do the change of basis of each of the demag tensor matrices
    !And make room for the very temporary complex output matrices
    allocate( Kxx_c(ntot,ntot), Kxy_c(ntot,ntot), Kxz_c(ntot,ntot) )
    allocate( Kyy_c(ntot,ntot), Kyz_c(ntot,ntot), Kzz_c(ntot,ntot) )
        
    alpha_c = cmplx(1.0,0.0)
    beta_c  = cmplx(0.0,0.0)
    Kxx_c = FT
    temp = cmplx(problem%Kxx)
        
    !temporary array to hold the results of the first calculation
    allocate(temp2(ntot,ntot) )

    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, FT   , ntot, temp, ntot, beta_c, temp2, ntot)
    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, temp2, ntot, IFT , ntot, beta_c, Kxx_c, ntot)
       
    Kxy_c = FT
    temp = cmplx(problem%Kxy)
    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, FT   , ntot, temp, ntot, beta_c, temp2, ntot)
    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, temp2, ntot, IFT , ntot, beta_c, Kxy_c, ntot)
        
    Kxz_c = FT
    temp = cmplx(problem%Kxz)
    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, FT   , ntot, temp, ntot, beta_c, temp2, ntot)
    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, temp2, ntot, IFT , ntot, beta_c, Kxz_c, ntot)
        
    Kyy_c = FT
    temp = cmplx(problem%Kyy)
    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, FT   , ntot, temp, ntot, beta_c, temp2, ntot)
    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, temp2, ntot, IFT , ntot, beta_c, Kyy_c, ntot)
        
    Kyz_c = FT
    temp = cmplx(problem%Kyz)
    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, FT   , ntot, temp, ntot, beta_c, temp2, ntot)
    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, temp2, ntot, IFT , ntot, beta_c, Kyz_c, ntot)
        
    Kzz_c = FT
    temp = cmplx(problem%Kzz)
    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, FT   , ntot, temp, ntot, beta_c, temp2, ntot)
    call cgemm('N', 'N', ntot, ntot, ntot, alpha_c, temp2, ntot, IFT , ntot, beta_c, Kzz_c, ntot)
    !Kxx_c = matmul( matmul( FT, problem%Kxx ), IFT )
    !Kxy_c = matmul( matmul( FT, problem%Kxy ), IFT )
    !Kxz_c = matmul( matmul( FT, problem%Kxz ), IFT )
    !Kyy_c = matmul( matmul( FT, problem%Kyy ), IFT )
    !Kyz_c = matmul( matmul( FT, problem%Kyz ), IFT )
    !Kzz_c = matmul( matmul( FT, problem%Kzz ), IFT )
    deallocate(temp)
    deallocate(temp2)
        
    threshold_var = problem%demag_threshold
    !If a fraction of entries is specified that are to be removed, find the function value for this
    !The demag tensor is considered as a whole, and the fraction specified concern the number of elements greater than epsilon
    if ( problem%demag_approximation .eq. DemagApproximationFFTThresholdFraction ) then
            
        !Make a mask for each demag tensor with only the elements larger than zero
        !Make the masks only once - this is memory intensive, but computationally efficient
        nx_K = size(Kxx_c(:,1))
        ny_K = size(Kxx_c(1,:))
        allocate(mask_xx(nx_K,ny_K))
        allocate(mask_xy(nx_K,ny_K))
        allocate(mask_xz(nx_K,ny_K))
        allocate(mask_yy(nx_K,ny_K))
        allocate(mask_yz(nx_K,ny_K))
        allocate(mask_zz(nx_K,ny_K))
            
        mask_xx = abs(Kxx_c) .gt. 0
        mask_xy = abs(Kxy_c) .gt. 0
        mask_xz = abs(Kxz_c) .gt. 0
        mask_yy = abs(Kyy_c) .gt. 0
        mask_yz = abs(Kyz_c) .gt. 0
        mask_zz = abs(Kzz_c) .gt. 0
            
        !Make a copy (pack) of the tensors for speed so that the mask only has to be applied once
        allocate(Kxx_abs(count( mask_xx )))
        allocate(Kxy_abs(count( mask_xy )))
        allocate(Kxz_abs(count( mask_xz )))
        allocate(Kyy_abs(count( mask_yy )))
        allocate(Kyz_abs(count( mask_yz )))
        allocate(Kzz_abs(count( mask_zz )))
    
        !Pack the demag tensors into the temporary arrays. The intrinsic PACK routine overflows the stack, so use do loops
        k_xx = 1
        k_xy = 1
        k_xz = 1
        k_yy = 1
        k_yz = 1
        k_zz = 1
        do i=1,nx_K
            do j=1,ny_K
                if (mask_xx(i,j) .eqv. .true.) then
                    Kxx_abs(k_xx) = abs(Kxx_c(i,j)) 
                    k_xx = k_xx + 1
                endif
                if (mask_xy(i,j) .eqv. .true.) then
                    Kxy_abs(k_xy) = abs(Kxy_c(i,j))
                    k_xy = k_xy + 1
                endif
                if (mask_xz(i,j) .eqv. .true.) then
                    Kxz_abs(k_xz) = abs(Kxz_c(i,j))
                    k_xz = k_xz + 1
                endif
                if (mask_yy(i,j) .eqv. .true.) then
                    Kyy_abs(k_yy) = abs(Kyy_c(i,j)) 
                    k_yy = k_yy + 1
                endif
                if (mask_yz(i,j) .eqv. .true.) then
                    Kyz_abs(k_yz) = abs(Kyz_c(i,j))
                    k_yz = k_yz + 1
                endif
                if (mask_zz(i,j) .eqv. .true.) then
                    Kzz_abs(k_zz) = abs(Kzz_c(i,j)) 
                    k_zz = k_zz + 1
                endif
            enddo
        enddo
            
        call FindThresholdFraction(Kxx_abs, Kxy_abs, Kxz_abs, Kyy_abs, Kyz_abs, Kzz_abs, threshold_var)
            
        !cleanup
        deallocate(mask_xx)
        deallocate(mask_xy)
        deallocate(mask_xz)
        deallocate(mask_yy)
        deallocate(mask_yz)
        deallocate(mask_zz)
        deallocate(Kxx_abs)
        deallocate(Kxy_abs)
        deallocate(Kxz_abs)
        deallocate(Kyy_abs)
        deallocate(Kyz_abs)
        deallocate(Kzz_abs)
            
    endif
       
    !Apply the thresholding
    thres = cmplx(threshold_var,0.)
    call ConvertDenseToSparse_c( Kxx_c, problem%K_s_c(1), thres)
    call ConvertDenseToSparse_c( Kxy_c, problem%K_s_c(2), thres)
    call ConvertDenseToSparse_c( Kxz_c, problem%K_s_c(3), thres)
    call ConvertDenseToSparse_c( Kyy_c, problem%K_s_c(4), thres)
    call ConvertDenseToSparse_c( Kyz_c, problem%K_s_c(5), thres)
    call ConvertDenseToSparse_c( Kzz_c, problem%K_s_c(6), thres)
        
    !then use problem%K_s(1..6) with cuda or MKL to do the sparse matrix-vector product with the FFT(M) and subsequently the IFT on the whole thing to get H
        
    !clean-up
    status = DftiFreeDescriptor(desc_handle)
    deallocate(Kxx_c,Kxy_c,Kxz_c,Kyy_c,Kyz_c,Kzz_c)
    deallocate(eye,FT,IFT)
    !deallocate(problem%Kxx,problem%Kxy,problem%Kxz)
    !deallocate(problem%Kyy,problem%Kyz,problem%Kzz)
        
    !Make descriptor handles for the FFT of M and IFFT of H
    status = DftiCreateDescriptor( problem%desc_hndl_FFT_M_H, DFTI_SINGLE, DFTI_COMPLEX, 3, L )
        
    status = DftiCommitDescriptor( problem%desc_hndl_FFT_M_H )
    
    end subroutine ApplyThresholdFFT
    
    !>-----------------------------------------
    !> @author Rasmus Bjoerk, rabj@dtu.dk, DTU, 2020
    !> @brief
    !> Find the threshold value that corresponds to a certain fraction of the non-zero elements in the demag tensor
    !> Its assumed that the absolute value of the matrix is passed. This means that the subroutine will work for both single and complex tensors 
    !> Do this using bisection algorithm
    !> @params[in] threshold a faction specifying the fraction of the smallest elements that are to be removed
    !>-----------------------------------------
    subroutine FindThresholdFraction(Kxx, Kxy, Kxz, Kyy, Kyz, Kzz, threshold_var)
    real(SP),dimension(:),intent(in) :: Kxx, Kxy, Kxz, Kyy, Kyz, Kzz                !> The absolute of the demag tensors 
    real(SP),intent(inout) :: threshold_var
    real(SP) :: f_large, f_small, f_middle
    integer,dimension(6) :: count_ind
    integer :: count_middle, n_ele_nonzero, k_do  
    character*(10) :: prog_str
    
    if (threshold_var .ge. 1) then ! special case. 
        call displayGUIMessage( 'Threshold fraction >= 1 selected. Setting demagnetization to zero.' )
        !f_large = max(maxval(Kxx), maxval(Kxy), maxval(Kxz), maxval(Kyy), maxval(Kyz), maxval(Kzz))
        !2*f_large because f_large sometimes leaves elements in demag tensors, possibly due to numerical errors.
        !threshold_var = 2*f_large
        threshold_var = huge(threshold_var)
    else
        call displayGUIMessage( 'Starting threshold calculation.' )
    
        !The total number of nonzero elements in the demag tensor
        n_ele_nonzero = size( Kxx )
        n_ele_nonzero = n_ele_nonzero + size( Kxy )
        n_ele_nonzero = n_ele_nonzero + size( Kxz )
        n_ele_nonzero = n_ele_nonzero + size( Kyy )
        n_ele_nonzero = n_ele_nonzero + size( Kyz )
        n_ele_nonzero = n_ele_nonzero + size( Kzz )
    
        !Find the maximum and minimum value in the individual demag tensors than is greater than zero
        f_large = max(maxval(Kxx), maxval(Kxy), maxval(Kxz), maxval(Kyy), maxval(Kyz), maxval(Kzz))    
        f_small = min(minval(Kxx), minval(Kxy), minval(Kxz), minval(Kyy), minval(Kyz), minval(Kzz))
    
        !Bisection algoritm for find the function value for a corresponding fraction
        k_do = 1
        do 
            f_middle = (f_large-f_small)/2+f_small
                
            !Count only the elements larger than epsilon in each of the matrices
            count_ind(1) = count( Kxx .le. f_middle )
            count_ind(2) = count( Kxy .le. f_middle )
            count_ind(3) = count( Kxz .le. f_middle )
            count_ind(4) = count( Kyy .le. f_middle )
            count_ind(5) = count( Kyz .le. f_middle )
            count_ind(6) = count( Kzz .le. f_middle )
                
            count_middle = sum(count_ind)
                   
            if ( count_middle .gt. threshold_var*n_ele_nonzero ) then
                f_large = f_middle
            else
                f_small = f_middle
            endif            
            
            !check if we have found the value that defines the threshold
            if ( ( count_middle .ge. (threshold_var-0.01)*n_ele_nonzero ) .and. ( count_middle .le. (threshold_var+0.01)*n_ele_nonzero ) ) then
                threshold_var = f_middle
                
                exit
            endif
            
            if ( k_do .gt. 1000 ) then
                call displayGUIMessage( 'Iterations exceeded in finding threshold value. Stopping iterations.' )
                
                threshold_var = f_middle
                exit
            endif
                
            k_do = k_do+1
        enddo  
    
        call displayGUIMessage( 'Using a threshold value of :' )
        write (prog_str,'(F10.9)') threshold_var
        call displayGUIMessage( prog_str )
        call displayGUIMessage( 'i.e. a fraction of:' )
        write (prog_str,'(F6.4)') real(count_middle)/real(n_ele_nonzero)
        call displayGUIMessage( prog_str )
    endif
    
    end subroutine FindThresholdFraction
    
    
    !>-----------------------------------------
    !> @author Rasmus Bj�rk, rabj@dtu.dk, DTU, 2024
    !> @brief
    !> Change the demag field based on a error drawn from a standard distribution  
    !> @param[inout] problem the data structure containing the problem
    !---------------------------------------------------------------------------   
    subroutine AddUncertaintyToDemagField( problem, solution)
    type(MicroMagProblem),intent(in) :: problem         !> Problem data structure    
    type(MicroMagSolution),intent(inout) :: solution    !> Solution data structure
    
    integer :: nx,ny,nz,ntot
    
        nx = problem%grid%nx
        ny = problem%grid%ny
        nz = problem%grid%nz
        ntot = nx * ny * nz
        
        !For each field value, use this as the mean for a normal distribution and draw random numbers from this
        !Use the Box-Muller transformation to generate the random numbers
        allocate(solution%u1(ntot), solution%u2(ntot),solution%u3(ntot), solution%u4(ntot),solution%u5(ntot), solution%u6(ntot))
        call random_number(solution%u1)
        call random_number(solution%u2)
        call random_number(solution%u3)
        call random_number(solution%u4)
        call random_number(solution%u5)
        call random_number(solution%u6)
        solution%u1 = 1d0 - solution%u1
        solution%u2 = 1d0 - solution%u2
        solution%u3 = 1d0 - solution%u3
        solution%u4 = 1d0 - solution%u4
        solution%u5 = 1d0 - solution%u5
        solution%u6 = 1d0 - solution%u6
        
    end subroutine AddUncertaintyToDemagField
    
    
    !>-----------------------------------------
    !> @author Rasmus Bj�rk, rabj@dtu.dk, DTU, 2024
    !> @brief
    !> Reduce the size of the demag tensor by figuring out which tiles will have the same demag tensor 
    !> @param[inout] problem the data structure containing the problem
    !---------------------------------------------------------------------------   
    subroutine ConstructDemagTensorMap( problem)
    type(MicroMagProblem),intent(inout) :: problem          !> Problem data structure    
    integer, dimension(:), allocatable :: hash_map          !> Map of hash values to corresponding 3x3 demag tensors
    logical,dimension(:,:),allocatable :: globalHash_sign   !> The sign of the respective 3x3 demag tensors
    integer, dimension(:), allocatable :: hash_arr          !> List of computed hash values for a single tile to different evaluation points
    logical,dimension(:,:),allocatable :: sign_arr          !> The sign of the respective 3x3 demag tensors computed for a single tile to different evaluation points
    type(MagTile),dimension(1) :: tile                      !> Tile representing the current tile under consideration
    integer :: nx, ny, nz, ntot, index, i, j, k, m
    real(DP), dimension(:,:),allocatable :: pts_arr         !> The list of points being evaluated
    integer(kind=int32) :: i32
    character*(100) :: prog_str
    
        nx = problem%grid%nx
        ny = problem%grid%ny
        nz = problem%grid%nz
        ntot = nx * ny * nz             ! Number of tiles
        m = size(problem%grid%pts(:,1)) ! Number of evaluation points. Here the same as ntot
        
        ! The points at which the demag tensor is evaluated
        allocate(pts_arr(m,3))
        pts_arr(:,1) =  problem%grid%pts(:,1)
        pts_arr(:,2) =  problem%grid%pts(:,2)
        pts_arr(:,3) =  problem%grid%pts(:,3)
        
        ! Allocate the respective arrays used
        allocate(hash_map(huge(i32)))   !> A very simply hash_map comprising all possible integers is used. This is memory intensive, but fast
        hash_map(:) = 0
        
        allocate(globalHash_sign(ntot*m,3))
        globalHash_sign(:,:) = .false.
        
        allocate(problem%tensorMap(ntot,m))     !> The tensor map
        problem%tensorMap(:,:) = 0
        
        allocate(problem%tensorMapX(ntot,m))    !> The sign of the respective 3x3 demagnetization tensors x-components
        problem%tensorMapX(:,:) = .false.
        
        allocate(problem%tensorMapY(ntot,m))    !> The sign of the respective 3x3 demagnetization tensors y-components
        problem%tensorMapY(:,:) = .false.
        
        allocate(problem%tensorMapZ(ntot,m))    !> The sign of the respective 3x3 demagnetization tensors z-components
        problem%tensorMapZ(:,:) = .false.
        
        allocate(hash_arr(ntot))
        allocate(sign_arr(ntot,3))
       
        ! For each tile find the hashes related to each point
        index = 0
        do i=1,ntot
            
            !write(prog_str,'(A10, I7.2, A8, I7.2, A6, F6.2, A7)') 'Tile nr.: ', i, ' out of ', ntot, ' i.e. ', real(i)/real(ntot)*100,'% done'
            !call displayGUIMessage( trim(prog_str) )
                        
            ! Setup the current tile that we are looking at
            tile(1)%tileType = 2 !(for prism)
            tile(1)%a = problem%grid%abc(i,1)
            tile(1)%b = problem%grid%abc(i,2)
            tile(1)%c = problem%grid%abc(i,3)
            tile(1)%offset(1) = problem%grid%pts(i,1)
            tile(1)%offset(2) = problem%grid%pts(i,2)
            tile(1)%offset(3) = problem%grid%pts(i,3)
            
            ! Get the hash values for all evaluations points for the given tile
            call getHashFromTile(tile, pts_arr, hash_arr, sign_arr);
                                 
            ! For each hash ( (tile,pts) set) figure out if it has been treated previously
            do j=1,m
                if (hash_map(hash_arr(j)) /= 0) then
                    ! Then set tensorMap(j,i) to the relevant already calculated tensor
                    k = hash_map(hash_arr(j))
                    problem%tensorMap(i,j) = k
                    
                    ! Check the signs of the different 3x3 tensor components
                    if (sign_arr(j,1) /= globalHash_sign(k,1)) then
                        problem%tensorMapX(i,j) = .true.
                    endif
                    if (sign_arr(j,2) /= globalHash_sign(k,2)) then
                        problem%tensorMapY(i,j) = .true.
                    endif
                    if (sign_arr(j,3) /= globalHash_sign(k,3)) then
                        problem%tensorMapZ(i,j) = .true.
                    endif
                else 
                    ! If the hash value does not exist then add it to the tensor and the hash_map
                    index = index + 1
                    hash_map(hash_arr(j)) = index;
                    globalHash_sign(index,:) = sign_arr(j,:);

                    problem%tensorMap(i,j) = index;
                endif
            enddo
        enddo
        
        ! Deallocate
        deallocate(pts_arr)
        deallocate(hash_map)
        deallocate(globalHash_sign)
        deallocate(hash_arr)
        deallocate(sign_arr)
        
        !problem%tensorMap = Transpose(problem%tensorMap)
        !problem%tensorMapX = Transpose(problem%tensorMapX)
        !problem%tensorMapY = Transpose(problem%tensorMapY)
        !problem%tensorMapZ = Transpose(problem%tensorMapZ)
        
        !open(21,file='tensorMap.txt',status='unknown',form='formatted',action='write')
        !open(22,file='tensorMapX.txt',status='unknown',form='formatted',action='write')
        !open(23,file='tensorMapY.txt',status='unknown',form='formatted',action='write')
        !open(24,file='tensorMapZ.txt',status='unknown',form='formatted',action='write')
        !do i=1,ntot
        !    do j=1,m
	    !        write(21,*)  problem%tensorMap(i,j)
        !        write(22,*)  problem%tensorMapX(i,j)
        !        write(23,*)  problem%tensorMapY(i,j)
        !        write(24,*)  problem%tensorMapZ(i,j)
        !    enddo
        !enddo
        !close(21)
        !close(22)
        !close(23)
        !close(24)
        
    end subroutine ConstructDemagTensorMap
    
    
    !>-----------------------------------------
    !> @author Rasmus Bj�rk, rabj@dtu.dk, DTU, 2024
    !> @brief
    !> Calculate a hash array for each tile to all field evaluation points
    !> @param[inout] problem the data structure containing the problem
    !---------------------------------------------------------------------------   
    subroutine getHashFromTile( tile, pts, hash_arr, sign_arr)
        type(MagTile),dimension(1),intent(in) :: tile                       !> Tile representing the current tile under consideration
        real(DP), dimension(:,:),allocatable,intent(in) :: pts
        integer, dimension(:),intent(inout) :: hash_arr
        logical, dimension(:,:),intent(inout) :: sign_arr
        integer :: m, j, hash_single
        real(DP) :: scaling
        real(DP), dimension(6) :: hash_arr_temp
        
            m = size(pts(:,1));
            scaling = sqrt(tile(1)%a**2+tile(1)%b**2+tile(1)%c**2)
            hash_arr_temp(1) = tile(1)%a/scaling
            hash_arr_temp(2) = tile(1)%b/scaling
            hash_arr_temp(3) = tile(1)%c/scaling
                        
            do j=1,m
                hash_arr_temp(4:6)  = (tile(1)%offset(1:3)-pts(j,1:3))/scaling               
                call permutationVariantHash(abs(hash_arr_temp), hash_single)
                hash_arr(j) = hash_single
                if (hash_arr_temp(4) >= 0) then
                    sign_arr(j,1) = .true.
                else
                    sign_arr(j,1) = .false.
                endif
                if (hash_arr_temp(5) >= 0) then
                    sign_arr(j,2) = .true.
                else
                    sign_arr(j,2) = .false.
                endif
                if (hash_arr_temp(6) >= 0) then
                    sign_arr(j,3) = .true.
                else
                    sign_arr(j,3) = .false.
                endif
                !sign(sign_arr(j,1:3), int(hash_arr_temp(4:6)))
            enddo
            
    end subroutine getHashFromTile

    
    !>-----------------------------------------
    !> @author Rasmus Bj�rk, rabj@dtu.dk, DTU, 2024
    !> @brief
    !> Calculate a hash value for an array of floats
    !> @param[inout] problem the data structure containing the problem
    !---------------------------------------------------------------------------   
    subroutine permutationVariantHash( array, hash)
        
        real(DP),dimension(:),intent(in) :: array                      !> The array of floats from which to find a hash
        integer, intent(inout) :: hash
        integer :: p, i, large_num
        
            ! A prime number as the base of the polynomial
            p = 31
            
            ! A large number
            large_num = 1e9 + 9
            
            ! Initial hash value
            hash = 0
            
            ! Evaluate polynomial for each element
            do i = 1,size(array)
                hash = hash + array(i) * p**(i - 1)
            enddo
            
            ! Modulo operation to keep the hash value within a range
            hash = mod(hash, large_num)
        
    end subroutine permutationVariantHash
    
end module DemagAuxFunctions