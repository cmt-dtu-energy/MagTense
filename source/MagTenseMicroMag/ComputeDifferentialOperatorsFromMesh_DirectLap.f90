module DifferentialOperators
  use MKL_SPBLAS
  use BLAS95
  use LAPACK95
  use MicroMagParameters
  use IO_GENERAL
  
  implicit none

contains

    subroutine computeDifferentialOperatorsFromMesh_DirectLap(GridInfo, interpn, weight, method, Aexch, ExtW)
        ! Subroutine to compute differential operators from mesh using Direct Laplacian method
        ! Inputs:
        ! GridInfo - Contains mesh information
        ! interpn - Interpolation scheme
        ! weight - Weighting scheme
        ! method - Calculation method
        ! Aexch - Exchange interaction strength
        ! ExtW - Optional external weights
        ! Outputs:
        ! DX, DY, DZ - Sparse matrices for derivatives
        ! W - Sparse matrix for averages over faces

        type(MicroMagGridInfo), intent(inout) :: GridInfo
        character(len=*), intent(inout), optional :: interpn
        real(dp), intent(inout), optional :: weight
        character(len=*), intent(inout), optional :: method
        real(dp), dimension(:), intent(in), optional :: Aexch
        real(dp), dimension(:,:), intent(in), optional :: ExtW
        !real(dp), dimension(:,:), intent(out) :: DX, DY, DZ
        !real(dp), dimension(:,:), intent(out) :: W
        type(sparse_matrix_t) :: W_matrix, FX_matrix, DDXA_matrix, DX_matrix

        ! Unpack the variables
        real(dp), dimension(:),allocatable :: NX, NY, NZ, Areas, Volumes, Xel, Yel, Zel, Xf, Yf, Zf
        integer, dimension(:,:),allocatable :: Signs, T, D, el2fa
        real(dp), dimension(:), allocatable :: Aexch_local
        integer :: i,j,dims, n, k, kk, k_mask1D, counter, k_i, lind, sjask, indx, k_j, k_row, info
        real(dp), dimension(:),allocatable :: VolCoeff, AX, AY, AZ, s, w
        integer, dimension(:),allocatable :: ss, ns, ks, inds1, inds2, inds2_vals, ks_sorted, ns_sorted, ks_CSC_start, ks_CSC_end, ks_CSC_start_temp, ks_CSC_end_temp
        integer, dimension(:),allocatable :: ks_sorted_1, ns_sorted_1, ns_CSC_start, ns_CSC_end, ns_CSC_start_temp, ns_CSC_end_temp
        integer, dimension(:),allocatable :: ks_sorted_reduced, ns_sorted_reduced, ks_sorted_reduced_vx, ns_sorted_reduced_vx, ks_sorted_non_zeroes, ns_sorted_non_zeroes, ks_sorted_reduced_1, ns_sorted_reduced_1
        integer, dimension(:),allocatable :: ns_packed, k_log
        integer, dimension(:),allocatable :: ks_non_zeroes
        real(dp), dimension(:), allocatable :: Amat
        integer, dimension(:), allocatable :: tmp, tmp2, ind, IPIV
        logical, allocatable :: mask1D(:), mask_int2log(:), mask1D_2(:), mask1D_3(:)
        real(dp), dimension(:), allocatable :: ddxA, ddyA, ddzA, ddx, ddy, ddz, dx, dy, dz, vx, vy, vz, vw, dks, dxk, dyk, dzk, nns
        real(dp), dimension(:), allocatable :: ddxA_non_zeroes, vx_non_zeroes, ddxA_sorted, ddxA_sorted_reduced, vx_non_zeroes_reduced, vx_sorted_reduced_1
        real(dp), dimension(:), allocatable :: dxk2, dyk2, dzk2, Wk, Wk2, e, sjaskarr
        real(dp), dimension(:,:), allocatable :: Gkl1, Gk, Hk, Gkl1_temp, Gkl1_T, GkRed, HkRed, Wktmp
        real(dp) :: wm, infinity, scale, scale_local, eps_criteria
        real(dp), allocatable :: extra(:,:), DX_matrix_dense(:,:)
        integer :: stat                                   !> Status value for the various sparse matrix operations 
        character*(40) :: prog_str
        !integer, dimension(:), allocatable :: cols_start(:)
        integer :: rows, cols, nnz
        integer, pointer :: row_ind(:) => null()
        integer, pointer :: col_ind(:) => null()
        integer, pointer :: cols_start(:) => null()
        integer, pointer :: cols_end(:) => null()
        real(dp), pointer :: values(:) => null()
        integer, dimension(:),allocatable :: A_rows_CSC, A_cols_start_CSC, A_cols_end_CSC, B_rows_CSC, B_cols_start_CSC, B_cols_end_CSC
        integer, dimension(:),allocatable :: A_cols_CSR, A_rows_start_CSR, A_rows_end_CSR, B_cols_CSR, B_rows_start_CSR, B_rows_end_CSR
        real(dp), dimension(:),allocatable :: A_values_CSC, B_values_CSC, A_values_CSR, B_values_CSR
        type(sparse_matrix_t) :: descr
        type(sparse_matrix_t) :: A_sparse_CSC, A_sparse_CSR, C_sparse_CSC
        type(sparse_matrix_t) :: B_sparse_CSC, B_sparse_CSR, C_sparse_CSR
        type(sparse_matrix_t) :: Identity_matrix_K_sparse,Identity_matrix_N_sparse
        real(dp), dimension(:,:),allocatable :: C_dense_CSC, C_dense_CSR
        integer, dimension(:),allocatable :: Identity_matrix_K_CSC_start, Identity_matrix_K_CSC_end, Identity_matrix_K_rows, Identity_matrix_N_CSC_start, Identity_matrix_N_CSC_end, Identity_matrix_N_rows
        real(dp), dimension(:), allocatable :: Identity_matrix_K_values, Identity_matrix_N_values
        real(dp), dimension(:,:), allocatable :: DDXA_matrix_dense, FX_matrix_dense
        real(dp) :: alpha, beta
        
        
        
          !INTEGER          Nmkl, NRHS
          !PARAMETER        ( Nmkl = 5, NRHS = 3 )
          !INTEGER          LDA, LDB
          !PARAMETER        ( LDA = Nmkl, LDB = Nmkl )
          !INTEGER          INFO
          !INTEGER          IPIV( Nmkl )
          !real(dp), allocatable :: A( :, : ), B( :, : )
          
          !allocate(A(5,5))
          !A(1,:) = [6.80,  -6.05,  -0.45,   8.32,  -9.67]
          !A(2,:) = [-2.11,  -3.30,   2.58,   2.71,  -5.14]
          !A(3,:) = [5.66,   5.36,  -2.70,   4.35,  -7.26]
          !A(4,:) = [5.97,  -4.44,   0.27,  -7.17,   6.08]
          !A(5,:) = [8.23,   1.08,   9.04,   2.14,  -6.87]
          
          !allocate(B(5,3))
          !B(1,:) = [4.02,  -1.56,   9.81]
          !B(2,:) = [6.19,   4.00,  -4.09]
          !B(3,:) = [-8.22,  -8.67,  -4.57]
          !B(4,:) = [-7.57,   1.75,  -8.61]
          !B(5,:) = [-3.03,   2.86,   8.99]
          !B(1,:) = [4.02]
          !B(2,:) = [6.19]
          !B(3,:) = [-8.22]
          !B(4,:) = [-7.57]
          !B(5,:) = [-3.03]
         
          !CALL DGESV( Nmkl, NRHS, A, LDA, IPIV, B, LDB, INFO )
        
          !call displayGUIMessage( 'LAPACK' )
          !write (prog_str,'(I10)') (Nmkl)
          !call displayGUIMessage( prog_str )
          !write (prog_str,'(I10)') (NRHS)
          !call displayGUIMessage( prog_str )
          !write (prog_str,'(I10)') (LDA)
          !call displayGUIMessage( prog_str )
          !write (prog_str,'(I10)') (LDB)
          !call displayGUIMessage( prog_str )
          !write (prog_str,'(I10)') (INFO)
          !call displayGUIMessage( prog_str )
        
          !open(21,file='B.txt',status='unknown',form='formatted',action='write')
          !do i=1,size(B,2)
          !    write(21,*)  B(1,i)
          !    write(21,*)  B(2,i)
          !    write(21,*)  B(3,i)
          !    write(21,*)  B(4,i)
          !    write(21,*)  B(5,i)
          !enddo
          !close(21)
        
        allocate(A_values_CSC(13), A_rows_CSC(13), A_cols_start_CSC(5), A_cols_end_CSC(5))
        allocate(B_values_CSC(13), B_rows_CSC(13), B_cols_start_CSC(5), B_cols_end_CSC(5))
        A_values_CSC(:)	=	[1,	-2,	-4,	-1,	5,	8,	4,	2,	-3,	6,	7,	4,	-5]
        A_rows_CSC(:)	=	[1,	2,	4,	1,	2,	5,	3,	4,	1,	3,	4,	3,	5]
        A_cols_start_CSC(:)	=	[1,	4,	7,	9,	12]	 	 	 	 	 	 	 	 
        A_cols_end_CSC(:) =	[4,	7,	9,	12,	14]
        
        B_values_CSC(:)	=	[1,	-2,	-4,	-1,	5,	8,	4,	2,	-3,	6,	7,	4,	-5]
        B_rows_CSC(:)	=	[1,	2,	4,	1,	2,	5,	3,	4,	1,	3,	4,	3,	5]
        B_cols_start_CSC(:)	=	[1,	4,	7,	9,	12]	 	 	 	 	 	 	 	 
        B_cols_end_CSC(:) =	[4,	7,	9,	12,	14]
        
        N = 5
        K = 5
        stat = mkl_sparse_d_create_csc (A_sparse_CSC, SPARSE_INDEX_BASE_ONE, N, K, A_cols_start_CSC, A_cols_end_CSC, A_rows_CSC, A_values_CSC)
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        stat = mkl_sparse_d_create_csc (B_sparse_CSC, SPARSE_INDEX_BASE_ONE, N, K, B_cols_start_CSC, B_cols_end_CSC, B_rows_CSC, B_values_CSC)
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        allocate(C_dense_CSC(N,N))
        stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, A_sparse_CSC, B_sparse_CSC, SPARSE_LAYOUT_COLUMN_MAJOR, C_dense_CSC, N)
        
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        open(21,file='C_dense_CSC.txt',status='unknown',form='formatted',action='write')
        do i=1,size(C_dense_CSC,1)
            do j=1,size(C_dense_CSC,2)
                write(21,*)  C_dense_CSC(i,j)
            enddo
        enddo
        close(21)
       
        
        call displayGUIMessage( 'Starting sparse product' )
        
        stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, A_sparse_CSC, B_sparse_CSC, C_sparse_CSC)
        
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        call displayGUIMessage( 'Sparse product done' )
        
        
        
        allocate(A_values_CSR(13), A_cols_CSR(13), A_rows_start_CSR(5), A_rows_end_CSR(5))
        allocate(B_values_CSR(13), B_cols_CSR(13), B_rows_start_CSR(5), B_rows_end_CSR(5))
        A_values_CSR(:)	=	[1,	-1,	-3,	-2,	5,	4,	6,	4,	-4,	2,	7,	8,	-5]
        A_cols_CSR(:)	=	[1,	2,	4,	1,	2,	3,	4,	5,	1,	3,	4,	2,	5]
        A_rows_start_CSR(:)	=	[1,	4,	6,	9,	12]	 	 	 	 	 	 	 	 
        A_rows_end_CSR(:) =	[4,	6,	9,	12,	14]
        
        B_values_CSR(:)	=	[1,	-1,	-3,	-2,	5,	4,	6,	4,	-4,	2,	7,	8,	-5]
        B_cols_CSR(:)	=	[1,	2,	4,	1,	2,	3,	4,	5,	1,	3,	4,	2,	5]
        B_rows_start_CSR(:)	=	[1,	4,	6,	9,	12]	 	 	 	 	 	 	 	 
        B_rows_end_CSR(:) =	[4,	6,	9,	12,	14]
        
        N = 5
        K = 5
        stat = mkl_sparse_d_create_csr (A_sparse_CSR, SPARSE_INDEX_BASE_ONE, N, K, A_rows_start_CSR, A_rows_end_CSR, A_cols_CSR, A_values_CSR)
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        stat = mkl_sparse_d_create_csr (B_sparse_CSR, SPARSE_INDEX_BASE_ONE, N, K, B_rows_start_CSR, B_rows_end_CSR, B_cols_CSR, B_values_CSR)
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        allocate(C_dense_CSR(N,N))
        stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, A_sparse_CSR, B_sparse_CSR, SPARSE_LAYOUT_COLUMN_MAJOR, C_dense_CSR, N)
        
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        open(21,file='C_dense_CSR.txt',status='unknown',form='formatted',action='write')
        do i=1,size(C_dense_CSR,1)
            do j=1,size(C_dense_CSR,2)
                write(21,*)  C_dense_CSR(i,j)
            enddo
        enddo
        close(21)
        
        
        call displayGUIMessage( 'Starting sparse product' )
        
        stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, A_sparse_CSR, B_sparse_CSR, C_sparse_CSR)
        
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        call displayGUIMessage( 'Sparse product done' )
        
        
        
        
        
        
        
        
     
        eps_criteria = 1.0e-12
        
        NX = GridInfo%fNormX
        NY = GridInfo%fNormY
        NZ = GridInfo%fNormZ
        Areas = GridInfo%AreaFaces
        Volumes = GridInfo%Volumes
        Xel = GridInfo%Xel
        Yel = GridInfo%Yel
        Zel = GridInfo%Zel
        Xf = GridInfo%Xf
        Yf = GridInfo%Yf
        Zf = GridInfo%Zf
        T = GridInfo%TheTs
        D = GridInfo%TheDs
        Signs = GridInfo%TheSigns

        call displayGUIMessage( 'Test 1' )
        
        ! Initialize Aexch if not provided
        if (.not. present(Aexch)) then
            allocate(Aexch_local(size(Xel)))
            Aexch_local = 1.0
        else
            Aexch_local = Aexch
        end if

        ! Error checking
        if (size(Aexch_local) /= size(Xel)) then
            error stop 'Aexch must have length ns, where ns is the number of tiles'
        end if

        if (method == "GGNeumann" .and. any(Aexch_local /= 1.0)) then
            error stop 'Green-Gauss method not available for heterogeneous exchange stiffness'
        end if

        call displayGUIMessage( 'Test 2' )
        write (prog_str,'(I10)') (size(Zel))
        call displayGUIMessage( prog_str )
            
        ! Determine dimensionality
        if (any(abs(Zel-Zel(1)) .gt. 1e-12)) then
            dims = 3
        elseif (any(abs(Yel-Yel(1)) .gt. 1e-12)) then
            dims = 2
        else
            dims = 1
        end if

        call displayGUIMessage( 'Test 2-1' )
        
        ! Dimensions
        N = maxval(Signs(:,1))  ! Number of tiles
        K = maxval(Signs(:,2))  ! Number of faces

        call displayGUIMessage( 'Test 2-2' )
        
        ! DDX, DDY, DDZ
        VolCoeff = 1.0 / Volumes
        AX = NX * Areas
        AY = NY * Areas
        AZ = NZ * Areas
        ss = Signs(:,3)
        ns = Signs(:,1)
        ks = Signs(:,2)

        call displayGUIMessage( 'Test 3' )

        allocate(mask1D(size(Signs(:,1))))
                    
        ! Constructing summing matrix according to reference
        ! Constructing N times K sparse matrix DDXA: DDXA*dphi(faces) = d2phi(elements)
        ! This takes the exchange stiffness Aexch into account to form the second part of the operator div(A grad(phi))
        if (method == "DirectLaplacianNeumann") then
            ! Setting up exchange interaction strength matrix for heterogeneous materials. Also works for homogeneous materials. [2]
            
            allocate(Amat(size(ns)))

            do kk = 1, size(ks)
                
                mask1D = (ks(kk) .eq. Signs(:,2))
                k_mask1D = count(mask1D)
                
                if (kk == 1) then
                    write (prog_str,'(I10)') (k_mask1D)
                    call displayGUIMessage( prog_str )
                end if
                
                if (k_mask1D == 1) then
                    Amat(kk) = Aexch_local(ns(kk))
                    
                    if (kk == 1) then
                        write (prog_str,'(f6.3)') (Aexch_local(ns(kk)))
                        call displayGUIMessage( prog_str )
                        write (prog_str,'(f6.3)') (Amat(kk))
                        call displayGUIMessage( prog_str )
                    end if
                    
                elseif (k_mask1D > 2) then
                    write(*,*) 'Warning: more than two cells share a face!'
                    write(*,*) 'Cells number ', tmp2, ' share face number ', ks(kk)
                else
                    ns_packed = pack(ns,mask1D)
                    if (kk == 1) then
                        open(21,file='ns_packed.txt',status='unknown',form='formatted',action='write')
                        do i=1,size(ns_packed)
                            write(21,*)  ns_packed(i)
                        enddo
                        close(21)
                    end if
                    
                    Amat(kk) = 2.0 * product(Aexch_local(ns_packed)) / sum(Aexch_local(ns_packed))
                    !if all(Aexch_local(mask1D) == 0.0_wp .and. Aexch_local(tmp2) == 0.0_wp) then
                    !    Amat(kk) = 0.0_wp
                    !end if
                end if
            end do

            ! Having constructed the exchange values for each face/tile pair, we build the summing matrix.
            !real(wp), dimension(:) :: ddxA, ddyA, ddzA
            allocate(ddxA(size(ns)))
            ddxA = ss * AX(ks) * Amat * VolCoeff(ns)
            
                    ! Sort the indices of the faces and tiles
                    allocate(mask1D_3(size(ns)))    
                    mask1D_3(:) = .true.
                    allocate(ks_sorted_1(size(ns)))
                    allocate(ns_sorted_1(size(ns)))
                    allocate(ddxA_sorted(size(ns)))
                    do i = 1, size(ks)
                        indx = minloc(ks, 1, mask1D_3)
                        ks_sorted_1(i) = ks(indx)
                        ns_sorted_1(i) = ns(indx)
                        ddxA_sorted(i) = ddxA(indx)
                        mask1D_3(MINLOC(ks,mask1D_3)) = .false.
                    end do    
                    deallocate(mask1D_3) 
        
                    !Then find the non-zero values of ddxA
                    allocate(mask1D_3(size(ddxA_sorted)))    
                    mask1D_3(:) = .false.
                    mask1D_3 = (abs(ddxA_sorted) .gt. eps_criteria)
                    
                    !Then reduce the size of the arrays
                    ks_sorted_reduced = pack(ks_sorted_1,mask1D_3)
                    ns_sorted_reduced = pack(ns_sorted_1,mask1D_3)
                    ddxA_sorted_reduced = pack(ddxA_sorted,mask1D_3)
                    deallocate(mask1D_3)
                    
                    ! Determine the two column arrays needed for the CSC sparse format
                    !allocate(ns_CSC_start(size(ns_sorted_reduced)),ns_CSC_end(size(ns_sorted_reduced)))
                    !ns_CSC_start(:) = 0
                    !ns_CSC_end(:) = 0
                    !ns_CSC_start(1) = 1
                    !k_i = 2
                    !do i=1,(size(ns_sorted_reduced)-1)
                    !    if (ns_sorted_reduced(i) .ne. ns_sorted_reduced(i+1)) then
                    !        ns_CSC_start(k_i) = i+1
                    !        ns_CSC_end(k_i-1) = i+1
                    !        k_i = k_i + 1
                    !    end if
                    !enddo
                    !ns_CSC_end(k_i-1) = size(ns_sorted_reduced)
        
                    !allocate(ns_CSC_start_temp(k_i-1),ns_CSC_end_temp(k_i-1))
                    !ns_CSC_start_temp(:) = ns_CSC_start(1:(k_i-1))
                    !ns_CSC_end_temp(:) = ns_CSC_end(1:(k_i-1))
                    !call move_alloc (ns_CSC_start_temp,ns_CSC_start)
                    !call move_alloc (ns_CSC_end_temp,ns_CSC_end)
                    
                    !stat = mkl_sparse_d_create_csc (DDXA_matrix, SPARSE_INDEX_BASE_ONE, N, K, ns_CSC_start, ns_CSC_end, ks_sorted_reduced, ddxA_sorted_reduced)
                    
                    
                    ! Determine the two column arrays needed for the CSC sparse format
                    allocate(ks_CSC_start(size(ks_sorted_reduced)),ks_CSC_end(size(ks_sorted_reduced)))
                    ks_CSC_start(:) = 0
                    ks_CSC_end(:) = 0
                    ks_CSC_start(1) = 1
                    k_i = 2
                    do i=1,(size(ks_sorted_reduced)-1)
                        if (ks_sorted_reduced(i) .ne. ks_sorted_reduced(i+1)) then
                            ks_CSC_start(k_i) = i+1
                            ks_CSC_end(k_i-1) = i+1
                            k_i = k_i + 1
                        end if
                    enddo
                    ks_CSC_end(k_i-1) = size(ks_sorted_reduced)+1
        
                    allocate(ks_CSC_start_temp(k_i-1),ks_CSC_end_temp(k_i-1))
                    ks_CSC_start_temp(:) = ks_CSC_start(1:(k_i-1))
                    ks_CSC_end_temp(:) = ks_CSC_end(1:(k_i-1))
                    call move_alloc (ks_CSC_start_temp,ks_CSC_start)
                    call move_alloc (ks_CSC_end_temp,ks_CSC_end)
                    
                    
        
                    call displayGUIMessage( 'Test 3-1' )
            
            !Sort the ns array here
            !Then sort the ks and ddxA array accordingly
            !Then find the non-zero values of ddxA
            !Then reduce the size of the arrays
            !Then find the CSC start and end indices
            
            !allocate(mask1D_2(size(ddxA)))
            !mask1D_2 = (abs(ddxA) .gt. eps_criteria)
            !ddxA_non_zeroes = pack(ddxA,mask1D_2)
            !deallocate(mask1D_2)
            
            stat = mkl_sparse_d_create_csc (DDXA_matrix, SPARSE_INDEX_BASE_ONE, N, K, ks_CSC_start, ks_CSC_end, ns_sorted_reduced, ddxA_sorted_reduced)
        
            
            ! Pseudocode plan:
            ! 1. mkl_sparse_d_export_csr expects pointer arrays for row_start, row_end, col_ind, values.
            ! 2. These arrays must be declared as pointers and not allocatable.
            ! 3. The variables must be of the correct type and rank (integer, pointer, 1D for indices; real(dp), pointer, 1D for values).
            ! 4. Use => NULL() to initialize pointers, then pass them to the export routine.
            ! 5. After export, the pointers will point to the internal MKL arrays.

            ! Replace the following declarations before the call:
            ! integer, allocatable :: row_ind(:), col_ind(:), rows_start(:), rows_end(:)
            ! real(dp), allocatable :: values(:)

            

            ! Then call the export routine as follows:
            !allocate(cols_start(size(ns_CSC_start)))
            !stat = mkl_sparse_d_export_csc(DDXA_matrix, SPARSE_INDEX_BASE_ONE, N, K, cols_start, cols_end, row_ind, values)

            !open(99, file='DDXA_matrix_coo.txt', status='unknown', action='write')
            !do i = 1, nnz
            !    write(99,*) row_ind(i), values(i)
            !end do
            !close(99)
            !call displayGUIMessage( 'STAT' )
            !write (prog_str,'(I10)') (stat)
            !call displayGUIMessage( prog_str )
            
            
            !DDXA = sparse(ns, ks, ddxA, N, K)
            if (dims > 1) then
                allocate(ddyA(size(ns)))
                ddyA = ss * AY(ks) * Amat * VolCoeff(ns)
            !    DDYA = sparse(ns, ks, ddyA, N, K)
                if (dims > 2) then
                    allocate(ddzA(size(ns)))
                    ddzA = ss * AZ(ks) * Amat * VolCoeff(ns)
            !        DDZA = sparse(ns, ks, ddzA, N, K)
                end if
            end if
            
            

        else if (method == "GGNeumann") then
            ! Constructing N times K sparse matrix DDX: DDX*phi(faces) = dphi(elements)
        !    real(wp), dimension(:) :: ddx, ddy, ddz
            allocate(ddx(size(ns)))
            ddx = ss * AX(ks) * VolCoeff(ns)
        !    DDX = sparse(ns, ks, ddx, N, K)
            allocate(ddy(size(ns)))
            ddy = ss * AY(ks) * VolCoeff(ns)
        !    DDY = sparse(ns, ks, ddy, N, K)
            if (dims > 2) then
                allocate(ddz(size(ns)))
                ddz = ss * AZ(ks) * VolCoeff(ns)
        !        DDZ = sparse(ns, ks, ddz, N, K)
            end if
        end if

        deallocate(mask1D)
        call displayGUIMessage( 'Test 4' )
        
        open(21,file='ks_sorted_reduced_ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ks_sorted_reduced)
            write(21,*)  ks_sorted_reduced(i)
        enddo
        close(21)
        
        open(21,file='ks_sorted_1_ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ks_sorted_1)
            write(21,*)  ks_sorted_1(i)
        enddo
        close(21)
        
        open(21,file='ns_sorted_reduced_ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ns_sorted_reduced)
            write(21,*)  ns_sorted_reduced(i)
        enddo
        close(21)
        
        open(21,file='ns_sorted_1_ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ns_sorted_1)
            write(21,*)  ns_sorted_1(i)
        enddo
        close(21)
        
        open(21,file='ddxA_sorted_reduced_ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ddxA_sorted_reduced)
            write(21,*)  ddxA_sorted_reduced(i)
        enddo
        close(21)
        
        open(21,file='ddxA_sorted_ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ddxA_sorted)
            write(21,*)  ddxA_sorted(i)
        enddo
        close(21)
        
        open(21,file='ks_1_ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ks)
            write(21,*)  ks(i)
        enddo
        close(21)
        
        open(21,file='ns_1_ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ks)
            write(21,*)  ns(i)
        enddo
        close(21)
        
        open(21,file='ks_CSC_start_ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ks_CSC_start)
            write(21,*)  ks_CSC_start(i)
        enddo
        close(21)
        
        open(21,file='ks_CSC_end_ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ks_CSC_end)
            write(21,*)  ks_CSC_end(i)
        enddo
        close(21)
        
        open(21,file='ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ddxA)
            write(21,*)  ddxA(i)
        enddo
        close(21) 
        
        open(21,file='ddxA_non_zeroes.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ddxA_non_zeroes)
            write(21,*)  ddxA_non_zeroes(i)
        enddo
        close(21) 
        
        
        
        
        
        
        ! Defaults
        if (.not. present(interpn)) interpn = 'extended'
        if (.not. present(weight)) then
            weight = 8.0
        !else if (.not. is_numeric(weight)) then
        !    weight = real(weight, wp)
        !    if (.not. allocated(weight)) error stop 'Supplied weight is not a number.'
        end if
        !if (.not. present(method)) method = 'DirectLaplacianNeumann'

        ! Interpolation schemes. Determines how many and which neighbours to use for interpolation.
        select case (trim(interpn))
        case ('extended')
            call displayGUIMessage( 'Test extended' )
            allocate(el2fa(size(D,1),2))
            el2fa(:,1) = D(:,2)
            el2fa(:,2) = D(:,1)
        case ('compact')
            call displayGUIMessage( 'Test compact' )
            allocate(el2fa(size(Signs,1),3))
            el2fa = Signs
            write(*,*) 'Warning: untested method "compact"'
        case default
            call displayGUIMessage( 'Test default' )
            write(*,*) 'Warning: unrecognized interpolation scheme "', interpn, '". Using extended scheme.'
            allocate(el2fa(size(D,1),2))
            el2fa(:,1) = D(:,2)
            el2fa(:,2) = D(:,1)
        end select

        !! Calculating weights
        deallocate(ns,ks)
        allocate(ns(size(el2fa,1)),ks(size(el2fa,1)))
        ns = el2fa(:,1)
        ks = el2fa(:,2)
       
        
        ! Sort the indices of the faces and tiles
        allocate(mask1D(size(ks)))    
        mask1D(:) = .true.
        allocate(ks_sorted(size(ks)))
        allocate(ns_sorted(size(ks)))
        do i = 1, size(ks)
            indx = minloc(ks, 1, mask1D)
            ks_sorted(i) = ks(indx)
            ns_sorted(i) = ns(indx)
            !ks_sorted(i) = MINVAL(ks,mask1D)
            mask1D(MINLOC(ks,mask1D)) = .false.
        end do    
        deallocate(mask1D) 
        
        
        
       ! Calculating weights
       ! Determines which weights are to be used in the first interpolation step 
        allocate(w(size(ns)))
        if (dims == 1) then
            w = ((Xel(ns_sorted) - Xf(ks_sorted))**2)**(-weight/2)
        else if (dims == 2) then
            w = ((Xel(ns_sorted) - Xf(ks_sorted))**2 + (Yel(ns_sorted) - Yf(ks_sorted))**2)**(-weight/2)
        else
            w = ((Xel(ns_sorted) - Xf(ks_sorted))**2 + (Yel(ns_sorted) - Yf(ks_sorted))**2 + (Zel(ns_sorted) - Zf(ks_sorted))**2)**(-weight/2)!
        end if

        ! Prepare distances for the interpolation.       
        allocate(inds2_vals(size(ns)))
        call unique_sort(ks_sorted, inds2_vals)
        
        call displayGUIMessage( 'Before inds2' )
        
        allocate(inds2(size(inds2_vals)))
        do i = 1, size(inds2_vals)
            inds2(i) = minloc(ks_sorted, 1, mask=ks_sorted .eq. inds2_vals(i), back=.true.)
        end do
        
        call displayGUIMessage( 'After inds2' )
        
        allocate(inds1(size(inds2)+1))
        inds1(1) = 1
        inds1(2:size(inds1)) = inds2(:) + 1
        dx = Xel(ns_sorted) - Xf(ks_sorted)
        dy = Yel(ns_sorted) - Yf(ks_sorted)
        dz = Zel(ns_sorted) - Zf(ks_sorted)
        allocate(vw(size(w)))
        allocate(vx(size(w)))
        allocate(vy(size(w)))
        allocate(vz(size(w)))
        vw(:) = 0
        vx(:) = 0
        vy(:) = 0
        vz(:) = 0
        
        
        ! Scale weights to avoid ill conditioning of the least squares interpolation.
        allocate(mask1D(size(w)))
        infinity = HUGE(w) 
        mask1D = w < infinity
        wm = maxval(w,mask1D)
        w = w * (wm**(1.0/weight)) / wm
        deallocate(mask1D)
        
        ! Main loop
        ! The creation of the matrix W for the first step:
        ! W*phi(elements) = phi(faces) [Green Gauss]
        ! W*phi(elements) = dphi(faces) [Direct Laplacian]
        ! Calculated by solving
        ! Gk * [phi(faces);dphi(faces)]_k = Hk ( * phi(elements) )
        ! for each face, ks. Details can be found in [2].
        allocate(mask1D(size(Signs(:,1))))
        
        call displayGUIMessage( 'Starting ind' )
        write (prog_str,'(I20)') (K)
        call displayGUIMessage( prog_str )
               
        sjask = 0
        counter = 0
        allocate(k_log(K))
        k_log(:) = 0
        do kk = 1, K
            
            allocate(ind(inds2(kk)-inds1(kk)+1))
            k_i = 1
            
            !write (prog_str,'(I20)') (kk)
            !call displayGUIMessage( prog_str )
          
            do i = inds1(kk), inds2(kk)
                ind(k_i) = i
                k_i = k_i + 1
            enddo
          
            Wk = w(ind)
            dxk = dx(ind)
            dyk = dy(ind)
            dzk = dz(ind)
            scale_local = sum(abs(dxk))/size(dxk)
            scale = 10.0 ** nint(log10(scale_local))
          
            if (dims > 1) then
                dks = [dxk, dyk]
                scale_local = sum(abs(dks))/size(dks)
                scale = 10.0 ** nint(log10(scale_local))
                if (dims > 2) then
                    dks = [dks, dzk]
                    scale_local = sum(abs(dks))/size(dks)
                    scale = 10.0 ** nint(log10(scale_local))
                    dzk = dzk / scale
                end if
                dyk = dyk / scale
            end if
            dxk = dxk / scale

            mask1D = (kk .eq. Signs(:,2))
           
            counter = counter + 1
            nns = [NX(kk), NY(kk), NZ(kk)]
            
            if (kk == 1) then
                call displayGUIMessage( 'Sum ind 1' )
                
                open(21,file='ind.txt',status='unknown',form='formatted',action='write')
                do i=1,size(ind)
                    write(21,*)  ind(i)
                enddo
                close(21)
                
                open(21,file='dxk_test.txt',status='unknown',form='formatted',action='write')
                do i=1,size(dxk)
                    write(21,*)  dxk(i)
                enddo
                close(21)
        
                write (prog_str,'(I20)') (sum(abs(pack(Signs(:,3),mask1D))))
                call displayGUIMessage( prog_str )
                
                allocate(sjaskarr(count(mask1D)))
                sjaskarr = abs(pack(Signs(:,2),mask1D))
                
                open(21,file='sjask_1.txt',status='unknown',form='formatted',action='write')
                do i=1,size(sjaskarr)
                    write(21,*)  sjaskarr(i)
                enddo
                close(21)
                
                deallocate(sjaskarr)
                
                open(21,file='mask1D_1.txt',status='unknown',form='formatted',action='write')
                do i=1,size(mask1D)
                    write(21,*)  mask1D(i)
                enddo
                close(21)
            end if
            
            if (kk == K) then
                call displayGUIMessage( 'Sum ind end' )
        
                write (prog_str,'(I20)') (sum(abs(pack(Signs(:,3),mask1D))))
                call displayGUIMessage( prog_str )
                
                open(21,file='nns_end.txt',status='unknown',form='formatted',action='write')
                do i=1,size(nns)
                    write(21,*)  nns(i)
                enddo
                close(21)
                
                open(21,file='mask1D_end.txt',status='unknown',form='formatted',action='write')
                do i=1,size(mask1D)
                    write(21,*)  mask1D(i)
                enddo
                close(21)
                
                allocate(sjaskarr(count(mask1D)))
                sjaskarr = abs(pack(Signs(:,2),mask1D))
                
                open(21,file='sjask_end.txt',status='unknown',form='formatted',action='write')
                do i=1,size(sjaskarr)
                    write(21,*)  sjaskarr(i)
                enddo
                close(21)
                
                deallocate(sjaskarr)
                
            end if
            
            lind = size(ind)
            
            k_log(kk) = sum(abs(pack(Signs(:,3),mask1D)))
            ! Mirror trick to enforce Neumann b.c. Creates a set of virtual nodes on the other side of an edge face.
            if (sum(abs(pack(Signs(:,3),mask1D))) == 1) then
                !call displayGUIMessage( 'Test loop 1' )
                k_log(kk) = 2
                
                if ((dims .eq. 1 .and. abs(nns(1)) < eps_criteria) .or. (dims .eq. 2 .and. abs(nns(1)) < eps_criteria .and. abs(nns(2)) < eps_criteria)) then
                    k_log(kk) = 3
                    deallocate(ind)
                    cycle
                endif
                
                allocate(e(2*lind))
                allocate(extra(size(dxk,1),3))
                extra(:,1) = dxk - 2 * nns(1) * (dxk * nns(1))
                extra(:,2) = dyk - 2 * nns(2) * (dyk * nns(2))
                extra(:,3) = dzk - 2 * nns(3) * (dzk * nns(3))
                dxk2 = [dxk, extra(:,1)]
                dyk2 = [dyk, extra(:,2)]
                dzk2 = [dzk, extra(:,3)]
                Wk2 = [Wk, Wk]
            
                deallocate(extra)
            else
                allocate(e(lind))
                
                dxk2 = dxk
                dyk2 = dyk
                dzk2 = dzk
                Wk2 = Wk
            end if
            
            e(:) = 1
            
            allocate(Gkl1(size(e,1),4))
            Gkl1(:,1) = e
            Gkl1(:,2) = dxk2
            Gkl1(:,3) = dyk2
            Gkl1(:,4) = dzk2
            allocate(Gk(4,4))
            Gk(1,:) = matmul(Wk2, Gkl1)
            allocate(Gkl1_temp(size(dxk2), size(Gkl1,2)))
            do i = 1, size(Gkl1, 2)
                Gkl1_temp(:,i) = dxk2(:) * Gkl1(:,i)
            end do
            Gk(2,:) = matmul(Wk2,Gkl1_temp)
            do i = 1, size(Gkl1, 2)
                Gkl1_temp(:,i) = dyk2(:) * Gkl1(:,i)
            end do
            Gk(3,:) = matmul(Wk2,Gkl1_temp)
            do i = 1, size(Gkl1, 2)
                Gkl1_temp(:,i) = dzk2(:) * Gkl1(:,i)
            end do
            Gk(4,:) = matmul(Wk2,Gkl1_temp)
                
            Gkl1_T = transpose(Gkl1)
            allocate(Hk(4,size(e,1)))
            do i = 1, 4
                Hk(i,:) = Wk2(:)*Gkl1_T(i,:)
            enddo
                
            allocate(mask_int2log(size(nns)+1))
            mask_int2log(1) = .true.
            do i = 2, size(mask_int2log)
                if (abs(nns(i-1)) > eps_criteria) then
                    mask_int2log(i) = .true.
                else
                    mask_int2log(i) = .false.
                endif
            enddo
                
            allocate(GkRed(count(mask_int2log),count(mask_int2log)))
            allocate(HkRed(count(mask_int2log),size(Hk,2)))
            k_i = 1
            do i = 1, size(mask_int2log)
                k_row = 0
                k_j = 1
                do j = 1, size(mask_int2log)
                    if (mask_int2log(i) .and. mask_int2log(j)) then
                            GkRed(k_i,k_j) = Gk(i,j)
                            k_j = k_j + 1
                            k_row = 1
                    end if                       
                enddo
                if (k_row > 0) then
                    k_i = k_i + 1
                endif 
            enddo
                
            k_i = 1
            do i = 1, size(mask_int2log)
                if (mask_int2log(i)) then
                    HkRed(k_i,:) = Hk(i, :)
                    k_i = k_i + 1
                endif
            enddo
         
                                   
            allocate(Wktmp(size(HkRed,1),size(HkRed,2)))
            Wktmp(:,:) = HkRed(:,:)
                
            !Taken from https://www.intel.com/content/www/us/en/docs/onemkl/code-samples-lapack/2025-0/dgesv-example-fortran.html
            !call displayGUIMessage( 'Solving linear system' ) 
            if (kk == 1) then
                allocate(IPIV(size(GkRed,1)))
            else
                deallocate(IPIV)
                allocate(IPIV(size(GkRed,1)))
            endif
            
            call dgesv( size(GkRed,1), size(Wktmp,2), GkRed, size(GkRed,1), IPIV, Wktmp, size(Wktmp,1), INFO )
            !call displayGUIMessage( 'System solved' ) 
         
            if (sum(abs(pack(Signs(:,3),mask1D))) == 1) then
                vw(ind) = Wktmp(1,1:lind)+Wktmp(1,lind+1:size(Wktmp,2)) ! Interpolated face values
                if (abs(nns(1)) > eps_criteria) then ! Interpolated x-components of face gradients
                    vx(ind)=(Wktmp(2,1:lind)+Wktmp(2,lind+1:size(Wktmp,2)))/scale ! rescaled
                    if (abs(nns(2)) > eps_criteria) then ! Interpolated y-components of face gradients
                        vy(ind) = (Wktmp(3,1:lind)+Wktmp(3,lind+1:size(Wktmp,2)))/scale ! rescaled
                    endif
                elseif (abs(nns(2)) > eps_criteria) then ! Interpolated y-components of face gradients
                    vy(ind)=(Wktmp(2,1:lind)+Wktmp(2,lind+1:size(Wktmp,2)))/scale ! rescaled
                endif
                if (abs(nns(3)) > eps_criteria) then ! Interpolated z-components of face gradients
                    vz(ind)=(Wktmp(size(Wktmp,1),1:lind)+Wktmp(size(Wktmp,1),lind+1:size(Wktmp,2)))/scale ! rescaled
                endif
            else
                vw(ind) = Wktmp(1,:); ! Interpolated face values
                if (abs(nns(1)) > eps_criteria) then ! Interpolated x-components of face gradients
                    vx(ind)=Wktmp(2,:)/scale ! rescaled
                    if (abs(nns(2)) > eps_criteria) then ! Interpolated y-components of face gradients
                        vy(ind) = Wktmp(3,:)/scale ! rescaled
                    endif
                elseif (abs(nns(2)) > eps_criteria) then ! Interpolated y-components of face gradients
                    vy(ind)=Wktmp(2,:)/scale ! rescaled
                endif
                if (abs(nns(3)) > eps_criteria) then ! Interpolated z-components of face gradients
                    vz(ind)=Wktmp(size(Wktmp,1),:)/scale ! rescaled
                endif
            endif
            
        
                 
                
                if (sjask .eq. 0) then
                    call displayGUIMessage( 'Saving files 0' )
                    
                    open(21,file='HkRed.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(HkRed,1)
                        do j=1,size(HkRed,2)
                            write(21,*)  HkRed(i,j)
                        enddo
                    enddo
                    close(21)
                    
                    open(21,file='Wktmp.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(Wktmp,1)
                        do j=1,size(Wktmp,2)
                            write(21,*)  Wktmp(i,j)
                        enddo
                    enddo
                    close(21)
                                        
                    open(21,file='vw1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(vw)
                        write(21,*)  vw(i)
                    enddo
                    close(21)
                                        
                    open(21,file='vx.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(vx)
                        write(21,*)  vx(i)
                    enddo
                    close(21)
                    
                    open(21,file='vy.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(vy)
                        write(21,*)  vy(i)
                    enddo
                    close(21)
                    
                    open(21,file='vz.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(vz)
                        write(21,*)  vz(i)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 6' )
                    
                    close(21)
                    open(21,file='Wk2.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(Wk2)
                        write(21,*)  Wk2(i)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 7' )
                    
                    open(21,file='dxk.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dxk)
                        write(21,*)  dxk(i)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 8' )
                    
                    open(21,file='dyk.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dyk)
                        write(21,*)  dyk(i)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 9' )
                    
                    open(21,file='dzk.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dzk)
                        write(21,*)  dzk(i)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 10' )
                    
                    open(21,file='dxk2.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dxk2)
                        write(21,*)  dxk2(i)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 11' )
                    
                    open(21,file='dyk2.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dyk2)
                        write(21,*)  dyk2(i)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 12' )
                    
                    open(21,file='dzk2.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dzk2)
                        write(21,*)  dzk2(i)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 13' )
                    
                    open(21,file='nns.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(nns)
                        write(21,*)  nns(i)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 14' )
                    
                    open(21,file='mask_int2log.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(mask_int2log)
                        write(21,*)  mask_int2log(i)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 15' )
                    
                    open(21,file='Gkl1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(Gkl1,1)
                        write(21,*)  Gkl1(i,1)
                        write(21,*)  Gkl1(i,2)
                        write(21,*)  Gkl1(i,3)
                        write(21,*)  Gkl1(i,4)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 16' )
                    
                    open(21,file='Gk.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(Gk,1)
                        write(21,*)  Gk(i,1)
                        write(21,*)  Gk(i,2)
                        write(21,*)  Gk(i,3)
                        write(21,*)  Gk(i,4)
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 17' )
                    
                    open(21,file='Hk.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(Hk,1)
                        do j=1,size(Hk,2)
                            write(21,*)  Hk(i,j)
                        enddo
                    enddo
                    close(21)
                
                    call displayGUIMessage( 'Saving 18' )
                    
                    !open(21,file='extra.txt',status='unknown',form='formatted',action='write')
                    !do i=1,size(extra,1)
                    !    write(21,*)  extra(i,1)
                    !    write(21,*)  extra(i,2)
                    !    write(21,*)  extra(i,3)
                    !enddo
                    !close(21)
                    
                    call displayGUIMessage( 'Saving 19' )
                    
                    open(21,file='GkRed.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(GkRed,1)
                        do j=1,size(GkRed,2)
                            write(21,*)  GkRed(i,j)
                        enddo
                    enddo
                    close(21)
            
                    sjask = sjask + 1
                endif
        
                            
                
                deallocate(ind,e,Gkl1,Gkl1_temp,Hk,Gk,mask_int2log,GkRed,HkRed,Wktmp)
                
        
        end do
        
        call displayGUIMessage( 'Saving final files' )
        
        open(21,file='Amat.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Amat)
            write(21,*)  Amat(i)
        enddo
        close(21)
            
        open(21,file='ddxA.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ddxA)
            write(21,*)  ddxA(i)
        enddo
        close(21)
            
        open(21,file='w.txt',status='unknown',form='formatted',action='write')
        do i=1,size(w)
            write(21,*)  w(i)
        enddo
        close(21)
        
        open(21,file='inds2.txt',status='unknown',form='formatted',action='write')
        do i=1,size(inds2)
            write(21,*)  inds2(i)
        enddo
        close(21)
        
        open(21,file='inds1.txt',status='unknown',form='formatted',action='write')
        do i=1,size(inds1)
            write(21,*)  inds1(i)
        enddo
        close(21)
        
        open(21,file='D1.txt',status='unknown',form='formatted',action='write')
        do i=1,size(D(:,1))
            write(21,*)  D(i,1)
        enddo
        close(21)
        
        open(21,file='el2fa1.txt',status='unknown',form='formatted',action='write')
        do i=1,size(el2fa(:,1))
            write(21,*)  el2fa(i,1)
        enddo
        close(21)     
        
        open(21,file='ns.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ns)
            write(21,*)  ns(i)
        enddo
        close(21)  
        
        open(21,file='ks.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ks)
            write(21,*)  ks(i)
        enddo
        close(21)
        
        open(21,file='ks_sorted.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ks_sorted)
            write(21,*)  ks_sorted(i)
        enddo
        close(21)
        
        open(21,file='ns_sorted.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ns_sorted)
            write(21,*)  ns_sorted(i)
        enddo
        close(21)  
        
        open(21,file='k_log.txt',status='unknown',form='formatted',action='write')
        do i=1,size(k_log)
            write(21,*)  k_log(i)
        enddo
        close(21)
        
        open(21,file='vw.txt',status='unknown',form='formatted',action='write')
        do i=1,size(vw)
            write(21,*)  vw(i)
        enddo
        close(21)
        
        call displayGUIMessage( 'Saving 3' )
                    
        open(21,file='vx.txt',status='unknown',form='formatted',action='write')
        do i=1,size(vx)
            write(21,*)  vx(i)
        enddo
        close(21)
        
        call displayGUIMessage( 'Saving 4' )
                    
        open(21,file='vy.txt',status='unknown',form='formatted',action='write')
        do i=1,size(vy)
            write(21,*)  vy(i)
        enddo
        close(21)
        
        call displayGUIMessage( 'Saving 5' )
                    
        open(21,file='vz.txt',status='unknown',form='formatted',action='write')
        do i=1,size(vz)
            write(21,*)  vz(i)
        enddo
        close(21)
        
        
        
        
        !stat = mkl_sparse_d_create_csc (W_matrix, SPARSE_INDEX_BASE_ONE, K, N, ks_CSC_start, ks_CSC_end, ns, vw)
        allocate(mask1D_2(size(vx)))
        mask1D_2 = (abs(vx) .gt. eps_criteria)
        vx_non_zeroes = pack(vx,mask1D_2)
        ks_sorted_non_zeroes = pack(ks_sorted,mask1D_2)
        ns_sorted_non_zeroes = pack(ns_sorted,mask1D_2)
        deallocate(mask1D_2)
            
        open(21,file='vx_non_zeroes.txt',status='unknown',form='formatted',action='write')
        do i=1,size(vx_non_zeroes)
            write(21,*)  vx_non_zeroes(i)
        enddo
        close(21)
        
        
        !Then find the non-zero values of ddxA
        allocate(mask1D_3(size(vx_non_zeroes)))    
        mask1D_3(:) = .false.
        mask1D_3 = (abs(vx_non_zeroes) .gt. eps_criteria)
                    
        !Then reduce the size of the arrays
        ks_sorted_reduced = pack(ks_sorted_non_zeroes,mask1D_3)
        ns_sorted_reduced = pack(ns_sorted_non_zeroes,mask1D_3)
        vx_non_zeroes_reduced = pack(vx_non_zeroes,mask1D_3)
        deallocate(mask1D_3)
        
        
        
        ! Sort the indices of the faces and tiles
        allocate(mask1D_3(size(ns)))    
        mask1D_3(:) = .true.
        allocate(ks_sorted_reduced_1(size(ks_sorted_reduced)))
        allocate(ns_sorted_reduced_1(size(ks_sorted_reduced)))
        allocate(vx_sorted_reduced_1(size(ks_sorted_reduced)))
        do i = 1, size(ns_sorted_reduced_1)
            indx = minloc(ns_sorted_reduced, 1, mask1D_3)
            ks_sorted_reduced_1(i) = ks_sorted_reduced(indx)
            ns_sorted_reduced_1(i) = ns_sorted_reduced(indx)
            vx_sorted_reduced_1(i) = vx_non_zeroes_reduced(indx)
            mask1D_3(MINLOC(ns_sorted_reduced,mask1D_3)) = .false.
        end do    
        deallocate(mask1D_3) 
        
        
        !Sort the ks array here
        !Then sort the ns and vx array accordingly
        !Then find the non-zero values of vx
        !Then reduce the size of the arrays
        !Then find the CSC start and end indices
        
        
        
        
        
        
        
        
        ! Determine the two column arrays needed for the CSC sparse format
        !deallocate(ns_CSC_start,ns_CSC_end)
        allocate(ns_CSC_start(size(ns_sorted_reduced_1)),ns_CSC_end(size(ns_sorted_reduced_1)))
        ns_CSC_start(:) = 0
        ns_CSC_end(:) = 0
        ns_CSC_start(1) = 1
        k_i = 2
        do i=1,(size(ns_sorted_reduced_1)-1)
            if (ns_sorted_reduced_1(i) .ne. ns_sorted_reduced_1(i+1)) then
                ns_CSC_start(k_i) = i+1
                ns_CSC_end(k_i-1) = i+1
                k_i = k_i + 1
            end if
        enddo
        ns_CSC_end(k_i-1) = size(ns_sorted_reduced_1)+1
        
        allocate(ns_CSC_start_temp(k_i-1),ns_CSC_end_temp(k_i-1))
        ns_CSC_start_temp(:) = ns_CSC_start(1:(k_i-1))
        ns_CSC_end_temp(:) = ns_CSC_end(1:(k_i-1))
        call move_alloc (ns_CSC_start_temp,ns_CSC_start)
        call move_alloc (ns_CSC_end_temp,ns_CSC_end)
        
        
        stat = mkl_sparse_d_create_csc (FX_matrix, SPARSE_INDEX_BASE_ONE, K, N, ns_CSC_start, ns_CSC_end, ks_sorted_reduced_1, vx_sorted_reduced_1)
        
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        call displayGUIMessage( 'POSSIBLE VALUES' )
        write (prog_str,'(I10)') (SPARSE_STATUS_SUCCESS)
        call displayGUIMessage( prog_str )
        
        write (prog_str,'(I10)') (SPARSE_STATUS_NOT_INITIALIZED)
        call displayGUIMessage( prog_str )
        
        write (prog_str,'(I10)') (SPARSE_STATUS_ALLOC_FAILED)
        call displayGUIMessage( prog_str )
        
        write (prog_str,'(I10)') (SPARSE_STATUS_INVALID_VALUE)
        call displayGUIMessage( prog_str )
        
        write (prog_str,'(I10)') (SPARSE_STATUS_EXECUTION_FAILED)
        call displayGUIMessage( prog_str )
        
        write (prog_str,'(I10)') (SPARSE_STATUS_INTERNAL_ERROR)
        call displayGUIMessage( prog_str )
        
        write (prog_str,'(I10)') (SPARSE_STATUS_NOT_SUPPORTED)
        call displayGUIMessage( prog_str )
        
        open(21,file='ns_CSC_start_vx.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ns_CSC_start)
            write(21,*)  ns_CSC_start(i)
        enddo
        close(21)
        
        open(21,file='ns_CSC_end_vx.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ns_CSC_end)
            write(21,*)  ns_CSC_end(i)
        enddo
        close(21)
        
        open(21,file='ks_sorted_reduced_vx.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ks_sorted_reduced_1)
            write(21,*)  ks_sorted_reduced_1(i)
        enddo
        close(21)
        
        open(21,file='ns_sorted_reduced_vx.txt',status='unknown',form='formatted',action='write')
        do i=1,size(ns_sorted_reduced_1)
            write(21,*)  ns_sorted_reduced_1(i)
        enddo
        close(21)
        
        open(21,file='vx_non_zeroes_reduced.txt',status='unknown',form='formatted',action='write')
        do i=1,size(vx_sorted_reduced_1)
            write(21,*)  vx_sorted_reduced_1(i)
        enddo
        close(21)

        call displayGUIMessage( 'Starting dense product' )
        allocate(DX_matrix_dense(N,N))
        stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDXA_matrix, FX_matrix, SPARSE_LAYOUT_COLUMN_MAJOR, DX_matrix_dense, N)
        
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        call displayGUIMessage( 'Dense product done' )
        
        open(21,file='DX_matrix_dense.txt',status='unknown',form='formatted',action='write')
        do i=1,size(DX_matrix_dense,1)
            do j=1,size(DX_matrix_dense,2)
                write(21,*)  DX_matrix_dense(i,j)
            enddo
        enddo
        close(21)
        
        !alpha = 1.0
        !beta = 0.0
        !stat = mkl_sparse_d_sp2md (SPARSE_OPERATION_NON_TRANSPOSE, SPARSE_MATRIX_TYPE_GENERAL, DDXA_matrix, SPARSE_OPERATION_NON_TRANSPOSE, SPARSE_MATRIX_TYPE_GENERAL, FX_matrix, alpha, beta, DX_matrix_dense, SPARSE_LAYOUT_COLUMN_MAJOR, N )
        
        !allocate(DX_matrix_dense(N, N))
        call displayGUIMessage( 'Starting sparse product' )
        
        stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDXA_matrix, FX_matrix, DX_matrix)

        !stat = mkl_sparse_sp2m (SPARSE_OPERATION_NON_TRANSPOSE, SPARSE_MATRIX_TYPE_GENERAL, DDXA_matrix, SPARSE_OPERATION_NON_TRANSPOSE, SPARSE_MATRIX_TYPE_GENERAL, FX_matrix, SPARSE_STAGE_FULL_MULT, DX_matrix)
        
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        call displayGUIMessage( 'Sparse product done' )

        !allocate(DX_matrix_dense(N, N))
        
        
        
        allocate(Identity_matrix_K_CSC_start(K), Identity_matrix_K_CSC_end(K), Identity_matrix_K_rows(K), Identity_matrix_K_values(K))
        do i = 1, K
            Identity_matrix_K_CSC_start(i) = i
            Identity_matrix_K_CSC_end(i) = i+1
            Identity_matrix_K_rows(i) = i
            Identity_matrix_K_values(i) = 1.0
        enddo
        
        open(21,file='Identity_matrix_K_CSC_start.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Identity_matrix_K_CSC_start)
            write(21,*)  Identity_matrix_K_CSC_start(i)
        enddo
        close(21)
        
        open(21,file='Identity_matrix_K_CSC_end.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Identity_matrix_K_CSC_end)
            write(21,*)  Identity_matrix_K_CSC_end(i)
        enddo
        close(21)
        
        open(21,file='Identity_matrix_K_rows.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Identity_matrix_K_rows)
            write(21,*)  Identity_matrix_K_rows(i)
        enddo
        close(21)
        
        open(21,file='Identity_matrix_K_values.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Identity_matrix_K_values)
            write(21,*)  Identity_matrix_K_values(i)
        enddo
        close(21)
        
        
        stat = mkl_sparse_d_create_csc (Identity_matrix_K_sparse, SPARSE_INDEX_BASE_ONE, K, K, Identity_matrix_K_CSC_start, Identity_matrix_K_CSC_end, Identity_matrix_K_rows, Identity_matrix_K_values)
        
        
        allocate(Identity_matrix_N_CSC_start(N), Identity_matrix_N_CSC_end(N), Identity_matrix_N_rows(N), Identity_matrix_N_values(N))
        do i = 1, N
            Identity_matrix_N_CSC_start(i) = i
            Identity_matrix_N_CSC_end(i) = i+1
            Identity_matrix_N_rows(i) = i
            Identity_matrix_N_values(i) = 1.0
        enddo
        
        stat = mkl_sparse_d_create_csc (Identity_matrix_N_sparse, SPARSE_INDEX_BASE_ONE, N, N, Identity_matrix_N_CSC_start, Identity_matrix_N_CSC_end, Identity_matrix_N_rows, Identity_matrix_N_values)
        
        open(21,file='Identity_matrix_N_CSC_start.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Identity_matrix_N_CSC_start)
            write(21,*)  Identity_matrix_N_CSC_start(i)
        enddo
        close(21)
        
        open(21,file='Identity_matrix_N_CSC_end.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Identity_matrix_N_CSC_end)
            write(21,*)  Identity_matrix_N_CSC_end(i)
        enddo
        close(21)
        
        open(21,file='Identity_matrix_N_rows.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Identity_matrix_N_rows)
            write(21,*)  Identity_matrix_N_rows(i)
        enddo
        close(21)
        
        open(21,file='Identity_matrix_N_values.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Identity_matrix_N_values)
            write(21,*)  Identity_matrix_N_values(i)
        enddo
        close(21)
        
        
        call displayGUIMessage( 'Starting dense product 1' )
        allocate(FX_matrix_dense(K,N))
        stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, FX_matrix, Identity_matrix_N_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, FX_matrix_dense, K)
        
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        open(21,file='FX_matrix_dense.txt',status='unknown',form='formatted',action='write')
        do i=1,size(FX_matrix_dense,1)
            do j=1,size(FX_matrix_dense,2)
                write(21,*)  FX_matrix_dense(i,j)
            enddo
        enddo
        close(21)
        
        call displayGUIMessage( 'Dense product done 1' )
        
        
        call displayGUIMessage( 'Starting dense product 2' )
        
        allocate(DDXA_matrix_dense(N,K))
        stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDXA_matrix, Identity_matrix_K_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, DDXA_matrix_dense, N)
        
        !WORKS
        !allocate(DDXA_matrix_dense(K,K))
        !stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_TRANSPOSE, DDXA_matrix, DDXA_matrix, SPARSE_LAYOUT_COLUMN_MAJOR, DDXA_matrix_dense, K)
        
        !allocate(DDXA_matrix_dense(K,K))
        !stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, Identity_matrix_K_sparse, Identity_matrix_K_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, DDXA_matrix_dense, K)
        
        
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        call displayGUIMessage( 'Dense product done 2' )
        
        open(21,file='DDXA_matrix_dense.txt',status='unknown',form='formatted',action='write')
        do i=1,size(DDXA_matrix_dense,1)
            do j=1,size(DDXA_matrix_dense,2)
                write(21,*)  DDXA_matrix_dense(i,j)
            enddo
        enddo
        close(21)
        
        
        !stat = mkl_sparse_d_create_csr ( d2dx2%A, SPARSE_INDEX_BASE_ONE, nx*ny*nz, nx*ny*nz, d2dx2%rows_start, d2dx2%rows_end, d2dx2%cols, d2dx2%values)
                
        ! Final operation, summing interpolated values according to either ...
        !if (method == "GGNeumann") then
            ! ... the Green-Gauss theorem, yielding an estimate for the gradient
        !    W = sparse(ks, ns, vw, K, N)
        !    DX = DDX * W
        !    if (dims > 1) then
        !        DY = DDY * W
        !        if (dims > 2) then
        !            DZ = DDZ * W
        !        end if
        !    end if
        !else if (method == "DirectLaplacianNeumann") then
            ! ... the divergence theorem, yielding an estimate for the Laplacian
        !    FX = sparse(ks, ns, vx, K, N)
        !    DX = DDXA * FX
        !    if (dims > 1) then
        !        FY = sparse(ks, ns, vy, K, N)
        !        DY = DDYA * FY
        !        if (dims > 2) then
        !            FZ = sparse(ks, ns, vz, K, N)
        !            DZ = DDZA * FZ
        !        end if
        !    end if
        !else
        !    error stop 'Unrecognized method "', method, '".'
        !end if

    end subroutine computeDifferentialOperatorsFromMesh_DirectLap


    
    
subroutine unique_sort(val, unique_val)
    implicit none
    integer :: i = 0, min_val, max_val
    integer, dimension(:), intent(in) :: val 
    integer, dimension(:), allocatable :: unique
    integer, dimension(:), allocatable, intent(out) :: unique_val

    allocate(unique(size(val)))
    min_val = minval(val)-1
    max_val = maxval(val)
    do while (min_val<max_val)
        i = i+1
        min_val = minval(val, mask=val>min_val)
        unique(i) = min_val
    enddo
    allocate(unique_val(i), source=unique(1:i))   !<-- Or, just use unique(1:i) 
    !call move_alloc(unique(1:i), unique_val)
end subroutine unique_sort
    
    
    
subroutine find_nonzero_indices(array, indices, count)
    implicit none
    integer, intent(out) :: count
    integer, dimension(:), intent(out),allocatable :: indices
    real(dp), dimension(:), intent(in) :: array

    integer :: i, ns

    ns = size(array)
    count = 0

    ! First, count the number of non-zero elements
    do i = 1, ns
        if (array(i) /= 0.0) then
            count = count + 1
        end if
    end do

    ! Allocate the indices array
    allocate(indices(count))

    ! Fill the indices array with the positions of non-zero elements
    count = 0
    do i = 1, ns
        if (array(i) /= 0.0) then
            count = count + 1
            indices(count) = i
        end if
    end do
end subroutine find_nonzero_indices


subroutine get_nonzeros(array, nonzeros, count)
    implicit none
    integer, intent(out) :: count
    real(dp), dimension(:), intent(out),allocatable :: nonzeros
    real(dp), dimension(:), intent(in) :: array

    integer :: i, ns

    ns = size(array)
    count = 0

    ! First, count the number of non-zero elements
    do i = 1, ns
        if (array(i) /= 0.0) then
            count = count + 1
        end if
    end do

    ! Allocate the nonzeros array
    allocate(nonzeros(count))

    ! Fill the nonzeros array with the non-zero elements
    count = 0
    do i = 1, ns
        if (array(i) /= 0.0) then
            count = count + 1
            nonzeros(count) = array(i)
        end if
    end do
end subroutine get_nonzeros

end module DifferentialOperators