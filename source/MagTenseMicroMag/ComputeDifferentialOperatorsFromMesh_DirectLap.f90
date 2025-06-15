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
        integer, dimension(:),allocatable :: ns_copy, ks_copy, CSC_start_copy, CSC_end_copy
        real(dp), dimension(:), allocatable :: Amat
        integer, dimension(:), allocatable :: tmp, tmp2, ind, IPIV
        logical, allocatable :: mask1D(:), mask_int2log(:), mask1D_2(:), mask1D_3(:)
        real(dp), dimension(:), allocatable :: ddxA, ddyA, ddzA, ddx, ddy, ddz, dx, dy, dz, vx, vy, vz, vw, dks, dxk, dyk, dzk, nns, ddxA_copy, vx_copy
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
        integer, dimension(:),allocatable :: D_rows_start_CSR, D_rows_end_CSR, D_cols_CSR
        real(dp), dimension(:),allocatable :: D_values_CSR, vx_reduced_COO
        type(sparse_matrix_t) :: D_sparse_CSR
        type(sparse_matrix_t) :: DDXA_sparse_COO, DDXA_sparse_CSR, FX_sparse_COO, FX_sparse_CSR
        type(sparse_matrix_t) :: Identity_matrix_K_sparse_CSR, Identity_matrix_N_sparse_CSR
        integer, dimension(:),allocatable :: Identity_matrix_K_rows_start_CSR, Identity_matrix_K_rows_end_CSR, Identity_matrix_K_cols, ks_reduced_COO, ns_reduced_COO
        real(dp), dimension(:), allocatable :: Identity_matrix_K_values_CSR, ddxA_reduced_COO
        integer, dimension(:),allocatable :: Identity_matrix_N_rows_start_CSR, Identity_matrix_N_rows_end_CSR, Identity_matrix_N_cols
        real(dp), dimension(:), allocatable :: Identity_matrix_N_values_CSR
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
        
        
        
        
        ! Create identity matrices for K and N
        allocate(Identity_matrix_K_CSC_start(K), Identity_matrix_K_CSC_end(K), Identity_matrix_K_rows(K), Identity_matrix_K_values(K))
        do i = 1, K
            Identity_matrix_K_CSC_start(i) = i
            Identity_matrix_K_CSC_end(i) = i+1
            Identity_matrix_K_rows(i) = i
            Identity_matrix_K_values(i) = 1.0
        enddo
        
        stat = mkl_sparse_d_create_csc (Identity_matrix_K_sparse, SPARSE_INDEX_BASE_ONE, K, K, Identity_matrix_K_CSC_start, Identity_matrix_K_CSC_end, Identity_matrix_K_rows, Identity_matrix_K_values)
        
        
        allocate(Identity_matrix_N_CSC_start(N), Identity_matrix_N_CSC_end(N), Identity_matrix_N_rows(N), Identity_matrix_N_values(N))
        do i = 1, N
            Identity_matrix_N_CSC_start(i) = i
            Identity_matrix_N_CSC_end(i) = i+1
            Identity_matrix_N_rows(i) = i
            Identity_matrix_N_values(i) = 1.0
        enddo
        
        stat = mkl_sparse_d_create_csc (Identity_matrix_N_sparse, SPARSE_INDEX_BASE_ONE, N, N, Identity_matrix_N_CSC_start, Identity_matrix_N_CSC_end, Identity_matrix_N_rows, Identity_matrix_N_values)
        
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
            
            
            open(21,file='ks_before.txt',status='unknown',form='formatted',action='write')
            do i=1,size(ks)
                write(21,*)  ks(i)
            enddo
            close(21)
            
            open(21,file='ns_before.txt',status='unknown',form='formatted',action='write')
            do i=1,size(ns)
                write(21,*)  ns(i)
            enddo
            close(21)
            
            open(21,file='ddxA_before.txt',status='unknown',form='formatted',action='write')
            do i=1,size(ddxA)
                write(21,*)  ddxA(i)
            enddo
            close(21)
        
            ! Determine the two column arrays needed for the CSC sparse format
            allocate(ns_copy(size(ns)), ks_copy(size(ns)), ddxA_copy(size(ns)))
            ns_copy(:) = ns(:)
            ks_copy(:) = ks(:)
            !ddxA_copy(:) = ddxA(:)
            allocate(CSC_start_copy(size(ks)),CSC_end_copy(size(ks)))
            call displayGUIMessage( 'Starting make_CSR_CSC_indices' )
            !call make_CSR_CSC_indices(ns,ks,ddxA,CSC_start_copy,CSC_end_copy,ns_copy,ddxA_copy)        
            call make_CSR_CSC_indices(ks,ns,ddxA,CSC_start_copy,CSC_end_copy,ks_copy,ddxA_copy)   
            call displayGUIMessage( 'Ending make_CSR_CSC_indices' )
            
            open(21,file='ns_copy_ddxA.txt',status='unknown',form='formatted',action='write')
            do i=1,size(ns_copy)
                write(21,*)  ns_copy(i)
            enddo
            close(21)
            
            open(21,file='ks_copy_ddxA.txt',status='unknown',form='formatted',action='write')
            do i=1,size(ks_copy)
                write(21,*)  ks_copy(i)
            enddo
            close(21)
            
            open(21,file='ddxA_copy_ddxA.txt',status='unknown',form='formatted',action='write')
            do i=1,size(ddxA_copy)
                write(21,*)  ddxA_copy(i)
            enddo
            close(21)
                        
            open(21,file='CSC_start_copy_ddxA.txt',status='unknown',form='formatted',action='write')
            do i=1,size(CSC_start_copy)
                write(21,*)  CSC_start_copy(i)
            enddo
            close(21)
                        
            open(21,file='CSC_end_copy_ddxA.txt',status='unknown',form='formatted',action='write')
            do i=1,size(CSC_end_copy)
                write(21,*)  CSC_end_copy(i)
            enddo
            close(21)
                                    
            stat = mkl_sparse_d_create_csr (DDXA_matrix, SPARSE_INDEX_BASE_ONE, N, K, CSC_start_copy, CSC_end_copy, ks_copy, ddxA_copy)
        
            
            
            
            
            
            allocate(D_values_CSR(8), D_cols_CSR(8), D_rows_start_CSR(4), D_rows_end_CSR(4))
        D_values_CSR(:) = [ -160000000.0, 160000000.0, -160000000.0, 160000000.0, -160000000.0, 160000000.0, -160000000.0, 160000000.0 ]
        !D_values_CSR(:) = [ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 ]
        D_rows_start_CSR(:) = [ 1, 3, 5, 7 ]
        D_rows_end_CSR(:) = [ 3, 5, 7, 9 ]
        D_cols_CSR(:) = [ 1, 2, 2, 5, 3, 4, 4, 6 ]
        N = 4
        K = 20
        stat = mkl_sparse_d_create_csr(D_sparse_CSR, SPARSE_INDEX_BASE_ONE, N, K, D_rows_start_CSR, D_rows_end_CSR, D_cols_CSR, D_values_CSR)
        call displayGUIMessage('STAT')
        write(prog_str, '(I10)') (stat)
        call displayGUIMessage(prog_str)
        
        
        ! Create identity matrices for K
        allocate(Identity_matrix_K_rows_start_CSR(K), Identity_matrix_K_rows_end_CSR(K), Identity_matrix_K_cols(K), Identity_matrix_K_values_CSR(K))
        do i = 1, K
            Identity_matrix_K_rows_start_CSR(i) = i
            Identity_matrix_K_rows_end_CSR(i) = i+1
            Identity_matrix_K_cols(i) = i
            Identity_matrix_K_values_CSR(i) = 1.0
        enddo
        
        stat = mkl_sparse_d_create_csr (Identity_matrix_K_sparse_CSR, SPARSE_INDEX_BASE_ONE, K, K, Identity_matrix_K_rows_start_CSR, Identity_matrix_K_rows_end_CSR, Identity_matrix_K_cols, Identity_matrix_K_values_CSR)
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        ! Create identity matrices for N
        allocate(Identity_matrix_N_rows_start_CSR(N), Identity_matrix_N_rows_end_CSR(N), Identity_matrix_N_cols(N), Identity_matrix_N_values_CSR(N))
        do i = 1, N
            Identity_matrix_N_rows_start_CSR(i) = i
            Identity_matrix_N_rows_end_CSR(i) = i+1
            Identity_matrix_N_cols(i) = i
            Identity_matrix_N_values_CSR(i) = 1.0
        enddo
        
        stat = mkl_sparse_d_create_csr (Identity_matrix_N_sparse_CSR, SPARSE_INDEX_BASE_ONE, N, N, Identity_matrix_N_rows_start_CSR, Identity_matrix_N_rows_end_CSR, Identity_matrix_N_cols, Identity_matrix_N_values_CSR)
        call displayGUIMessage( 'STAT' )
        write (prog_str,'(I10)') (stat)
        call displayGUIMessage( prog_str )
        
        
        deallocate(mask1D)
        
        
        
         !Then find the non-zero values of ddxA
        allocate(mask1D(size(ddxA)))    
        mask1D(:) = .false.
        mask1D = (abs(ddxA) .gt. eps_criteria)
        
        !Then reduce the size of the arrays
        ks_reduced_COO    = pack(ks,mask1D)
        ns_reduced_COO = pack(ns,mask1D)
        ddxA_reduced_COO  = pack(ddxA,mask1D)
        deallocate(mask1D)
     
        nnz = size(ddxA_reduced_COO)
        stat = mkl_sparse_d_create_coo (DDXA_sparse_COO, SPARSE_INDEX_BASE_ONE, N, K, nnz, ns_reduced_COO, ks_reduced_COO, ddxA_reduced_COO)
        
        stat = mkl_sparse_convert_csr (DDXA_sparse_COO, SPARSE_OPERATION_NON_TRANSPOSE, DDXA_sparse_CSR)
        
        
        
        
        
        
        
        
            
            call displayGUIMessage( 'Starting DDXA dense product' )
        
            allocate(DDXA_matrix_dense(N,K))
            DDXA_matrix_dense(:,:) = 0.0
            !stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, D_sparse_CSR, Identity_matrix_K_sparse_CSR, SPARSE_LAYOUT_COLUMN_MAJOR, DDXA_matrix_dense, N)
            !stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDXA_matrix, Identity_matrix_K_sparse_CSR, SPARSE_LAYOUT_COLUMN_MAJOR, DDXA_matrix_dense, N)
            stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDXA_sparse_CSR, Identity_matrix_K_sparse_CSR, SPARSE_LAYOUT_COLUMN_MAJOR, DDXA_matrix_dense, N)
            
            call displayGUIMessage( 'STAT' )
            write (prog_str,'(I10)') (stat)
            call displayGUIMessage( prog_str )
        
            call displayGUIMessage( 'Ending DDXA dense product' )
        
            open(21,file='DDXA_matrix_dense.txt',status='unknown',form='formatted',action='write')
            do i=1,size(DDXA_matrix_dense,1)
                do j=1,size(DDXA_matrix_dense,2)
                    write(21,*)  DDXA_matrix_dense(i,j)
                enddo
            enddo
            close(21)
            
            stat = mkl_sparse_destroy(D_sparse_CSR)
            stat = mkl_sparse_destroy(Identity_matrix_K_sparse_CSR)
            
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
        
            open(21,file='ns_sorted_reduced_ddxA.txt',status='unknown',form='formatted',action='write')
            do i=1,size(ns_sorted_reduced)
                write(21,*)  ns_sorted_reduced(i)
            enddo
            close(21)
               
            open(21,file='ddxA_sorted_reduced_ddxA.txt',status='unknown',form='formatted',action='write')
            do i=1,size(ddxA_sorted_reduced)
                write(21,*)  ddxA_sorted_reduced(i)
            enddo
            close(21)
        
            
            
            
            
            
            
        
            
            
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

        !deallocate(mask1D)
        call displayGUIMessage( 'Test 4' )
        
       
        
        
        
        
        
        
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
        
        
        
        
        ! Begin block comment
        !
        ! Determine the two column arrays needed for the CSC sparse format
        ! deallocate(ns_copy, ks_copy, ddxA_copy, CSC_start_copy ,CSC_end_copy)
        ! allocate(ns_copy(size(ns)), ks_copy(size(ns)), ddxA_copy(size(ns)))
        ! !ns_copy(:) = ns_sorted(:)
        ! !ks_copy(:) = ks_sorted(:)
        ! !vx_copy(:) = vx(:)
        ! allocate(CSC_start_copy(size(ks)),CSC_end_copy(size(ks)))
        ! call make_CSR_CSC_indices(ks,ns,vx,CSC_start_copy,CSC_end_copy,ks_copy,vx_copy) 
        ! 
        ! stat = mkl_sparse_d_create_csc (FX_matrix, SPARSE_INDEX_BASE_ONE, K, N, CSC_start_copy, CSC_end_copy, ks_copy, vx_copy)
        ! 
        ! call displayGUIMessage( 'STAT' )
        ! write (prog_str,'(I10)') (stat)
        ! call displayGUIMessage( prog_str )
        ! 
        ! open(21,file='CSC_start_copy_vx.txt',status='unknown',form='formatted',action='write')
        ! do i=1,size(CSC_start_copy)
        !     write(21,*)  CSC_start_copy(i)
        ! enddo
        ! close(21)
        ! 
        ! open(21,file='CSC_end_copy_vx.txt',status='unknown',form='formatted',action='write')
        ! do i=1,size(CSC_end_copy)
        !     write(21,*)  CSC_end_copy(i)
        ! enddo
        ! close(21)
        ! 
        ! open(21,file='ks_copy_vx.txt',status='unknown',form='formatted',action='write')
        ! do i=1,size(ks_copy)
        !     write(21,*)  ks_copy(i)
        ! enddo
        ! close(21)
        !        
        ! open(21,file='vx_copy_vx.txt',status='unknown',form='formatted',action='write')
        ! do i=1,size(vx_copy)
        !     write(21,*)  vx_copy(i)
        ! enddo
        ! close(21)
        !         
        
        
        
        deallocate(mask1D)
         !Then find the non-zero values of ddxA
        allocate(mask1D(size(vx)))    
        mask1D(:) = .false.
        mask1D = (abs(vx) .gt. eps_criteria)
        
        !Then reduce the size of the arrays
        ns_reduced_COO = pack(ns_sorted,mask1D)
        ks_reduced_COO = pack(ks_sorted,mask1D)
        vx_reduced_COO = pack(vx,mask1D)
        deallocate(mask1D)
        
        nnz = size(vx_reduced_COO)
        stat = mkl_sparse_d_create_coo (FX_sparse_COO, SPARSE_INDEX_BASE_ONE, K, N, nnz, ks_reduced_COO, ns_reduced_COO, vx_reduced_COO)
        
        stat = mkl_sparse_convert_csr (FX_sparse_COO, SPARSE_OPERATION_NON_TRANSPOSE, FX_sparse_CSR)
        
         call displayGUIMessage( 'Starting FX dense product' )
        
            allocate(FX_matrix_dense(K,N))
            FX_matrix_dense(:,:) = 0.0
            stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, FX_sparse_CSR, Identity_matrix_N_sparse_CSR, SPARSE_LAYOUT_COLUMN_MAJOR, FX_matrix_dense, K)
            
            call displayGUIMessage( 'STAT' )
            write (prog_str,'(I10)') (stat)
            call displayGUIMessage( prog_str )
        
            call displayGUIMessage( 'Ending FX dense product' )
        
            open(21,file='FX_matrix_dense.txt',status='unknown',form='formatted',action='write')
            do i=1,size(FX_matrix_dense,1)
                do j=1,size(FX_matrix_dense,2)
                    write(21,*)  FX_matrix_dense(i,j)
                enddo
            enddo
            close(21)
        
        
        ! 
        ! call displayGUIMessage( 'Starting FX_matrix dense product' )
        ! allocate(FX_matrix_dense(K,N))
        ! stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, FX_matrix, Identity_matrix_N_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, FX_matrix_dense, K)
        ! 
        ! call displayGUIMessage( 'STAT' )
        ! write (prog_str,'(I10)') (stat)
        ! call displayGUIMessage( prog_str )
        ! 
        ! open(21,file='FX_matrix_dense.txt',status='unknown',form='formatted',action='write')
        ! do i=1,size(FX_matrix_dense,1)
        !     do j=1,size(FX_matrix_dense,2)
        !         write(21,*)  FX_matrix_dense(i,j)
        !     enddo
        ! enddo
        ! close(21)
        ! 
        ! call displayGUIMessage( 'Ending FX_matrix dense product' )
        ! 
        ! 
        ! 
        ! 
        ! 
        ! 
        ! 
        ! 
        ! 
        ! 
         call displayGUIMessage( 'Starting dense product' )
         allocate(DX_matrix_dense(N,N))
         stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDXA_sparse_CSR, FX_sparse_CSR, SPARSE_LAYOUT_COLUMN_MAJOR, DX_matrix_dense, N)
         
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
         
        ! !alpha = 1.0
        ! !beta = 0.0
        ! !stat = mkl_sparse_d_sp2md (SPARSE_OPERATION_NON_TRANSPOSE, SPARSE_MATRIX_TYPE_GENERAL, DDXA_matrix, SPARSE_OPERATION_NON_TRANSPOSE, SPARSE_MATRIX_TYPE_GENERAL, FX_matrix, alpha, beta, DX_matrix_dense, SPARSE_LAYOUT_COLUMN_MAJOR, N )
        ! 
        ! call displayGUIMessage( 'Starting sparse product' )
        ! 
        ! stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDXA_matrix, FX_matrix, DX_matrix)
        ! 
        ! !stat = mkl_sparse_sp2m (SPARSE_OPERATION_NON_TRANSPOSE, SPARSE_MATRIX_TYPE_GENERAL, DDXA_matrix, SPARSE_OPERATION_NON_TRANSPOSE, SPARSE_MATRIX_TYPE_GENERAL, FX_matrix, SPARSE_STAGE_FULL_MULT, DX_matrix)
        ! 
        ! call displayGUIMessage( 'STAT' )
        ! write (prog_str,'(I10)') (stat)
        ! call displayGUIMessage( prog_str )
        ! 
        ! call displayGUIMessage( 'Sparse product done' )
        ! 
        ! !allocate(DX_matrix_dense(N, N))
        !
        ! End block comment
        
        
        
        
        
        
        
        
        
        
        
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
    

subroutine make_CSR_CSC_indices(rows,columns,values,CSC_start,CSC_end,rows_sorted,values_sorted)
    implicit none
    
    integer, dimension(:), intent(in) :: rows, columns
    real(dp), dimension(:), intent(in) :: values
    integer, dimension(:), allocatable, intent(inout) :: CSC_start, CSC_end
    integer, dimension(:), allocatable, intent(inout) :: rows_sorted
    real(dp), dimension(:), allocatable, intent(inout) :: values_sorted
    logical, allocatable :: mask1D(:)
    real(dp) :: eps_criteria
    integer :: i, j, indx, k_i
    integer, dimension(:), allocatable :: rows_reduced, columns_reduced, columns_sorted
    real(dp), dimension(:), allocatable :: values_reduced
    integer, dimension(:), allocatable :: rows_sorted_temp
    real(dp), dimension(:), allocatable :: values_sorted_temp
    integer, dimension(:), allocatable :: CSC_start_temp, CSC_end_temp
    integer, allocatable :: sort_idx(:)
    integer :: nvals
    character*(40) :: prog_str
    
    eps_criteria = 1.0e-12
    
    !Find the non-zero values of values
    !Then reduce the size of the arrays    
    !Then sort the colums array
    !Then sort the rows and values array accordingly
    !Then find the CSC start and end indices
    
    !Then find the non-zero values of ddxA
    allocate(mask1D(size(values)))    
    mask1D(:) = .false.
    mask1D = (abs(values) .gt. eps_criteria)
        
    !Then reduce the size of the arrays
    rows_reduced    = pack(rows,mask1D)
    columns_reduced = pack(columns,mask1D)
    values_reduced  = pack(values,mask1D)
    deallocate(mask1D)
        
    open(21,file='columns_reduced.txt',status='unknown',form='formatted',action='write')
    do i=1,size(columns_reduced)
        write(21,*)  columns_reduced(i)
    enddo
    close(21)
    
    open(21,file='rows_reduced.txt',status='unknown',form='formatted',action='write')
    do i=1,size(rows_reduced)
        write(21,*)  rows_reduced(i)
    enddo
    close(21)
    
    open(21,file='values_reduced.txt',status='unknown',form='formatted',action='write')
    do i=1,size(values_reduced)
        write(21,*)  values_reduced(i)
    enddo
    close(21)
            
    ! Sort the indices of the faces and tiles
    allocate(mask1D(size(columns_reduced)))    
    mask1D(:) = .true.
    allocate(columns_sorted(size(columns_reduced)))
    allocate(rows_sorted_temp(size(columns_reduced)))
    allocate(values_sorted_temp(size(columns_reduced)))
    !do i = 1, size(columns_reduced)
    !    indx = minloc(columns_reduced, 1, mask1D)
    !    columns_sorted(i)      = columns_reduced(indx)
    !    rows_sorted_temp(i)    = rows_reduced(indx)
    !    values_sorted_temp(i)  = values_reduced(indx)
    !    mask1D(indx) = .false.
        
    !    call displayGUIMessage( 'columns_sorted(i)' )
    !    write (prog_str,'(I10)') (columns_sorted(i))
    !    call displayGUIMessage( prog_str )
    !    write (prog_str,'(I10)') (indx)
    !    call displayGUIMessage( prog_str )
        
    !end do    
! Create an array of indices to sort by columns_reduced, then by rows_reduced for ties


    nvals = size(columns_reduced)
    allocate(sort_idx(nvals))
    sort_idx = [(i, i=1,nvals)]

    ! Custom insertion sort: sort by columns_reduced, then by rows_reduced for ties
    do i = 2, nvals
        indx = sort_idx(i)
        j = i - 1
        do while (j >= 1)
            if (columns_reduced(sort_idx(j)) > columns_reduced(indx)) then
                sort_idx(j+1) = sort_idx(j)
                j = j - 1
            else if (columns_reduced(sort_idx(j)) == columns_reduced(indx)) then
                if (rows_reduced(sort_idx(j)) > rows_reduced(indx)) then
                    sort_idx(j+1) = sort_idx(j)
                    j = j - 1
                else
                    exit
                end if
            else
                exit
            end if
        end do
        sort_idx(j+1) = indx
    end do

    ! Apply the sorted indices
    do i = 1, nvals
        columns_sorted(i)      = columns_reduced(sort_idx(i))
        rows_sorted_temp(i)    = rows_reduced(sort_idx(i))
        values_sorted_temp(i)  = values_reduced(sort_idx(i))
    end do
    deallocate(mask1D) 
    
    open(21,file='columns_sorted.txt',status='unknown',form='formatted',action='write')
    do i=1,size(columns_sorted)
        write(21,*)  columns_sorted(i)
    enddo
    close(21)
    
    open(21,file='rows_sorted_temp.txt',status='unknown',form='formatted',action='write')
    do i=1,size(rows_sorted_temp)
        write(21,*)  rows_sorted_temp(i)
    enddo
    close(21)
    
    open(21,file='values_sorted_temp.txt',status='unknown',form='formatted',action='write')
    do i=1,size(values_sorted_temp)
        write(21,*)  values_sorted_temp(i)
    enddo
    close(21)
    
        
    ! Determine the two column arrays needed for the CSC sparse format
    !allocate(CSC_start(size(columns_sorted)),CSC_end(size(columns_sorted)))
    CSC_start(:) = 0
    CSC_end(:) = 0
    CSC_start(1) = 1
    k_i = 2
    do i=1,(size(columns_sorted)-1)
        if (columns_sorted(i) .ne. columns_sorted(i+1)) then
            CSC_start(k_i) = i+1
            CSC_end(k_i-1) = i+1
            k_i = k_i + 1
        end if
    enddo
    CSC_end(k_i-1) = size(columns_sorted)+1
        
    allocate(CSC_start_temp(k_i-1),CSC_end_temp(k_i-1))
    CSC_start_temp(:) = CSC_start(1:(k_i-1))
    CSC_end_temp(:)   = CSC_end(1:(k_i-1))
    call move_alloc (CSC_start_temp,CSC_start)
    call move_alloc (CSC_end_temp,CSC_end)
        
    call move_alloc (rows_sorted_temp,rows_sorted)
    call move_alloc (values_sorted_temp,values_sorted)

end subroutine make_CSR_CSC_indices
    


end module DifferentialOperators