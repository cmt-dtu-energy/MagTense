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
        integer, intent(in) :: interpn
        integer, intent(in) :: method
        !character(len=*), intent(inout), optional :: interpn
        real(dp), intent(inout), optional :: weight
        !character(len=*), intent(inout), optional :: method
        real(dp), dimension(:), intent(in), optional :: Aexch
        real(dp), dimension(:,:), intent(in), optional :: ExtW
        type(sparse_matrix_t) :: DX_matrix, DY_matrix, DZ_matrix

        ! Unpack the variables
        real(dp), dimension(:),allocatable :: NX, NY, NZ, Areas, Volumes, Xel, Yel, Zel, Xf, Yf, Zf
        integer, dimension(:,:),allocatable :: Signs, T, D, el2fa
        real(dp), dimension(:), allocatable :: Aexch_local
        integer :: i,j,dims, N, K, kk, k_mask1D, k_i, lind, sjask, indx, k_j, k_row, info, n_unique
        real(dp), dimension(:),allocatable :: VolCoeff, AX, AY, AZ, w
        integer, dimension(:),allocatable :: ss, ns, ks, inds1, inds2, inds2_vals, inds2_vals_temp, ks_sorted, ns_sorted
        integer, dimension(:),allocatable :: ns_packed, k_log
        real(dp), dimension(:), allocatable :: Amat
        integer, dimension(:), allocatable :: ind, IPIV
        logical, allocatable :: mask1D(:), mask_int2log(:)
        real(dp), dimension(:), allocatable :: ddxA, ddyA, ddzA, ddx, ddy, ddz, dx, dy, dz, vx, vy, vz, vw, dks, dxk, dyk, dzk, nns
        real(dp), dimension(:), allocatable :: dxk2, dyk2, dzk2, Wk, Wk2, e
        real(dp), dimension(:,:), allocatable :: Gkl1, Gk, Hk, Gkl1_temp, Gkl1_T, GkRed, HkRed, Wktmp
        real(dp) :: wm, infinity, scale, scale_local, eps_criteria
        real(dp), allocatable :: extra(:,:), DX_matrix_dense(:,:), DY_matrix_dense(:,:), DZ_matrix_dense(:,:)
        integer :: stat                                   !> Status value for the various sparse matrix operations 
        character*(40) :: prog_str
        type(sparse_matrix_t) :: descr
        real(dp), dimension(:,:), allocatable :: DDXA_matrix_dense, FX_matrix_dense
        type(sparse_matrix_t) :: DDXA_sparse, FX_sparse, DDYA_sparse, FY_sparse, DDZA_sparse, FZ_sparse
        type(sparse_matrix_t) :: DDX_sparse, DDY_sparse, DDZ_sparse, W_sparse
        type(sparse_matrix_t) :: Identity_matrix_K_sparse, Identity_matrix_N_sparse
        integer, dimension(:),allocatable :: Identity_matrix_K_rows_start_CSR, Identity_matrix_K_rows_end_CSR, Identity_matrix_K_cols
        real(dp), dimension(:), allocatable :: Identity_matrix_K_values_CSR
        integer, dimension(:),allocatable :: Identity_matrix_N_rows_start_CSR, Identity_matrix_N_rows_end_CSR, Identity_matrix_N_cols
        real(dp), dimension(:), allocatable :: Identity_matrix_N_values_CSR
        
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

        ! Determine dimensionality
        if (any(abs(Zel-Zel(1)) .gt. 1e-12)) then
            dims = 3
        elseif (any(abs(Yel-Yel(1)) .gt. 1e-12)) then
            dims = 2
        else
            dims = 1
        end if
        
        ! Dimensions
        N = maxval(Signs(:,1))  ! Number of tiles
        K = maxval(Signs(:,2))  ! Number of faces
                
        ! DDX, DDY, DDZ
        VolCoeff = 1.0 / Volumes
        AX = NX * Areas
        AY = NY * Areas
        AZ = NZ * Areas
        ss = Signs(:,3)
        ns = Signs(:,1)
        ks = Signs(:,2)

        ! Create identity matrices for K
        allocate(Identity_matrix_K_rows_start_CSR(K), Identity_matrix_K_rows_end_CSR(K), Identity_matrix_K_cols(K), Identity_matrix_K_values_CSR(K))
        do i = 1, K
            Identity_matrix_K_rows_start_CSR(i) = i
            Identity_matrix_K_rows_end_CSR(i) = i+1
            Identity_matrix_K_cols(i) = i
            Identity_matrix_K_values_CSR(i) = 1.0
        enddo
        
        stat = mkl_sparse_d_create_csr (Identity_matrix_K_sparse, SPARSE_INDEX_BASE_ONE, K, K, Identity_matrix_K_rows_start_CSR, Identity_matrix_K_rows_end_CSR, Identity_matrix_K_cols, Identity_matrix_K_values_CSR)
        
        ! Create identity matrices for N
        allocate(Identity_matrix_N_rows_start_CSR(N), Identity_matrix_N_rows_end_CSR(N), Identity_matrix_N_cols(N), Identity_matrix_N_values_CSR(N))
        do i = 1, N
            Identity_matrix_N_rows_start_CSR(i) = i
            Identity_matrix_N_rows_end_CSR(i) = i+1
            Identity_matrix_N_cols(i) = i
            Identity_matrix_N_values_CSR(i) = 1.0
        enddo
        
        stat = mkl_sparse_d_create_csr (Identity_matrix_N_sparse, SPARSE_INDEX_BASE_ONE, N, N, Identity_matrix_N_rows_start_CSR, Identity_matrix_N_rows_end_CSR, Identity_matrix_N_cols, Identity_matrix_N_values_CSR)       
        
        ! Constructing summing matrix according to reference
        ! Constructing N times K sparse matrix DDXA: DDXA*dphi(faces) = d2phi(elements)
        ! This takes the exchange stiffness Aexch into account to form the second part of the operator div(A grad(phi))
        if ( method .eq. MicroMagExchMethodDirectLaplacianNeumann ) then
        !if (method == "DirectLaplacianNeumann") then
            ! Setting up exchange interaction strength matrix for heterogeneous materials. Also works for homogeneous materials. [2]
            allocate(Amat(size(ns)))
            allocate(mask1D(size(Signs(:,1))))
        
            do kk = 1, size(ks)
                
                mask1D = (ks(kk) .eq. Signs(:,2))
                k_mask1D = count(mask1D)
                                
                if (k_mask1D == 1) then
                    Amat(kk) = Aexch_local(ns(kk))
                elseif (k_mask1D > 2) then
                    write(*,*) 'Warning: more than two cells share a face!'
                else
                    ns_packed = pack(ns,mask1D)
                    
                    Amat(kk) = 2.0 * product(Aexch_local(ns_packed)) / sum(Aexch_local(ns_packed))
                    !IMPLEMENT THE CHECK BELOW
                    !if all(Aexch_local(mask1D) == 0.0_wp .and. Aexch_local(tmp2) == 0.0_wp) then
                    !    Amat(kk) = 0.0_wp
                    !end if
                end if
            end do

            deallocate(mask1D)
             
            ! Having constructed the exchange values for each face/tile pair, we build the summing matrix.
            allocate(ddxA(size(ns)))
            ddxA = ss * AX(ks) * Amat * VolCoeff(ns)
            !DDXA = sparse(ns, ks, ddxA, N, K)
            call create_CSR_matrix(ns, ks, ddxA, N, K, eps_criteria, DDXA_sparse)
           
            if (dims > 1) then
                allocate(ddyA(size(ns)))
                ddyA = ss * AY(ks) * Amat * VolCoeff(ns)
                !DDYA = sparse(ns, ks, ddyA, N, K)
                call create_CSR_matrix(ns, ks, ddYA, N, K, eps_criteria, DDYA_sparse)
                
                if (dims > 2) then
                    allocate(ddzA(size(ns)))
                    ddzA = ss * AZ(ks) * Amat * VolCoeff(ns)
                    !DDZA = sparse(ns, ks, ddzA, N, K)
                    call create_CSR_matrix(ns, ks, ddzA, N, K, eps_criteria, DDZA_sparse)
                end if
            end if
        
        else if ( method .eq. MicroMagExchMethodGGNeumann ) then
        !else if (method == "GGNeumann") then
            ! Constructing N times K sparse matrix DDX: DDX*phi(faces) = dphi(elements)
            allocate(ddx(size(ns)))
            ddx = ss * AX(ks) * VolCoeff(ns)
            !DDX = sparse(ns, ks, ddx, N, K)
            call create_CSR_matrix(ns, ks, ddx, N, K, eps_criteria, DDX_sparse)
            
            allocate(ddy(size(ns)))
            ddy = ss * AY(ks) * VolCoeff(ns)
            !DDY = sparse(ns, ks, ddy, N, K)
            call create_CSR_matrix(ns, ks, ddy, N, K, eps_criteria, DDY_sparse)
                        
            if (dims > 2) then
                allocate(ddz(size(ns)))
                ddz = ss * AZ(ks) * VolCoeff(ns)
                !DDZ = sparse(ns, ks, ddz, N, K)
                call create_CSR_matrix(ns, ks, ddz, N, K, eps_criteria, DDZ_sparse)
            end if
        end if


        call displayGUIMessage( 'DD-matrices constructed' )
        
        ! Defaults
        !if (.not. present(interpn)) interpn = 'extended'
        !if (.not. present(weight)) then
        !    weight = 8.0
        !else if (.not. is_numeric(weight)) then
        !    weight = real(weight, wp)
        !    if (.not. allocated(weight)) error stop 'Supplied weight is not a number.'
        !end if
        !if (.not. present(method)) method = 'DirectLaplacianNeumann'

        ! Interpolation schemes. Determines how many and which neighbours to use for interpolation.
        
        !select case (trim(interpn))
        if ( interpn .eq. MicroMagExchInterpnExtended ) then
        !case ('extended')
            call displayGUIMessage( 'Test extended' )
            allocate(el2fa(size(D,1),2))
            el2fa(:,1) = D(:,2)
            el2fa(:,2) = D(:,1)
        elseif ( interpn .eq. MicroMagExchInterpnCompact ) then
        !case ('compact')
            call displayGUIMessage( 'Test compact' )
            allocate(el2fa(size(Signs,1),3))
            el2fa = Signs
            call displayGUIMessage( 'Warning: untested method: compact' )
        !case default
        !    call displayGUIMessage( 'Test default' )
        !    write(*,*) 'Warning: unrecognized interpolation scheme "', interpn, '". Using extended scheme.'
        !    allocate(el2fa(size(D,1),2))
        !    el2fa(:,1) = D(:,2)
        !    el2fa(:,2) = D(:,1)
        !end select
        endif

        ! Calculating weights
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
            mask1D(MINLOC(ks,mask1D)) = .false.
        end do    
        deallocate(mask1D) 
                
       write (prog_str,'(F10.2)') (weight)
       call displayGUIMessage( prog_str )
        
       ! Calculating weights
       ! Determines which weights are to be used in the first interpolation step 
        allocate(w(size(ns)))
        if (dims == 1) then
            w = ((Xel(ns_sorted) - Xf(ks_sorted))**2)**(-weight/2.0)
        else if (dims == 2) then
            w = ((Xel(ns_sorted) - Xf(ks_sorted))**2 + (Yel(ns_sorted) - Yf(ks_sorted))**2)**(-weight/2.0)
        else
            w = ((Xel(ns_sorted) - Xf(ks_sorted))**2 + (Yel(ns_sorted) - Yf(ks_sorted))**2 + (Zel(ns_sorted) - Zf(ks_sorted))**2)**(-weight/2.0)
        end if
        
        open(21,file='w_early.txt',status='unknown',form='formatted',action='write')
        do i=1,size(w)
            write(21,*)  w(i)
        enddo
        close(21)
        
        ! Prepare distances for the interpolation.       
        allocate(inds2_vals_temp(size(ns)))
        call unique_sort(ks_sorted, inds2_vals_temp, n_unique)
        allocate(inds2_vals, source=inds2_vals_temp(1:n_unique))
        deallocate(inds2_vals_temp)
             
        allocate(inds2(size(inds2_vals)))
        do i = 1, size(inds2_vals)
            inds2(i) = minloc(ks_sorted, 1, mask=ks_sorted .eq. inds2_vals(i), back=.true.)
        end do    
        
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
        
        sjask = 0
        allocate(k_log(K))
        k_log(:) = 0
        do kk = 1, K
            
            allocate(ind(inds2(kk)-inds1(kk)+1))
            k_i = 1
            
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
           
            nns = [NX(kk), NY(kk), NZ(kk)]
            
            if (kk == 1) then
                call displayGUIMessage( 'Saving temporary files, loop 1' )
                
                open(21,file='ind_1.txt',status='unknown',form='formatted',action='write')
                do i=1,size(ind)
                    write(21,*)  ind(i)
                enddo
                close(21)
                
                open(21,file='dxk_test_1.txt',status='unknown',form='formatted',action='write')
                do i=1,size(dxk)
                    write(21,*)  dxk(i)
                enddo
                close(21)
                
                open(21,file='mask1D_1.txt',status='unknown',form='formatted',action='write')
                do i=1,size(mask1D)
                    write(21,*)  mask1D(i)
                enddo
                close(21)
            end if
            
            if (kk == K) then
                call displayGUIMessage( 'Saving temporary files, loop end' )
        
                !write (prog_str,'(I20)') (sum(abs(pack(Signs(:,3),mask1D))))
                !call displayGUIMessage( prog_str )
                
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
                                
            end if
            
            lind = size(ind)
            
            k_log(kk) = sum(abs(pack(Signs(:,3),mask1D)))
            ! Mirror trick to enforce Neumann b.c. Creates a set of virtual nodes on the other side of an edge face.
            if (sum(abs(pack(Signs(:,3),mask1D))) == 1) then
                k_log(kk) = 2
                
                if ((dims .eq. 1 .and. abs(nns(1)) < eps_criteria) .or. (dims .eq. 2 .and. abs(nns(1)) < eps_criteria .and. abs(nns(2)) < eps_criteria)) then
                    k_log(kk) = 3
                    deallocate(ind)
                    
                    cycle
                endif
                
                allocate(extra(size(dxk,1),3))
                extra(:,1) = dxk - 2 * nns(1) * (dxk * nns(1))
                extra(:,2) = dyk - 2 * nns(2) * (dyk * nns(2))
                extra(:,3) = dzk - 2 * nns(3) * (dzk * nns(3))
                dxk2 = [dxk, extra(:,1)]
                dyk2 = [dyk, extra(:,2)]
                dzk2 = [dzk, extra(:,3)]
                Wk2 = [Wk, Wk]               
            
                deallocate(extra)
                
                allocate(e(2*lind))
            else
                dxk2 = dxk
                dyk2 = dyk
                dzk2 = dzk
                Wk2 = Wk
                
                allocate(e(lind))
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
                    
                    open(21,file='HkRed_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(HkRed,1)
                        do j=1,size(HkRed,2)
                            write(21,*)  HkRed(i,j)
                        enddo
                    enddo
                    close(21)
                    
                    open(21,file='Wktmp_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(Wktmp,1)
                        do j=1,size(Wktmp,2)
                            write(21,*)  Wktmp(i,j)
                        enddo
                    enddo
                    close(21)
                                        
                    open(21,file='vw_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(vw)
                        write(21,*)  vw(i)
                    enddo
                    close(21)
                                        
                    open(21,file='vx_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(vx)
                        write(21,*)  vx(i)
                    enddo
                    close(21)
                    
                    open(21,file='vy_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(vy)
                        write(21,*)  vy(i)
                    enddo
                    close(21)
                    
                    open(21,file='vz_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(vz)
                        write(21,*)  vz(i)
                    enddo
                    close(21)
                    
                    close(21)
                    open(21,file='Wk2_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(Wk2)
                        write(21,*)  Wk2(i)
                    enddo
                    close(21)
                    
                    open(21,file='dxk_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dxk)
                        write(21,*)  dxk(i)
                    enddo
                    close(21)
                    
                    open(21,file='dyk_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dyk)
                        write(21,*)  dyk(i)
                    enddo
                    close(21)
                    
                    open(21,file='dzk_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dzk)
                        write(21,*)  dzk(i)
                    enddo
                    close(21)
                    
                    open(21,file='dxk2_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dxk2)
                        write(21,*)  dxk2(i)
                    enddo
                    close(21)
                    
                    open(21,file='dyk2_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dyk2)
                        write(21,*)  dyk2(i)
                    enddo
                    close(21)
                    
                    open(21,file='dzk2_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(dzk2)
                        write(21,*)  dzk2(i)
                    enddo
                    close(21)
                    
                    open(21,file='nns_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(nns)
                        write(21,*)  nns(i)
                    enddo
                    close(21)
                    
                    open(21,file='mask_int2log_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(mask_int2log)
                        write(21,*)  mask_int2log(i)
                    enddo
                    close(21)
                    
                    open(21,file='Gkl1_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(Gkl1,1)
                        write(21,*)  Gkl1(i,1)
                        write(21,*)  Gkl1(i,2)
                        write(21,*)  Gkl1(i,3)
                        write(21,*)  Gkl1(i,4)
                    enddo
                    close(21)
                                        
                    open(21,file='Gk_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(Gk,1)
                        write(21,*)  Gk(i,1)
                        write(21,*)  Gk(i,2)
                        write(21,*)  Gk(i,3)
                        write(21,*)  Gk(i,4)
                    enddo
                    close(21)
                                        
                    open(21,file='Hk_1.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(Hk,1)
                        do j=1,size(Hk,2)
                            write(21,*)  Hk(i,j)
                        enddo
                    enddo
                    close(21)
                                   
                    open(21,file='GkRed_1.txt',status='unknown',form='formatted',action='write')
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
                
        ! Final operation, summing interpolated values according to either ...
        if ( method .eq. MicroMagExchMethodGGNeumann ) then
        !if (method == "GGNeumann") then
            ! ... the Green-Gauss theorem, yielding an estimate for the gradient
            !    W = sparse(ks, ns, vw, K, N)
            call create_CSR_matrix(ks_sorted, ns_sorted, vw, K, N, eps_criteria, W_sparse)
            
            !    DX = DDX * W
            stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDX_sparse, W_sparse, DX_matrix)
            
                    !Debug code to output the matrices
                    call displayGUIMessage( 'Starting DX dense product' )
                    allocate(DX_matrix_dense(N,N))
                    stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDX_sparse, W_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, DX_matrix_dense, N)
         
                    call displayGUIMessage( 'Ending DX dense product' )
         
                    open(21,file='DX_matrix_dense.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(DX_matrix_dense,1)
                        do j=1,size(DX_matrix_dense,2)
                            write(21,*)  DX_matrix_dense(i,j)
                        enddo
                    enddo
                    close(21)
            
            if (dims > 1) then
                !    DY = DDY * W
                stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDY_sparse, W_sparse, DY_matrix)
                
                        !Debug code to output the matrices
                        call displayGUIMessage( 'Starting DY dense product' )
                        allocate(DY_matrix_dense(N,N))
                        stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDY_sparse, W_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, DY_matrix_dense, N)
         
                        call displayGUIMessage( 'Ending DY dense product' )
         
                        open(21,file='DY_matrix_dense.txt',status='unknown',form='formatted',action='write')
                        do i=1,size(DY_matrix_dense,1)
                            do j=1,size(DY_matrix_dense,2)
                                write(21,*)  DY_matrix_dense(i,j)
                            enddo
                        enddo
                        close(21)
                
                if (dims > 2) then
                    !    DZ = DDZ * W
                    stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDZ_sparse, W_sparse, DZ_matrix)
                    
                            !Debug code to output the matrices
                            call displayGUIMessage( 'Starting DZ dense product' )
                            allocate(DZ_matrix_dense(N,N))
                            stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDZ_sparse, W_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, DZ_matrix_dense, N)
         
                            call displayGUIMessage( 'Ending DZ dense product' )
         
                            open(21,file='DZ_matrix_dense.txt',status='unknown',form='formatted',action='write')
                            do i=1,size(DZ_matrix_dense,1)
                                do j=1,size(DZ_matrix_dense,2)
                                    write(21,*)  DZ_matrix_dense(i,j)
                                enddo
                            enddo
                            close(21)
                end if
            end if
        else if ( method .eq. MicroMagExchMethodDirectLaplacianNeumann ) then
        !else if (method == "DirectLaplacianNeumann") then
            ! ... the divergence theorem, yielding an estimate for the Laplacian          
            !    FX = sparse(ks, ns, vx, K, N)
            call create_CSR_matrix(ks_sorted, ns_sorted, vx, K, N, eps_criteria, FX_sparse)
        
                    !Debug code to output the matrices
                    call displayGUIMessage( 'Starting DDXA dense product' )
        
                    ! For debugging, save the sparse matrix as a dense by multiplying it with the identity matrix
                    ! See https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2024-0/mkl-sparse-spmmd.html
                    !allocate(DDXA_matrix_dense(N,K))
                    !DDXA_matrix_dense(:,:) = 0.0
                    !stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDXA_sparse, Identity_matrix_K_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, DDXA_matrix_dense, N)
                    !
                    !call displayGUIMessage( 'Ending DDXA dense product' )
                    !
                    !open(21,file='DDXA_matrix_dense.txt',status='unknown',form='formatted',action='write')
                    !do i=1,size(DDXA_matrix_dense,1)
                    !    do j=1,size(DDXA_matrix_dense,2)
                    !        write(21,*)  DDXA_matrix_dense(i,j)
                    !    enddo
                    !enddo
                    !close(21)
            
            
                    !call displayGUIMessage( 'Starting FX dense product' )
                    !
                    !allocate(FX_matrix_dense(K,N))
                    !FX_matrix_dense(:,:) = 0.0
                    !stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, FX_sparse, Identity_matrix_N_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, FX_matrix_dense, K)
                    !
                    !call displayGUIMessage( 'Ending FX dense product' )
                    !
                    !open(21,file='FX_matrix_dense.txt',status='unknown',form='formatted',action='write')
                    !do i=1,size(FX_matrix_dense,1)
                    !    do j=1,size(FX_matrix_dense,2)
                    !        write(21,*)  FX_matrix_dense(i,j)
                    !    enddo
                    !enddo
                    !close(21)
        
        
                    call displayGUIMessage( 'Starting DX dense product' )
                    allocate(DX_matrix_dense(N,N))
                    stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDXA_sparse, FX_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, DX_matrix_dense, N)
         
                    !call displayGUIMessage( 'STAT' )
                    !write (prog_str,'(I10)') (stat)
                    !call displayGUIMessage( prog_str )
         
                    call displayGUIMessage( 'Ending DX dense product' )
         
                    open(21,file='DX_matrix_dense.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(DX_matrix_dense,1)
                        do j=1,size(DX_matrix_dense,2)
                            write(21,*)  DX_matrix_dense(i,j)
                        enddo
                    enddo
                    close(21)
                    
                    
            !    DX = DDXA * FX
            ! Compute the matrix product of two sparse matrices
            ! See https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2025-1/mkl-sparse-spmm.html
            stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDXA_sparse, FX_sparse, DX_matrix)
        
            !call displayGUIMessage( 'STAT' )
            !write (prog_str,'(I10)') (stat)
            !call displayGUIMessage( prog_str )
        
            call displayGUIMessage( 'Sparse product done' )
        
            if (dims > 1) then
                !FY = sparse(ks, ns, vy, K, N)
                call create_CSR_matrix(ks_sorted, ns_sorted, vy, K, N, eps_criteria, FY_sparse)
                
                    !Debug code to output the matrices
                    call displayGUIMessage( 'Starting DY dense product' )
                    allocate(DY_matrix_dense(N,N))
                    stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDYA_sparse, FY_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, DY_matrix_dense, N)
         
                    call displayGUIMessage( 'Ending DY dense product' )
         
                    open(21,file='DY_matrix_dense.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(DY_matrix_dense,1)
                        do j=1,size(DY_matrix_dense,2)
                            write(21,*)  DY_matrix_dense(i,j)
                        enddo
                    enddo
                    close(21)
                
                !DY = DDYA * FY
                stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDYA_sparse, FY_sparse, DY_matrix)
                if (dims > 2) then
                    !FZ = sparse(ks, ns, vz, K, N)
                    call create_CSR_matrix(ks_sorted, ns_sorted, vz, K, N, eps_criteria, FZ_sparse)
                    
                        !Debug code to output the matrices
                        call displayGUIMessage( 'Starting DZ dense product' )
                        allocate(DZ_matrix_dense(N,N))
                        stat = mkl_sparse_d_spmmd (SPARSE_OPERATION_NON_TRANSPOSE, DDZA_sparse, FZ_sparse, SPARSE_LAYOUT_COLUMN_MAJOR, DZ_matrix_dense, N)
         
                        call displayGUIMessage( 'Ending DZ dense product' )
         
                        open(21,file='DZ_matrix_dense.txt',status='unknown',form='formatted',action='write')
                        do i=1,size(DZ_matrix_dense,1)
                            do j=1,size(DZ_matrix_dense,2)
                                write(21,*)  DZ_matrix_dense(i,j)
                            enddo
                        enddo
                        close(21)
                    
                    !DZ = DDZA * FZ
                    stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDZA_sparse, FZ_sparse, DZ_matrix)
                end if
            end if
        end if

        
        
        
        !Debug code to output the files
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
        
        
    end subroutine computeDifferentialOperatorsFromMesh_DirectLap


    
    
subroutine unique_sort(val, unique_val, n_unique)
    implicit none
    integer :: i, min_val, max_val
    integer, dimension(:), intent(in) :: val 
    integer, dimension(:), intent(out) :: unique_val
    integer, intent(out) :: n_unique
    integer, dimension(:), allocatable :: unique
    logical, allocatable :: mask1D(:)
      
    allocate(unique(size(val)))
    min_val = minval(val)-1
    max_val = maxval(val)
    
    allocate(mask1D(size(val)))
    
    i = 0
    do while (min_val < max_val)
        i = i+1
        mask1D = (val > min_val)
        min_val = minval(val, mask1D)
        unique(i) = min_val
    enddo
    unique_val(1:i) = unique(1:i)
    
    n_unique = i
    
    deallocate(unique,mask1D)
    
end subroutine unique_sort


subroutine create_CSR_matrix(rows, columns, values, K, N, eps_criteria, CSR_matrix)
    implicit none
    integer, dimension(:), intent(in) :: rows, columns
    real(dp), dimension(:), intent(in) :: values
    integer, intent(in) :: K, N
    real(dp), intent(in) :: eps_criteria
    type(sparse_matrix_t), intent(out) :: CSR_matrix
    integer :: nnz, stat
    integer, dimension(:), allocatable :: rows_reduced_COO, columns_reduced_COO
    real(dp), dimension(:), allocatable :: values_reduced_COO
    logical, allocatable :: mask1D(:)
    type(sparse_matrix_t) :: COO_matrix
    
    !Find the non-zero values of values
    allocate(mask1D(size(values)))    
    mask1D(:) = .false.
    mask1D = (abs(values) .gt. eps_criteria)
        
    !Then reduce the size of the arrays
    rows_reduced_COO    = pack(rows,mask1D)
    columns_reduced_COO = pack(columns,mask1D)
    values_reduced_COO  = pack(values,mask1D)
    deallocate(mask1D)

    ! Create a sparse matrix in COO format from the reduced arrays
    ! See https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-1/mkl-sparse-create-coo.html
    nnz = size(values_reduced_COO)
    stat = mkl_sparse_d_create_coo (COO_matrix, SPARSE_INDEX_BASE_ONE, K, N, nnz, rows_reduced_COO, columns_reduced_COO, values_reduced_COO)

    ! Create a sparse matrix in CSR format from COO format
    ! See https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2025-1/mkl-sparse-convert-csr.html
    stat = mkl_sparse_convert_csr (COO_matrix, SPARSE_OPERATION_NON_TRANSPOSE, CSR_matrix)
    
    stat = mkl_sparse_destroy(COO_matrix)
        
end subroutine create_CSR_matrix

end module DifferentialOperators