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

        ! Unpack the variables
        real(dp), dimension(:),allocatable :: NX, NY, NZ, Areas, Volumes, Xel, Yel, Zel, Xf, Yf, Zf
        integer, dimension(:,:),allocatable :: Signs, T, D, el2fa
        real(dp), dimension(:), allocatable :: Aexch_local
        integer :: i,j,dims, n, k, kk, k_mask1D, counter, k_i, lind, sjask, indx, k_j, k_row 
        real(dp), dimension(:),allocatable :: VolCoeff, AX, AY, AZ, s, w
        integer, dimension(:),allocatable :: ss, ns, ks, inds1, inds2, inds2_vals, ks_sorted, ns_sorted
        integer, dimension(:),allocatable :: ns_packed, k_log
        real(dp), dimension(:), allocatable :: Amat
        integer, dimension(:), allocatable :: tmp, tmp2, ind
        logical, allocatable :: mask1D(:), mask_int2log(:)
        real(dp), dimension(:), allocatable :: ddxA, ddyA, ddzA, ddx, ddy, ddz, dx, dy, dz, vx, vy, vz, vw, dks, dxk, dyk, dzk, nns
        real(dp), dimension(:), allocatable :: dxk2, dyk2, dzk2, Wk, Wk2, e, sjaskarr
        real(dp), dimension(:,:), allocatable :: Gkl1, Gk, Hk, Gkl1_temp, Gkl1_T, GkRed, HkRed, Wktmp
        real(dp) :: wm, infinity, scale, scale_local, eps_criteria
        real(dp), allocatable :: extra(:,:)
        character*(40) :: prog_str
        
        
        
          INTEGER          Nmkl, NRHS
          PARAMETER        ( Nmkl = 5, NRHS = 3 )
          INTEGER          LDA, LDB
          PARAMETER        ( LDA = Nmkl, LDB = Nmkl )
          INTEGER          INFO
          INTEGER          IPIV( Nmkl )
          real(dp), allocatable :: A( :, : ), B( :, : )
          
          allocate(A(5,5))
          A(1,:) = [6.80,  -6.05,  -0.45,   8.32,  -9.67]
          A(2,:) = [-2.11,  -3.30,   2.58,   2.71,  -5.14]
          A(3,:) = [5.66,   5.36,  -2.70,   4.35,  -7.26]
          A(4,:) = [5.97,  -4.44,   0.27,  -7.17,   6.08]
          A(5,:) = [8.23,   1.08,   9.04,   2.14,  -6.87]
          
          allocate(B(5,3))
         B(1,:) = [4.02,  -1.56,   9.81]
         B(2,:) = [6.19,   4.00,  -4.09]
         B(3,:) = [-8.22,  -8.67,  -4.57]
         B(4,:) = [-7.57,   1.75,  -8.61]
         B(5,:) = [-3.03,   2.86,   8.99]
         !B(1,:) = [4.02]
         !B(2,:) = [6.19]
         !B(3,:) = [-8.22]
         !B(4,:) = [-7.57]
         !B(5,:) = [-3.03]
         
         CALL DGESV( Nmkl, NRHS, A, LDA, IPIV, B, LDB, INFO )
        
         call displayGUIMessage( 'LAPACK' )
        write (prog_str,'(I10)') (Nmkl)
        call displayGUIMessage( prog_str )
        write (prog_str,'(I10)') (NRHS)
        call displayGUIMessage( prog_str )
        write (prog_str,'(I10)') (LDA)
        call displayGUIMessage( prog_str )
        write (prog_str,'(I10)') (LDB)
        call displayGUIMessage( prog_str )
        write (prog_str,'(I10)') (INFO)
        call displayGUIMessage( prog_str )
        
        open(21,file='B.txt',status='unknown',form='formatted',action='write')
        do i=1,size(B,2)
            write(21,*)  B(1,i)
            write(21,*)  B(2,i)
            write(21,*)  B(3,i)
            write(21,*)  B(4,i)
            write(21,*)  B(5,i)
        enddo
        close(21)
     
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
            
        call displayGUIMessage( 'Test 2-0' )
        
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
        !N = size(Signs, 1)  ! Number of tiles
        !K = size(Signs, 2)  ! Number of faces
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
        !ns = pack([(i, i=1, size(Signs, 1))], mask=(Signs/=0))
        !ks = pack([(j, j=1, size(Signs, 2))], mask=(Signs/=0))

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

            open(21,file='Amat.txt',status='unknown',form='formatted',action='write')
            do i=1,size(Amat)
                write(21,*)  Amat(i)
            enddo
            close(21)
            
            ! Having constructed the exchange values for each face/tile pair, we build the summing matrix.
            !real(wp), dimension(:) :: ddxA, ddyA, ddzA
            allocate(ddxA(size(ns)))
            ddxA = ss * AX(ks) * Amat * VolCoeff(ns)
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
            
            open(21,file='ddxA.txt',status='unknown',form='formatted',action='write')
            do i=1,size(ddxA)
                write(21,*)  ddxA(i)
            enddo
            close(21)

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
        !! Determines which weights are to be used in the first interpolation step
        !if (.not. present(ExtW)) then
        !    [ns,ks] = find_nonzero_indices(el2fa)
        ns = el2fa(:,1)
        ks = el2fa(:,2)
        !    ! Tile volumes may be used, but have not been here. Example:
        !    ! Volumes(ns).*((Xel(ns)-Xf(ks)).^2+(Yel(ns)-Yf(ks)).^2+(Zel(ns)-Zf(ks)).^2).^(-weight/2)
        
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
        
            allocate(w(size(ns)))
            if (dims == 1) then
                !wm = ((Xel(1) - Xf(1))**2)**(-weight/2)
                w = ((Xel(ns_sorted) - Xf(ks_sorted))**2)**(-weight/2)!*(wm**(1.0/weight)) / wm
            else if (dims == 2) then
                !wm = ((Xel(1) - Xf(1))**2 + (Yel(1) - Yf(1))**2)**(-weight/2)
                w = ((Xel(ns_sorted) - Xf(ks_sorted))**2 + (Yel(ns_sorted) - Yf(ks_sorted))**2)**(-weight/2)!*(wm**(1.0/weight)) / wm
            else
                !wm = ((Xel(1) - Xf(1))**2 + (Yel(1) - Yf(1))**2 + (Zel(1) - Zf(1))**2)**(-weight/2)
                w = ((Xel(ns_sorted) - Xf(ks_sorted))**2 + (Yel(ns_sorted) - Yf(ks_sorted))**2 + (Zel(ns_sorted) - Zf(ks_sorted))**2)**(-weight/2)!*(wm**(1.0/weight)) / wm
            end if
        !    W = sparse(ks, ns, w, K, N)
        !else
        !    if (associated(ExtW)) then
        !        [ks,ns] = find_nonzero_indices(el2fa)
        !        w = ExtW(Xel(ns), Yel(ns), Zel(ns), Xf(ks), Yf(ks), Zf(ks))
        !        W = sparse(ks, ns, w, K, N)
        !    else
        !        [ks,ns] = find_nonzero_indices(ExtW)
        !        w = nonzeros(ExtW)
        !        W = ExtW
        !    end if
        !end if

        ! Prepare distances for the interpolation.
        !integer, dimension(:) :: inds1, inds2
        !real(wp), dimension(:) :: dx, dy, dz, vw, vx, vy, vz

         
        
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
        !inds2 = [find_nonzero_indices(diff(ks)), size(ns, 1)]
        !inds1 = [1, inds2 + 1]
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
        !real(wp) :: wm
        
        allocate(mask1D(size(w)))
        
        infinity = HUGE(w) 
        mask1D = w < infinity
        wm = maxval(w,mask1D)
        w = w * (wm**(1.0/weight)) / wm
        
        deallocate(mask1D)
        
        write (prog_str,'(I20)') (size(ns))
        call displayGUIMessage( prog_str )
        
        write (prog_str,'(I20)') (size(w))
        call displayGUIMessage( prog_str )
        
        write (prog_str,'(f20.6)') (wm)
        call displayGUIMessage( prog_str )
        
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

        ! Main loop
        ! The creation of the matrix W for the first step:
        ! W*phi(elements) = phi(faces) [Green Gauss]
        ! W*phi(elements) = dphi(faces) [Direct Laplacian]
        ! Calculated by solving
        ! Gk * [phi(faces);dphi(faces)]_k = Hk ( * phi(elements) )
        ! for each face, ks. Details can be found in [2].
        !integer :: kk, ind, lind, counter
        !real(wp) :: scale
        !real(wp), dimension(:) :: dxk, dyk, dzk, Wk, Gkl1, Gk, Hk, GkRed, HkRed, Wktmp, e, nns

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
            
            !call displayGUIMessage( 'Doing loop number:' )
            write (prog_str,'(I20)') (kk)
            call displayGUIMessage( prog_str )
            
            !write (prog_str,'(I20)') (inds1(kk))
            !call displayGUIMessage( prog_str )
            
            !write (prog_str,'(I20)') (inds2(kk))
            !call displayGUIMessage( prog_str )
                
            do i = inds1(kk), inds2(kk)
                !write (prog_str,'(I20)') (i)
                !call displayGUIMessage( prog_str )
        
                ind(k_i) = i
                k_i = k_i + 1
            enddo
            Wk = w(ind)
            dxk = dx(ind)
            dyk = dy(ind)
            dzk = dz(ind)
            scale_local = sum(abs(dxk))/size(dxk)
            scale = 10.0 ** nint(log10(scale_local))
            
            if (kk .eq. 1) then
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
            endif

            if (dims > 1) then
        !        real(wp), dimension(:) :: dks
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

            !write (prog_str,'(I10)') (sum(abs(pack(Signs(:,2),mask1D))))
            !call displayGUIMessage( prog_str )
            
            counter = counter + 1
            nns = [NX(kk), NY(kk), NZ(kk)]
            
            if (kk == 1) then
                call displayGUIMessage( 'Sum ind 1' )
        
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
                call displayGUIMessage( 'Test loop 1' )
                k_log(kk) = 2
                
                if ((dims .eq. 1 .and. abs(nns(1)) < eps_criteria) .or. (dims .eq. 2 .and. abs(nns(1)) < eps_criteria .and. abs(nns(2)) < eps_criteria)) then
                    k_log(kk) = 3
                    deallocate(ind)
                    cycle
                endif
                !!!CHECK THAT THE SECOND ARGUMENT HERE IS EQUAL TO prod(~nns(1:2)
                
                !call displayGUIMessage( 'Test loop 2' )
                
                allocate(e(2*lind))
                allocate(extra(size(dxk,1),3))
                !extra = [dxk, dyk, dzk] - 2.0 * nns * [dxk, dyk, dzk] * transpose(nns)
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
            
            !call displayGUIMessage( 'Sizes:' )
            
            !write (prog_str,'(I20)') (size(e,1))
            !call displayGUIMessage( prog_str )
            
            !write (prog_str,'(I20)') (size(Gkl1,1))
            !call displayGUIMessage( prog_str )
            
            !write (prog_str,'(I20)') (size(dxk2,1))
            !call displayGUIMessage( prog_str )

            !write (prog_str,'(I20)') (size(Gkl1_temp,1))
            !call displayGUIMessage( prog_str )
            
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
         
                
!            open(21,file='HkRed.txt',status='unknown',form='formatted',action='write')
!            do i=1,size(HkRed,1)
!                do j=1,size(HkRed,2)
!                    write(21,*)  HkRed(i,j)
!                enddo
!            enddo
!            close(21)
                    
            allocate(Wktmp(size(HkRed,1),size(HkRed,2)))
            Wktmp(:,:) = HkRed(:,:)
                
            !call displayGUIMessage( 'Solving linear system' ) 
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
                    call displayGUIMessage( 'Saving files' )
                    
                    open(21,file='Wktmp.txt',status='unknown',form='formatted',action='write')
                    do i=1,size(Wktmp,1)
                        do j=1,size(Wktmp,2)
                            write(21,*)  Wktmp(i,j)
                        enddo
                    enddo
                    close(21)
                    
                    call displayGUIMessage( 'Saving 2' )
                    
                    open(21,file='vw1.txt',status='unknown',form='formatted',action='write')
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
                        write(21,*)  Hk(i,1)
                        write(21,*)  Hk(i,2)
                        write(21,*)  Hk(i,3)
                        write(21,*)  Hk(i,4)
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