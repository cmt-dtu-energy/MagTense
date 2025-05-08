module DifferentialOperators
  use MKL_SPBLAS
  use BLAS95
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
        integer :: i,j,dims, n, k, kk, k_mask1D
        real(dp), dimension(:),allocatable :: VolCoeff, AX, AY, AZ, s, w
        integer, dimension(:),allocatable :: ss, ns, ks, inds1, inds2
        integer, dimension(:),allocatable :: ns_packed
        real(dp), dimension(:), allocatable :: Amat
        integer, dimension(:), allocatable :: tmp, tmp2
        logical, allocatable :: mask1D(:)
        real(dp), dimension(:), allocatable :: ddxA, ddyA, ddzA, ddx, ddy, ddz, dx, dy, dz, vx, vy, vz, vw
        real(dp) :: wm, infinity
        character*(40) :: prog_str
        
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
        N = size(Signs, 1)  ! Number of tiles
        K = size(Signs, 2)  ! Number of faces

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
        
            allocate(w(size(ns)))
            if (dims == 1) then
                !wm = ((Xel(1) - Xf(1))**2)**(-weight/2)
                w = ((Xel(ns) - Xf(ks))**2)**(-weight/2)!*(wm**(1.0/weight)) / wm
            else if (dims == 2) then
                !wm = ((Xel(1) - Xf(1))**2 + (Yel(1) - Yf(1))**2)**(-weight/2)
                w = ((Xel(ns) - Xf(ks))**2 + (Yel(ns) - Yf(ks))**2)**(-weight/2)!*(wm**(1.0/weight)) / wm
            else
                !wm = ((Xel(1) - Xf(1))**2 + (Yel(1) - Yf(1))**2 + (Zel(1) - Zf(1))**2)**(-weight/2)
                w = ((Xel(ns) - Xf(ks))**2 + (Yel(ns) - Yf(ks))**2 + (Zel(ns) - Zf(ks))**2)**(-weight/2)!*(wm**(1.0/weight)) / wm
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

        !allocate(inds2(size(ns)))
        call unique_sort(ks, inds2)
        allocate(inds1(size(inds2)+1))
        inds1(1) = 1
        inds1(2:size(inds1)) = inds2(:) + 1
        !inds2 = [find_nonzero_indices(diff(ks)), size(ns, 1)]
        !inds1 = [1, inds2 + 1]
        dx = Xel(ns) - Xf(ks)
        allocate(vw(size(w)))
        allocate(vx(size(w)))
        dy = Yel(ns) - Yf(ks)
        allocate(vy(size(w)))
        dz = Zel(ns) - Zf(ks)
        allocate(vz(size(w)))

        ! Scale weights to avoid ill conditioning of the least squares interpolation.
        !real(wp) :: wm
        
        infinity = HUGE(w) 
        mask1D = w < infinity
        wm = maxval(w,mask1D)
        w = w * (wm**(1.0/weight)) / wm
        
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

        !counter = 0
        !do kk = 1, K
        !    ind = inds1(kk):inds2(kk)
        !    Wk = w(ind)
        !    dxk = dx(ind)
        !    dyk = dy(ind)
        !    dzk = dz(ind)
        !    scale = 10.0_wp ** nint(log10(meanval(abs(dxk))))

        !    if (dims > 1) then
        !        real(wp), dimension(:) :: dks
        !        dks = [dxk, dyk]
        !        scale = 10.0_wp ** nint(log10(meanval(abs(dks))))
        !        if (dims > 2) then
        !            dks = [dks, dzk]
        !            scale = 10.0_wp ** nint(log10(meanval(abs(dks))))
        !            dzk = dzk / scale
        !        end if
        !        dyk = dyk / scale
        !    end if
        !    dxk = dxk / scale

            ! Mirror trick to enforce Neumann b.c. Creates a set of virtual nodes on the other side of an edge face.
        !    if (sum(abs(Signs(:,kk))) == 1) then
        !        counter = counter + 1
        !        lind = size(ind)
        !        e = ones(2*lind)
        !        nns = [NX(kk), NY(kk), NZ(kk)]
        !        if ((dims == 1 .and. nns(1) == 0) .or. (dims == 2 .and. all(nns(1:2) == 0))) cycle
        !        extra = [dxk, dyk, dzk] - 2.0_wp * nns * transpose([dxk, dyk, dzk]) * nns
        !        dxk = [dxk, extra(:,1)]
        !        dyk = [dyk, extra(:,2)]
        !        dzk = [dzk, extra(:,3)]
        !        Wk = [Wk, Wk]
        !        Gkl1 = [e, dxk, dyk, dzk]
        !        Gk = [Wk * Gkl1, Wk * (dxk * Gkl1), Wk * (dyk * Gkl1), Wk * (dzk * Gkl1)]
        !        Hk = Wk * transpose(Gkl1)

                ! Pick out only the components used in the face-sum.
                ! This means e.g. only x gradient of phi if dims == 1 and/or the norm of the face is in the x-direction.
        !        GkRed = Gk([.true., nns] /= 0, [.true., nns] /= 0)
        !        HkRed = Hk([.true., nns] /= 0, :)

        !        try
        !            R = chol(GkRed)
        !            Wktmp = solve(R, solve(transpose(R), HkRed))
        !        catch
        !            Wktmp = solve(GkRed, HkRed)
        !        end try

        !        vw(ind) = Wktmp(1, 1:lind) + Wktmp(1, lind+1:)
        !        if (nns(1) /= 0) vx(ind) = (Wktmp(2, 1:lind) + Wktmp(2, lind+1:)) / scale
        !        if (nns(2) /= 0) vy(ind) = (Wktmp(3, 1:lind) + Wktmp(3, lind+1:)) / scale
        !        if (nns(3) /= 0) vz(ind) = (Wktmp(end, 1:lind) + Wktmp(end, lind+1:)) / scale

        !    else
        !        e = ones(lind)
        !        nns = [NX(kk), NY(kk), NZ(kk)]
        !        Gkl1 = [e, dxk, dyk, dzk]
        !        Gk = [Wk * Gkl1, Wk * (dxk * Gkl1), Wk * (dyk * Gkl1), Wk * (dzk * Gkl1)]
        !        Hk = Wk * transpose(Gkl1)

        !        GkRed = Gk([.true., nns] /= 0, [.true., nns] /= 0)
        !        HkRed = Hk([.true., nns] /= 0, :)

        !        try
        !            R = chol(GkRed)
        !            Wktmp = solve(R, solve(transpose(R), HkRed))
        !        catch
        !            Wktmp = solve(GkRed, HkRed)
        !        end try

        !        vw(ind) = Wktmp(1, :)
        !        if (nns(1) /= 0) vx(ind) = Wktmp(2, :) / scale
        !        if (nns(2) /= 0) vy(ind) = Wktmp(3, :) / scale
        !        if (nns(3) /= 0) vz(ind) = Wktmp(end, :) / scale

        !    end if
        !end do

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