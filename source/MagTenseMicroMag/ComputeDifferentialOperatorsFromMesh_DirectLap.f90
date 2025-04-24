module DifferentialOperators
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    private
    public :: computeDifferentialOperatorsFromMesh_DirectLap

    type :: GridInfo_Type
        real(wp), dimension(:) :: fNormX, fNormY, fNormZ
        real(wp), dimension(:) :: AreaFaces, Volumes
        integer, dimension(:,:) :: TheSigns, TheTs, TheDs
        real(wp), dimension(:) :: Xel, Yel, Zel, Xf, Yf, Zf
    end type GridInfo_Type

contains

    subroutine computeDifferentialOperatorsFromMesh_DirectLap(DX, DY, DZ, W, GridInfo, interpn, weight, method, Aexch, ExtW)
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

        type(GridInfo_Type), intent(in) :: GridInfo
        character(len=*), intent(in) :: interpn
        real(wp), intent(in) :: weight
        character(len=*), intent(in) :: method
        real(wp), dimension(:), intent(in), optional :: Aexch
        real(wp), dimension(:,:), intent(in), optional :: ExtW
        real(wp), dimension(:,:), intent(out) :: DX, DY, DZ
        real(wp), dimension(:,:), intent(out) :: W

        ! Unpack the variables
        real(wp), dimension(:) :: NX, NY, NZ, Areas, Volumes, Xel, Yel, Zel, Xf, Yf, Zf
        integer, dimension(:,:) :: Signs, T, D
        real(wp), dimension(:), allocatable :: Aexch_local
        integer :: dims, N, K
        real(wp), dimension(:) :: VolCoeff, AX, AY, AZ, s
        integer, dimension(:) :: n, k

        NX = GridInfo%fNormX
        NY = GridInfo%fNormY
        NZ = GridInfo%fNormZ
        Areas = GridInfo%AreaFaces
        Volumes = GridInfo%Volumes
        Signs = GridInfo%TheSigns
        Xel = GridInfo%Xel
        Yel = GridInfo%Yel
        Zel = GridInfo%Zel
        Xf = GridInfo%Xf
        Yf = GridInfo%Yf
        Zf = GridInfo%Zf
        T = GridInfo%TheTs
        D = GridInfo%TheDs

        ! Initialize Aexch if not provided
        if (.not. present(Aexch)) then
            allocate(Aexch_local(size(Xel)))
            Aexch_local = 1.0_wp
        else
            Aexch_local = Aexch
        end if

        ! Error checking
        if (size(Aexch_local) /= size(Xel)) then
            error stop 'Aexch must have length n, where n is the number of tiles'
        end if

        if (method == "GGNeumann" .and. any(Aexch_local /= 1.0_wp)) then
            error stop 'Green-Gauss method not available for heterogeneous exchange stiffness'
        end if

        ! Determine dimensionality
        if (size(unique(GridInfo%Zel)) > 1) then
            dims = 3
        elseif (size(unique(GridInfo%Yel)) > 1) then
            dims = 2
        else
            dims = 1
        end if

        ! Dimensions
        N = size(Signs, 1)  ! Number of tiles
        K = size(Signs, 2)  ! Number of faces

        ! DDX, DDY, DDZ
        VolCoeff = 1.0_wp / Volumes
        AX = NX * Areas
        AY = NY * Areas
        AZ = NZ * Areas
        s = pack(Signs, mask=(Signs/=0))
        n = pack([(i, i=1, size(Signs, 1))], mask=(Signs/=0))
        k = pack([(j, j=1, size(Signs, 2))], mask=(Signs/=0))

        ! Constructing summing matrix according to reference
        ! Constructing N times K sparse matrix DDXA: DDXA*dphi(faces) = d2phi(elements)
        ! This takes the exchange stiffness Aexch into account to form the second part of the operator div(A grad(phi))
        if (method == "DirectLaplacianNeumann") then
            ! Setting up exchange interaction strength matrix for heterogeneous materials. Also works for homogeneous materials. [2]
            real(wp), dimension(:), allocatable :: Amat
            allocate(Amat(size(n)))

            do kk = 1, size(k)
                integer, dimension(:), allocatable :: tmp, tmp2
                tmp = find_nonzero_indices(Signs(:, k(kk)) /= 0)
                tmp2 = pack(tmp, mask=(tmp /= n(kk)))

                if (size(tmp2) == 0) then
                    Amat(kk) = Aexch_local(n(kk))
                elseif (size(tmp2) > 1) then
                    write(*,*) 'Warning: more than two cells share a face!'
                    write(*,*) 'Cells number ', tmp2, ' share face number ', k(kk)
                else
                    Amat(kk) = 2.0_wp * Aexch_local(n(kk)) * Aexch_local(tmp2) / (Aexch_local(n(kk)) + Aexch_local(tmp2))
                    if (Aexch_local(n(kk)) == 0.0_wp .and. Aexch_local(tmp2) == 0.0_wp) then
                        Amat(kk) = 0.0_wp
                    end if
                end if
            end do

            ! Having constructed the exchange values for each face/tile pair, we build the summing matrix.
            real(wp), dimension(:) :: ddxA, ddyA, ddzA
            ddxA = s * AX(k) * Amat * VolCoeff(n)
            DDXA = sparse(n, k, ddxA, N, K)
            if (dims > 1) then
                ddyA = s * AY(k) * Amat * VolCoeff(n)
                DDYA = sparse(n, k, ddyA, N, K)
                if (dims > 2) then
                    ddzA = s * AZ(k) * Amat * VolCoeff(n)
                    DDZA = sparse(n, k, ddzA, N, K)
                end if
            end if

        else if (method == "GGNeumann") then
            ! Constructing N times K sparse matrix DDX: DDX*phi(faces) = dphi(elements)
            real(wp), dimension(:) :: ddx, ddy, ddz
            ddx = s * AX(k) * VolCoeff(n)
            DDX = sparse(n, k, ddx, N, K)
            ddy = s * AY(k) * VolCoeff(n)
            DDY = sparse(n, k, ddy, N, K)
            if (dims > 2) then
                ddz = s * AZ(k) * VolCoeff(n)
                DDZ = sparse(n, k, ddz, N, K)
            end if
        end if

        ! Defaults
        if (.not. present(interpn)) interpn = 'extended'
        if (.not. present(weight)) then
            weight = 8.0_wp
        else if (.not. is_numeric(weight)) then
            weight = real(weight, wp)
            if (.not. allocated(weight)) error stop 'Supplied weight is not a number.'
        end if
        if (.not. present(method)) method = 'DirectLaplacianNeumann'

        ! Interpolation schemes. Determines how many and which neighbours to use for interpolation.
        select case (trim(interpn))
        case ('extended')
            el2fa = transpose(D)
        case ('compact')
            el2fa = Signs
            write(*,*) 'Warning: untested method "compact"'
        case default
            write(*,*) 'Warning: unrecognized interpolation scheme "', interpn, '". Using extended scheme.'
            el2fa = transpose(D)
        end select

        ! Calculating weights
        ! Determines which weights are to be used in the first interpolation step
        if (.not. present(ExtW)) then
            [n,k] = find_nonzero_indices(el2fa)
            ! Tile volumes may be used, but have not been here. Example:
            ! Volumes(n).*((Xel(n)-Xf(k)).^2+(Yel(n)-Yf(k)).^2+(Zel(n)-Zf(k)).^2).^(-weight/2)
            if (dims == 1) then
                w = ((Xel(n) - Xf(k))**2)**(-weight/2)
            else if (dims == 2) then
                w = ((Xel(n) - Xf(k))**2 + (Yel(n) - Yf(k))**2)**(-weight/2)
            else
                w = ((Xel(n) - Xf(k))**2 + (Yel(n) - Yf(k))**2 + (Zel(n) - Zf(k))**2)**(-weight/2)
            end if
            W = sparse(k, n, w, K, N)
        else
            if (associated(ExtW)) then
                [k,n] = find_nonzero_indices(el2fa)
                w = ExtW(Xel(n), Yel(n), Zel(n), Xf(k), Yf(k), Zf(k))
                W = sparse(k, n, w, K, N)
            else
                [k,n] = find_nonzero_indices(ExtW)
                w = nonzeros(ExtW)
                W = ExtW
            end if
        end if

        ! Prepare distances for the interpolation.
        integer, dimension(:) :: inds1, inds2
        real(wp), dimension(:) :: dx, dy, dz, vw, vx, vy, vz

        inds2 = [find_nonzero_indices(diff(k)), size(n, 1)]
        inds1 = [1, inds2 + 1]
        dx = Xel(n) - Xf(k)
        vw = zeros(size(w))
        vx = zeros(size(w))
        dy = Yel(n) - Yf(k)
        vy = zeros(size(w))
        dz = Zel(n) - Zf(k)
        vz = zeros(size(w))

        ! Scale weights to avoid ill conditioning of the least squares interpolation.
        real(wp) :: wm
        wm = maxval(w, mask=(.not. isinf(w)))
        w = w * (wm**(1.0_wp/weight)) / wm

        ! Main loop
        ! The creation of the matrix W for the first step:
        ! W*phi(elements) = phi(faces) [Green Gauss]
        ! W*phi(elements) = dphi(faces) [Direct Laplacian]
        ! Calculated by solving
        ! Gk * [phi(faces);dphi(faces)]_k = Hk ( * phi(elements) )
        ! for each face, k. Details can be found in [2].
        integer :: kk, ind, lind, counter
        real(wp) :: scale
        real(wp), dimension(:) :: dxk, dyk, dzk, Wk, Gkl1, Gk, Hk, GkRed, HkRed, Wktmp, e, nns

        counter = 0
        do kk = 1, K
            ind = inds1(kk):inds2(kk)
            Wk = w(ind)
            dxk = dx(ind)
            dyk = dy(ind)
            dzk = dz(ind)
            scale = 10.0_wp ** nint(log10(meanval(abs(dxk))))

            if (dims > 1) then
                real(wp), dimension(:) :: dks
                dks = [dxk, dyk]
                scale = 10.0_wp ** nint(log10(meanval(abs(dks))))
                if (dims > 2) then
                    dks = [dks, dzk]
                    scale = 10.0_wp ** nint(log10(meanval(abs(dks))))
                    dzk = dzk / scale
                end if
                dyk = dyk / scale
            end if
            dxk = dxk / scale

            ! Mirror trick to enforce Neumann b.c. Creates a set of virtual nodes on the other side of an edge face.
            if (sum(abs(Signs(:,kk))) == 1) then
                counter = counter + 1
                lind = size(ind)
                e = ones(2*lind)
                nns = [NX(kk), NY(kk), NZ(kk)]
                if ((dims == 1 .and. nns(1) == 0) .or. (dims == 2 .and. all(nns(1:2) == 0))) cycle
                extra = [dxk, dyk, dzk] - 2.0_wp * nns * transpose([dxk, dyk, dzk]) * nns
                dxk = [dxk, extra(:,1)]
                dyk = [dyk, extra(:,2)]
                dzk = [dzk, extra(:,3)]
                Wk = [Wk, Wk]
                Gkl1 = [e, dxk, dyk, dzk]
                Gk = [Wk * Gkl1, Wk * (dxk * Gkl1), Wk * (dyk * Gkl1), Wk * (dzk * Gkl1)]
                Hk = Wk * transpose(Gkl1)

                ! Pick out only the components used in the face-sum.
                ! This means e.g. only x gradient of phi if dims == 1 and/or the norm of the face is in the x-direction.
                GkRed = Gk([.true., nns] /= 0, [.true., nns] /= 0)
                HkRed = Hk([.true., nns] /= 0, :)

                try
                    R = chol(GkRed)
                    Wktmp = solve(R, solve(transpose(R), HkRed))
                catch
                    Wktmp = solve(GkRed, HkRed)
                end try

                vw(ind) = Wktmp(1, 1:lind) + Wktmp(1, lind+1:)
                if (nns(1) /= 0) vx(ind) = (Wktmp(2, 1:lind) + Wktmp(2, lind+1:)) / scale
                if (nns(2) /= 0) vy(ind) = (Wktmp(3, 1:lind) + Wktmp(3, lind+1:)) / scale
                if (nns(3) /= 0) vz(ind) = (Wktmp(end, 1:lind) + Wktmp(end, lind+1:)) / scale

            else
                e = ones(lind)
                nns = [NX(kk), NY(kk), NZ(kk)]
                Gkl1 = [e, dxk, dyk, dzk]
                Gk = [Wk * Gkl1, Wk * (dxk * Gkl1), Wk * (dyk * Gkl1), Wk * (dzk * Gkl1)]
                Hk = Wk * transpose(Gkl1)

                GkRed = Gk([.true., nns] /= 0, [.true., nns] /= 0)
                HkRed = Hk([.true., nns] /= 0, :)

                try
                    R = chol(GkRed)
                    Wktmp = solve(R, solve(transpose(R), HkRed))
                catch
                    Wktmp = solve(GkRed, HkRed)
                end try

                vw(ind) = Wktmp(1, :)
                if (nns(1) /= 0) vx(ind) = Wktmp(2, :) / scale
                if (nns(2) /= 0) vy(ind) = Wktmp(3, :) / scale
                if (nns(3) /= 0) vz(ind) = Wktmp(end, :) / scale

            end if
        end do

        ! Final operation, summing interpolated values according to either ...
        if (method == "GGNeumann") then
            ! ... the Green-Gauss theorem, yielding an estimate for the gradient
            W = sparse(k, n, vw, K, N)
            DX = DDX * W
            if (dims > 1) then
                DY = DDY * W
                if (dims > 2) then
                    DZ = DDZ * W
                end if
            end if
        else if (method == "DirectLaplacianNeumann") then
            ! ... the divergence theorem, yielding an estimate for the Laplacian
            FX = sparse(k, n, vx, K, N)
            DX = DDXA * FX
            if (dims > 1) then
                FY = sparse(k, n, vy, K, N)
                DY = DDYA * FY
                if (dims > 2) then
                    FZ = sparse(k, n, vz, K, N)
                    DZ = DDZA * FZ
                end if
            end if
        else
            error stop 'Unrecognized method "', method, '".'
        end if

    end subroutine computeDifferentialOperatorsFromMesh_DirectLap


subroutine find_nonzero_indices(array, indices, count)
    implicit none
    integer, intent(out) :: count
    integer, dimension(:), intent(out) :: indices
    real(wp), dimension(:), intent(in) :: array

    integer :: i, n

    n = size(array)
    count = 0

    ! First, count the number of non-zero elements
    do i = 1, n
        if (array(i) /= 0.0_wp) then
            count = count + 1
        end if
    end do

    ! Allocate the indices array
    allocate(indices(count))

    ! Fill the indices array with the positions of non-zero elements
    count = 0
    do i = 1, n
        if (array(i) /= 0.0_wp) then
            count = count + 1
            indices(count) = i
        end if
    end do
end subroutine find_nonzero_indices


subroutine get_nonzeros(array, nonzeros, count)
    implicit none
    integer, intent(out) :: count
    real(wp), dimension(:), intent(out) :: nonzeros
    real(wp), dimension(:), intent(in) :: array

    integer :: i, n

    n = size(array)
    count = 0

    ! First, count the number of non-zero elements
    do i = 1, n
        if (array(i) /= 0.0_wp) then
            count = count + 1
        end if
    end do

    ! Allocate the nonzeros array
    allocate(nonzeros(count))

    ! Fill the nonzeros array with the non-zero elements
    count = 0
    do i = 1, n
        if (array(i) /= 0.0_wp) then
            count = count + 1
            nonzeros(count) = array(i)
        end if
    end do
end subroutine get_nonzeros

end module DifferentialOperators