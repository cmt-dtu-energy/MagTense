module DifferentialOperators
  use MKL_SPBLAS
  use BLAS95
  use LAPACK95
  use MicroMagParameters
  use IO_GENERAL
  use ISO_C_BINDING
  use UTIL_CALL
  use UTIL_MICROMAG

  use sort_mod
  use trace_mod
  
  implicit none

    contains

    !>-----------------------------------------
    !> @author Rasmus Bj�rk, rabj@dtu.dk, DTU, 2025
    !> Original Matlab implementation by Emil J�rgensen
    !> @brief
    !> computeDifferentialOperatorsFromMesh_DirectLap computes the exchange operator for an unstructured mesh
    !> 
    !> References: [1] "Gradient Calculation Methods on Arbitrary Polyhedral Unstructured Meshes for Cell-Centered CFD Solvers"
    !>             [2] "DifferentialOperatorMeshes" in the documentation folder in MagTense git repository
    !> 
    !> Creates differential operators for an unstructured mesh
    !> The mesh is composed by N elements (or tiles, or cells) and K faces 
    !> 
    !> There are two methods in this function: Green Gauss and Direct Laplacian.
    !> For Green Gauss (Direct Laplacian), the goal is to evaluate the (second) 
    !> derivative of a function Phi in the center of each element n from the 
    !> values that Phi assumes in the center of all the elements.
    !> The goal corresponds to the creation of an N times N sparse matrix DX
    !> 
    !> This task is decomposed in two steps:
    !>
    !> 1) The value that (the gradient of) Phi assumes on the center of each face k
    !> is evaluated from the value that Phi assumes in the center of each element
    !> "close" the face k. This could mean e.g. the simple average of the two
    !> elements having the face k on their boundaries (see Eq. 15) 
    !> or the Inverse Distance Weighted (IDW) Face Interpolation among all the elements 
    !> sharing an edge with the face k (see Eq. 16)
    !> 
    !> 2) The (second) derivative of Phi in the center of a given elment n is
    !> evaluated from the value that (the gradient of) Phi assumes on the center
    !> of each face that is on the boundary of the element n (see Eq. 14)
    !> 
    !> For Green Gauss,
    !> The first step is performed by a K times N sparse matrix W:    W*phi(elements) =  phi(faces)
    !> The second step is performed by N times K sparse matrix DDX: DDX*phi(faces)    = dphi(elements)
    !> The matrix DX is thus given by DDX*W: DX*phi(elements) = (DDX*W)*phi(elements) = dphi(elements)
    !>
    !> For Direct Laplacian,
    !> The first step is performed by a K times N sparse matrix W:    W*phi(elements) =  dphi(faces)
    !> The second step is performed by N times K sparse matrix DDXA: DDXA*phi(faces)    = d2phi(elements)
    !> The matrix DX is thus given by DDX*W: DX*phi(elements) = (DDXA*W)*phi(elements) = d2phi(elements)
    !> 
    !> DDXA differs from DDX in that the former takes the exchange stiffness
    !> Aexch into account to form the second part of the operator
    !> div(A grad(phi)),
    !> whereas the latter forms the second part of the operator grad(phi).
    !> 
    !> INPUT ARGUMENTS:
    !>
    !> GridInfo contains several elements:
    !> Volumes,Xel,Yel,Zel, are N times 1 arrays with the corresponding properties 
    !> of the N elements (volumes, and centers)
    !>
    !> Areas,NX,NY,Nz,Xf,Yf,Zf are K times 1 arrays with the corresponding properties
    !> of the K faces (areas, unitg normals, and centers). 
    !>
    !> T is a K times N sparse matrix which is used to create the matrix W.
    !> The entry T(k,n) is 1 if the n-th element shares at least one edge
    !> with the k-th face 
    !> 
    !> D is a K times N sparse matrix which is used to create the matrix W.
    !> The entry D(k,n) is 1 if the n-th element shares at least one vertex
    !> with the k-th face 
    !>
    !> Signs is a N times K sparse matrix.
    !> The entry Signs(n,k) is equal to 1 if the normal to the k-th faces points outwards 
    !> with respect to the n-th element, and -1 if the normal points inwards.
    !> If the k-th face is not on the boundary of the n-th element Signs(n,k) is zero.
    !> 
    !> interpn is the interpolation scheme used to calculate (the gradient of) Phi.
    !>
    !> weight is the weighting scheme used in the interpolation. Any number is
    !> treated as the exponent of reverse distance weighting.
    !>
    !> method can be either Green Gauss (GG) or DirectLaplacian (DL).
    !> 
    !> Jfact_local is an n-vector containing the exchange interaction strength at each tile.
    !> 
    !> OUTPUT ARGUMENTS:
    !>  
    !> Aexch_matrix - The sparse CSR matrix of the exchange coupling
    !---------------------------------------------------------------------------
    subroutine computeDifferentialOperatorsFromMesh_DirectLap(GridInfo, interpn, weight, method, Jfact_local, phase_local, A_int_local, Aexch_matrix)
        type(MicroMagGridInfo), intent(inout) :: GridInfo
        integer, intent(in) :: interpn
        integer, intent(in) :: method
        real(dp), intent(inout) :: weight
        real(dp), dimension(:), intent(in) :: Jfact_local
        integer, dimension(:), intent(in) :: phase_local     !> Material index of each tile
        real(dp), dimension(:,:), intent(in) :: A_int_local  !> Interface exchange per pair of
                                                             !> phases, normalised the same way
                                                             !> as Jfact_local. Negative means
                                                             !> 'use the harmonic mean'.
        type(sparse_matrix_t),intent(out) :: Aexch_matrix
        type(sparse_matrix_t) :: DX_matrix, DY_matrix, DZ_matrix, A2

        ! Unpack the variables
        real(dp), dimension(:),allocatable :: NX, NY, NZ, Areas, Volumes, Xel, Yel, Zel, Xf, Yf, Zf
        integer, dimension(:,:),allocatable :: Signs, D, el2fa
        real(dp), dimension(:), allocatable :: Aexch_local
        integer :: i,j,dims, N, K, kk, k_mask1D, k_i, lind, indx, k_j, k_row, info, n_unique
        real(dp), dimension(:),allocatable :: VolCoeff, AX, AY, AZ, w
        integer, dimension(:),allocatable :: ss, ns, ks, inds1, inds2, inds2_vals, inds2_vals_temp, ks_sorted, ns_sorted
        integer, dimension(:),allocatable :: ns_packed
        real(dp), dimension(:), allocatable :: Amat
        integer, dimension(:), allocatable :: ind, IPIV
        logical, allocatable :: mask1D(:), mask_int2log(:)
        logical, allocatable :: liveNode(:)
        integer :: n_singular
        !--- face -> incidence index, built once in O(M) ---
        integer, allocatable :: fptr(:), felem(:), fsign(:), fcnt(:)
        logical, allocatable :: edgeFace(:)
        real(dp), allocatable :: Amat_face(:)
        integer :: kf, ip, nface_cells, e1, e2
        real(dp) :: asum, aprd, A_face_ovr
        real(dp), dimension(:), allocatable :: ddxA, ddyA, ddzA, ddx, ddy, ddz, dx, dy, dz, vx, vy, vz, vw, dks, dxk, dyk, dzk, nns
        real(dp), dimension(:), allocatable :: dxk2, dyk2, dzk2, Wk, Wk2, e
        real(dp), dimension(:,:), allocatable :: Gkl1, Gk, Hk, Gkl1_temp, Gkl1_T, GkRed, HkRed, Wktmp
        real(dp) :: wm, infinity, scale, scale_local, eps_criteria, const
        real(dp), allocatable :: extra(:,:)
        integer :: stat                                   !> Status value for the various sparse matrix operations 
        character*(40) :: prog_str
        type(MATRIX_DESCR) :: descr_copy
        type(sparse_matrix_t) :: DDXA_sparse, FX_sparse, DDYA_sparse, FY_sparse, DDZA_sparse, FZ_sparse
        type(sparse_matrix_t) :: DDX_sparse, DDY_sparse, DDZ_sparse, W_sparse        
        integer, save :: itimer=0
        !>------------------------------------------------------------------------------
        integer, dimension(:), allocatable :: sorted_indices
        integer :: u 

        call trace%begin( "computeDifferentialOperatorsFromMesh_DirectLap", itimer=itimer )
        
        eps_criteria = 1.0e-12
        const = 1.
        
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
        !GridInfo%TheTs is not read anywhere in this routine - only TheDs (through el2fa)
        !and TheSigns are. Copying it cost a full allocation and copy of a 2 x 29*Nel
        !integer array, about 26 MB at Nel = 110592, for nothing.
        D = GridInfo%TheDs
        Signs = GridInfo%TheSigns
        Aexch_local = Jfact_local
        
        if (method == MicroMagExchMethodGGNeumann .and. ((maxval(Aexch_local)-minval(Aexch_local))/maxval(Aexch_local)) > 1e-6) then
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

        !>-----------------------------------------
        ! Face -> incidence index.
        !
        ! Two loops below used to answer "which cells touch face k?" by scanning the whole
        ! of Signs, once per incidence and once per face. Both were therefore O(M^2) and
        ! O(K*M), and together they were ~97% of this routine, growing as N^2.1. Signs(:,2)
        ! is a face index in [1,K], so one counting sort turns the same question into an
        ! O(1) lookup: the incidences of face kf are fptr(kf) .. fptr(kf+1)-1, holding the
        ! cell in felem and the sign in fsign. Building it is a single O(M) pass.
        allocate(fcnt(K), fptr(K+1))
        fcnt = 0
        do i = 1, size(ks)
            fcnt(ks(i)) = fcnt(ks(i)) + 1
        end do
        fptr(1) = 1
        do i = 1, K
            fptr(i+1) = fptr(i) + fcnt(i)
        end do
        allocate(felem(size(ks)), fsign(size(ks)))
        fcnt = 0
        do i = 1, size(ks)
            kf = ks(i)
            felem(fptr(kf) + fcnt(kf)) = ns(i)
            fsign(fptr(kf) + fcnt(kf)) = ss(i)
            fcnt(kf) = fcnt(kf) + 1
        end do

        ! A face lies on the boundary when exactly one cell claims it. The original test,
        ! sum(abs(Signs(:,3)) over the incidences of the face) == 1, is the same statement
        ! because every sign is +1 or -1.
        allocate(edgeFace(K))
        do kf = 1, K
            indx = 0
            do ip = fptr(kf), fptr(kf+1)-1
                indx = indx + abs(fsign(ip))
            end do
            edgeFace(kf) = (indx .eq. 1)
        end do


        !>-----------------------------------------
        ! Constructing summing matrix according to reference
        ! Constructing N times K sparse matrix DDXA: DDXA*dphi(faces) = d2phi(elements)
        ! This takes the exchange stiffness Jfact_local into account to form the second part of the operator div(A grad(phi))
        if ( method .eq. MicroMagExchMethodDirectLaplacianNeumann ) then
            ! Setting up exchange interaction strength matrix for heterogeneous materials. Also works for homogeneous materials. [2]
            allocate(Amat(size(ns)))

            ! The exchange value depends only on the face, so it is computed once per face
            ! (K of them) and then scattered over the incidences, rather than recomputed
            ! identically for every incidence of the same face.
            allocate(Amat_face(K))
            Amat_face = 0.0_dp
            do kf = 1, K
                nface_cells = fptr(kf+1) - fptr(kf)
                if (nface_cells == 1) then
                    Amat_face(kf) = Aexch_local(felem(fptr(kf)))
                elseif (nface_cells > 2) then
                    ! Reported once per face now, where the scan reported it once per
                    ! incidence of that face. The mesh is invalid either way.
                    write(*,*) 'Warning: more than two cells share a face!'
                elseif (nface_cells == 2) then
                    e1 = felem(fptr(kf))
                    e2 = felem(fptr(kf)+1)

                    ! A face between two different materials can carry an exchange value the
                    ! user specified for that pair, instead of the harmonic mean of the two
                    ! cells. Every internal face joins exactly two cells, so the pair is
                    ! unambiguous however many materials the mesh contains.
                    !
                    ! A cell with A0 = 0 is not magnetic and carries no exchange, so a face
                    ! touching one stays uncoupled whatever the table says. The harmonic mean
                    ! already behaves that way - 2*A1*A2/(A1+A2) vanishes as soon as either
                    ! side is zero - and matching it here keeps the override consistent both
                    ! with the default and with the uniform grid.
                    A_face_ovr = -1.0_dp
                    if ( size(A_int_local,1) .gt. 1 ) then
                        if ( phase_local(e1) .ne. phase_local(e2) .and. &
                             Aexch_local(e1) .gt. 0.0_dp .and. Aexch_local(e2) .gt. 0.0_dp ) &
                            A_face_ovr = A_int_local(phase_local(e1),phase_local(e2))
                    endif
                    if ( A_face_ovr .ge. 0.0_dp ) then
                        Amat_face(kf) = A_face_ovr
                    else
                    ! product() and sum() over the two-element pack, written out. The pack
                    ! preserved the order of Signs, which the counting sort above also does,
                    ! so the two operands are multiplied and added in the same order.
                    aprd = Aexch_local(e1) * Aexch_local(e2)
                    asum = Aexch_local(e1) + Aexch_local(e2)
                    if ( asum .gt. 0.0_dp ) then
                        Amat_face(kf) = 2.0 * aprd / asum
                    else
                        !Both cells are non-magnetic, where the harmonic mean is 0/0 - an
                        !invalid operation under /fpe:0, not a quiet NaN. Two cells with no
                        !exchange are not exchange coupled, so the face carries nothing.
                        Amat_face(kf) = 0.0_dp
                    endif
                    endif
                end if
            end do

            do kk = 1, size(ks)
                Amat(kk) = Amat_face(ks(kk))
            end do
            deallocate(Amat_face)
             
            ! Having constructed the exchange values for each face/tile pair, we build the summing matrix.
            allocate(ddxA(size(ns)))
            ddxA = ss * AX(ks) * Amat * VolCoeff(ns)
            call create_CSR_matrix(ns, ks, ddxA, N, K, eps_criteria, DDXA_sparse)
           
            if (dims > 1) then
                allocate(ddyA(size(ns)))
                ddyA = ss * AY(ks) * Amat * VolCoeff(ns)
                call create_CSR_matrix(ns, ks, ddYA, N, K, eps_criteria, DDYA_sparse)
                
                if (dims > 2) then
                    allocate(ddzA(size(ns)))
                    ddzA = ss * AZ(ks) * Amat * VolCoeff(ns)
                    call create_CSR_matrix(ns, ks, ddzA, N, K, eps_criteria, DDZA_sparse)
                end if
            end if
        
        else if ( method .eq. MicroMagExchMethodGGNeumann ) then
            ! Constructing N times K sparse matrix DDX: DDX*phi(faces) = dphi(elements)
            allocate(ddx(size(ns)))
            ddx = ss * AX(ks) * VolCoeff(ns)
            call create_CSR_matrix(ns, ks, ddx, N, K, eps_criteria, DDX_sparse)
            
            allocate(ddy(size(ns)))
            ddy = ss * AY(ks) * VolCoeff(ns)
            call create_CSR_matrix(ns, ks, ddy, N, K, eps_criteria, DDY_sparse)
                        
            if (dims > 2) then
                allocate(ddz(size(ns)))
                ddz = ss * AZ(ks) * VolCoeff(ns)
                call create_CSR_matrix(ns, ks, ddz, N, K, eps_criteria, DDZ_sparse)
            end if
        end if
        
        if ( interpn .eq. MicroMagExchInterpnExtended ) then
            allocate(el2fa(size(D,1),2))
            el2fa(:,1) = D(:,2)
            el2fa(:,2) = D(:,1)
        elseif ( interpn .eq. MicroMagExchInterpnCompact ) then
            allocate(el2fa(size(Signs,1),3))
            el2fa = Signs
            call displayGUIMessage( 'Warning: untested method: compact' )
        endif

        ! !>-----------------------------------------
        ! ! Calculating weights
        deallocate(ns,ks)
        allocate(ns(size(el2fa,1)),ks(size(el2fa,1)))
        allocate(ks_sorted(size(ks)))
        allocate(ns_sorted(size(ks)))
        ns = el2fa(:,1)
        ks = el2fa(:,2)

        allocate(sorted_indices(size(ks)))
        call argsort(ks, sorted_indices, algo="quicksort")
        call apply_perm(ks, sorted_indices, ks_sorted)
        call apply_perm(ns, sorted_indices, ns_sorted)
        deallocate(sorted_indices)

        ! Distances from the centers of the tiles to the centers of the faces
        dx = Xel(ns_sorted) - Xf(ks_sorted)
        dy = Yel(ns_sorted) - Yf(ks_sorted)
        dz = Zel(ns_sorted) - Zf(ks_sorted)

        ! Along a periodic direction the tiles at the two ends of the domain have been linked together
        ! in the mesh analysis. The distance between such a linked pair of tiles has to be measured
        ! through the boundary, i.e. using the minimum image convention, and not across the entire
        ! domain. A tile in the interpolation stencil of a face is never further away from that face
        ! than the largest tile size, and the mesh analysis has verified that the period is larger
        ! than twice the largest tile, so only the linked pairs are affected by the wrapping below.
        if ( GridInfo%exchPBC(1) ) then
            where ( dx .gt.  0.5*GridInfo%Lper(1) ) dx = dx - GridInfo%Lper(1)
            where ( dx .lt. -0.5*GridInfo%Lper(1) ) dx = dx + GridInfo%Lper(1)
        endif
        if ( GridInfo%exchPBC(2) ) then
            where ( dy .gt.  0.5*GridInfo%Lper(2) ) dy = dy - GridInfo%Lper(2)
            where ( dy .lt. -0.5*GridInfo%Lper(2) ) dy = dy + GridInfo%Lper(2)
        endif
        if ( GridInfo%exchPBC(3) ) then
            where ( dz .gt.  0.5*GridInfo%Lper(3) ) dz = dz - GridInfo%Lper(3)
            where ( dz .lt. -0.5*GridInfo%Lper(3) ) dz = dz + GridInfo%Lper(3)
        endif

       ! A cell with no exchange stiffness must not take part in the interpolation either.
       ! Its own faces already carry no coupling, because the harmonic mean vanishes, but it
       ! would otherwise still sit in the least squares stencil of the faces around it and
       ! pull on their gradient estimates. Zero weight removes it from those stencils exactly:
       ! it drops out of the weighted Gram matrix Gk and its column of Hk is zero, so the
       ! solve returns a zero interpolation weight for it and it never enters the operator.
       ! A mesh whose cells all have A0 > 0 is untouched.
       ! Determines which weights are to be used in the first interpolation step
        allocate(w(size(ns)))
        allocate(liveNode(size(ns)))
        if (dims == 1) then
            w = (dx**2)**(-weight/2.0)
        else if (dims == 2) then
            w = (dx**2 + dy**2)**(-weight/2.0)
        else
            w = (dx**2 + dy**2 + dz**2)**(-weight/2.0)
        end if


        
        do i = 1, size(w)
            liveNode(i) = ( Aexch_local(ns_sorted(i)) .gt. 0.0_dp )
            if ( .not. liveNode(i) ) w(i) = 0.0_dp
        end do

        !>-----------------------------------------
        ! Prepare distances for the interpolation.       
        !-------------------------------------------------------------------------------------------
        ! Count unique values
        !-------------------------------------------------------------------------------------------
        n_unique = 0
        if (size(ks_sorted) > 0) then
            n_unique = 1
            do i = 2, size(ks_sorted)
                if (ks_sorted(i) /= ks_sorted(i-1)) n_unique = n_unique + 1
            end do
        end if
        !-------------------------------------------------------------------------------------------
        allocate(inds2_vals(n_unique))
        allocate(inds2(n_unique))
        !-------------------------------------------------------------------------------------------
        ! Fill unique values and record the last index for each (equivalent to minloc(..., back=.true.))
        !-------------------------------------------------------------------------------------------
        if (n_unique > 0) then
            u = 1
            inds2_vals(u) = ks_sorted(1)

            do i = 2, size(ks_sorted)
                if (ks_sorted(i) /= ks_sorted(i-1)) then
                    inds2(u) = i - 1
                    u = u + 1
                    inds2_vals(u) = ks_sorted(i)
                end if
            end do
            inds2(u) = size(ks_sorted)
        end if

        allocate(inds1(size(inds2)+1))
        inds1(1) = 1
        inds1(2:size(inds1)) = inds2(:) + 1
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
        ! Global distance scaling currently avoided in favour of local distance
        ! scaling (see below). This is a fractional cost of computing time for a
        ! fractional gain of numerical behaviour. Indeed, without scaling at all
        ! the same results should be acchieved, but the least squares problem would
        ! generally involve solving nearly singular matrix systems.
        ! dist_scale=mean(abs([dx;dy;dz]));
        ! dx=dx/dist_scale;
        ! dy=dy/dist_scale;
        ! dz=dz/dist_scale;

        !>-----------------------------------------
        ! Main loop
        ! The creation of the matrix W for the first step:
        ! W*phi(elements) = phi(faces) [Green Gauss]
        ! W*phi(elements) = dphi(faces) [Direct Laplacian]
        ! Calculated by solving
        ! Gk * [phi(faces);dphi(faces)]_k = Hk ( * phi(elements) )
        ! for each face, ks. Details can be found in [2].
        

        n_singular = 0
        do kk = 1, K
            
            allocate(ind(inds2(kk)-inds1(kk)+1))
            k_i = 1
            
            do i = inds1(kk), inds2(kk)
                ind(k_i) = i             ! list of indices for tiles used in interpolation for face kk.
                k_i = k_i + 1
            enddo
          
            Wk = w(ind)                  ! relevant subset of weights
            dxk = dx(ind)                ! relevant subset of x-distances
            dyk = dy(ind)
            dzk = dz(ind)
            ! Local distance scaling (see above for the global alternative). The set of distances
            ! the scale is derived from is assembled once for the actual dimensionality, rather
            ! than computing log10 for x, then for xy, then for xyz and keeping only the last.
            ! That ordering evaluated log10(scale_local) before the dims tests, so a stencil whose
            ! x-distances all vanish - possible for a face whose neighbours share its x coordinate -
            ! raised a divide-by-zero and handed -Inf to nint, whose result is undefined.
            if (dims == 1) then
                dks = dxk
            else if (dims == 2) then
                dks = [dxk, dyk]
            else
                dks = [dxk, dyk, dzk]
            end if
            scale_local = sum(abs(dks))/size(dks)

            if (scale_local > 0.0_dp) then
                scale = 10.0 ** nint(log10(scale_local))
            else
                ! Every distance in the stencil is zero, so there is no length to normalise by and
                ! nothing to gain from scaling. Leaving the distances alone keeps the least squares
                ! system exactly as it was rather than multiplying it by an arbitrary factor.
                scale = 1.0_dp
            end if

            ! Scale distances to avoid ill conditioning
            if (dims > 2) dzk = dzk / scale
            if (dims > 1) dyk = dyk / scale
            dxk = dxk / scale

           
            nns = [NX(kk), NY(kk), NZ(kk)]   ! normal vector at the edge
            
            lind = size(ind)
            
            ! Mirror trick to enforce Neumann b.c. Creates a set of virtual nodes on the other side of an edge face.
            if (edgeFace(kk)) then     ! edge face. 
                
                if ((dims .eq. 1 .and. abs(nns(1)) < eps_criteria) .or. (dims .eq. 2 .and. abs(nns(1)) < eps_criteria .and. abs(nns(2)) < eps_criteria)) then
                    deallocate(ind)
                    
                    cycle    ! If the face is pointing exclusively in an extradimensional direction, ignore it.
                endif        ! This check is only necessary for edge faces.
                
                allocate(extra(size(dxk,1),3))                       ! mirror reflection on face plane
                extra(:,1) = dxk - 2 * nns(1) * (dxk * nns(1))
                extra(:,2) = dyk - 2 * nns(2) * (dyk * nns(2))
                extra(:,3) = dzk - 2 * nns(3) * (dzk * nns(3))
                dxk2 = [dxk, extra(:,1)]                             ! Virtual x-distances appended to the real ones
                dyk2 = [dyk, extra(:,2)]
                dzk2 = [dzk, extra(:,3)]
                Wk2 = [Wk, Wk]                                       ! Virtual nodes inherit their mirror weights.
            
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
                Hk(i,:) = Wk2(:)*Gkl1_T(i,:)                !  least squares vector. See [1] or [2].
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
            
            ! Pick out only the components used in the face-sum.
            ! This means e.g. only x gradient of phi if dims==1 and/or the norm
            ! of the face is in the x-direction.
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
            !Keyed off the allocation status rather than off kk == 1. The loop can reach the cycle
            !above before ever allocating IPIV - an edge face pointing purely in an extradimensional
            !direction is skipped - so if that happened on the first face, the kk /= 1 branch would
            !deallocate an unallocated array. Face ordering happens to put an x-normal face first
            !today, which is the only reason this never fired.
            if (allocated(IPIV)) deallocate(IPIV)
            allocate(IPIV(size(GkRed,1)))
            
            call dgesv( size(GkRed,1), size(Wktmp,2), GkRed, size(GkRed,1), IPIV, Wktmp, size(Wktmp,1), INFO )

            !A face whose stencil has too few live cells leaves the least squares system
            !singular, and dgesv then returns whatever it has. That only happens where every
            !useful neighbour has A0 = 0, in which case the face carries no exchange anyway
            !and these weights are multiplied by zero; stopping here keeps a NaN out of the
            !sparse product rather than salvaging the face.
            if ( INFO .ne. 0 ) then
                n_singular = n_singular + 1
                deallocate(ind,e,Gkl1,Gkl1_temp,Hk,Gk,mask_int2log,GkRed,HkRed,Wktmp)
                cycle
            endif
            
            if (edgeFace(kk)) then
                vw(ind) = Wktmp(1,1:lind)+Wktmp(1,lind+1:size(Wktmp,2))                 ! Interpolated face values
                if (abs(nns(1)) > eps_criteria) then                                    ! Interpolated x-components of face gradients
                    vx(ind)=(Wktmp(2,1:lind)+Wktmp(2,lind+1:size(Wktmp,2)))/scale       ! rescaled
                    if (abs(nns(2)) > eps_criteria) then                                ! Interpolated y-components of face gradients
                        vy(ind) = (Wktmp(3,1:lind)+Wktmp(3,lind+1:size(Wktmp,2)))/scale ! rescaled
                    endif
                elseif (abs(nns(2)) > eps_criteria) then                                ! Interpolated y-components of face gradients
                    vy(ind)=(Wktmp(2,1:lind)+Wktmp(2,lind+1:size(Wktmp,2)))/scale       ! rescaled
                endif
                if (abs(nns(3)) > eps_criteria) then                                    ! Interpolated z-components of face gradients
                    vz(ind)=(Wktmp(size(Wktmp,1),1:lind)+Wktmp(size(Wktmp,1),lind+1:size(Wktmp,2)))/scale ! rescaled
                endif
            else                                            ! NOT an edge face. No mirror tricks necessary.
                vw(ind) = Wktmp(1,:);                       ! Interpolated face values
                if (abs(nns(1)) > eps_criteria) then        ! Interpolated x-components of face gradients
                    vx(ind)=Wktmp(2,:)/scale                ! rescaled
                    if (abs(nns(2)) > eps_criteria) then    ! Interpolated y-components of face gradients
                        vy(ind) = Wktmp(3,:)/scale          ! rescaled
                    endif
                elseif (abs(nns(2)) > eps_criteria) then    ! Interpolated y-components of face gradients
                    vy(ind)=Wktmp(2,:)/scale                ! rescaled
                endif
                if (abs(nns(3)) > eps_criteria) then        ! Interpolated z-components of face gradients
                    vz(ind)=Wktmp(size(Wktmp,1),:)/scale    ! rescaled
                endif
            endif
                
            deallocate(ind,e,Gkl1,Gkl1_temp,Hk,Gk,mask_int2log,GkRed,HkRed,Wktmp)        
        
        end do

        if ( n_singular .gt. 0 ) then
            write (prog_str,'(I10)') n_singular
            call displayGUIMessage( 'Faces left uninterpolated (no live neighbours): '//trim(adjustl(prog_str)) )
        endif
                
        ! Final operation, summing interpolated values according to either ...
        if ( method .eq. MicroMagExchMethodGGNeumann ) then
        !if (method == "GGNeumann") then
            ! ... the Green-Gauss theorem, yielding an estimate for the gradient
            !    W = sparse(ks, ns, vw, K, N)
            call create_CSR_matrix(ks_sorted, ns_sorted, vw, K, N, eps_criteria, W_sparse)
            
            !    DX = DDX * W
            stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDX_sparse, W_sparse, DX_matrix)
            
            if (dims > 1) then
                !    DY = DDY * W
                stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDY_sparse, W_sparse, DY_matrix)
                                
                if (dims > 2) then
                    !    DZ = DDZ * W
                    stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDZ_sparse, W_sparse, DZ_matrix)
                
                end if
            end if
        else if ( method .eq. MicroMagExchMethodDirectLaplacianNeumann ) then
            ! ... the divergence theorem, yielding an estimate for the Laplacian          
            !    FX = sparse(ks, ns, vx, K, N)
            call create_CSR_matrix(ks_sorted, ns_sorted, vx, K, N, eps_criteria, FX_sparse)
                    
            !    DX = DDXA * FX
            ! Compute the matrix product of two sparse matrices
            ! See https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2025-1/mkl-sparse-spmm.html
            stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDXA_sparse, FX_sparse, DX_matrix)
            
            !call displayGUIMessage( 'STAT' )
            !write (prog_str,'(I10)') (stat)
            !call displayGUIMessage( prog_str )
        
            call displayGUIMessage( 'Sparse product done' )
        
            if (dims > 1) then
                call create_CSR_matrix(ks_sorted, ns_sorted, vy, K, N, eps_criteria, FY_sparse)
                                    
                !DY = DDYA * FY
                stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDYA_sparse, FY_sparse, DY_matrix)
                                
                if (dims > 2) then
                    call create_CSR_matrix(ks_sorted, ns_sorted, vz, K, N, eps_criteria, FZ_sparse)
                    
                    !DZ = DDZA * FZ
                    stat = mkl_sparse_spmm (SPARSE_OPERATION_NON_TRANSPOSE, DDZA_sparse, FZ_sparse, DZ_matrix)
                end if
            end if
        end if

        if (dims == 1) then
            descr_copy%type = SPARSE_MATRIX_TYPE_GENERAL 

            stat = mkl_sparse_copy (DX_matrix, descr_copy, Aexch_matrix)
        endif
        if (dims == 2) then
            stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, DX_matrix, const, DY_matrix, Aexch_matrix)
            stat = mkl_sparse_destroy(DY_matrix)
        endif
        if (dims == 3) then
            stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, DX_matrix, const, DY_matrix, A2)
            stat = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, A2, const, DZ_matrix, Aexch_matrix)
            stat = mkl_sparse_destroy(A2)
            stat = mkl_sparse_destroy(DY_matrix)
            stat = mkl_sparse_destroy(DZ_matrix)
        endif
        
        call create_COO_values_from_CSR(Aexch_matrix,GridInfo)
        
        stat = mkl_sparse_destroy(DX_matrix)

    
        call trace%end( "computeDifferentialOperatorsFromMesh_DirectLap", itimer=itimer )
        
    end subroutine computeDifferentialOperatorsFromMesh_DirectLap


    subroutine passDifferentialOperators(problem)
        implicit none
        type(MicroMagProblem),intent(inout) :: problem     !> The problem data structure
        real(dp) :: eps_criteria
    
        eps_criteria = 1.0e-12
        call displayGUIMessage( 'Creating Exchange matrix start' )
        call create_CSR_matrix(problem%grid%A_exch_load%rows, problem%grid%A_exch_load%cols, problem%grid%A_exch_load%values, problem%grid%A_exch_load%nrows, problem%grid%A_exch_load%ncols, eps_criteria, problem%A_exch)
        call displayGUIMessage( 'Creating Exchange matrix end' )
    
    end subroutine passDifferentialOperators

    
    

    
    ! For debugging, save the sparse matrix as a dense by multiplying it with the identity matrix
    ! See https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2024-0/mkl-sparse-spmmd.html

    !real(dp), dimension(:,:), allocatable :: DDXA_matrix_dense, FX_matrix_dense
    !integer, dimension(:),allocatable :: Identity_matrix_K_rows_start_CSR, Identity_matrix_K_rows_end_CSR, Identity_matrix_K_cols
    !real(dp), dimension(:), allocatable :: Identity_matrix_K_values_CSR
    !real(dp), allocatable :: DX_matrix_dense(:,:), DY_matrix_dense(:,:), DZ_matrix_dense(:,:)
    
    !Create identity matrices for K
    !allocate(Identity_matrix_K_rows_start_CSR(K), Identity_matrix_K_rows_end_CSR(K), Identity_matrix_K_cols(K), Identity_matrix_K_values_CSR(K))
    !do i = 1, K
    !    Identity_matrix_K_rows_start_CSR(i) = i
    !    Identity_matrix_K_rows_end_CSR(i) = i+1
    !    Identity_matrix_K_cols(i) = i
    !    Identity_matrix_K_values_CSR(i) = 1.0
    !enddo
    !
    !stat = mkl_sparse_d_create_csr (Identity_matrix_K_sparse, SPARSE_INDEX_BASE_ONE, K, K, Identity_matrix_K_rows_start_CSR, Identity_matrix_K_rows_end_CSR, Identity_matrix_K_cols, Identity_matrix_K_values_CSR)
       
    !Debug code to output the matrices
    !call displayGUIMessage( 'Starting DDXA dense product' )
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
    
    !open(21,file='rows.txt',status='unknown',form='formatted',action='write')
    !    do i=1,size(rows(:))
    !        write(21,*)  rows(i)
    !    enddo
    !    close(21)
    
end module DifferentialOperators
