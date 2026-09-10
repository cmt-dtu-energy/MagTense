module TetrahedralMeshAnalysis
  use MicroMagParameters
  use IO_GENERAL
  use trace_mod

  implicit none

    contains

    !>-----------------------------------------
    !> @author Rasmus Bjoerk, rabj@dtu.dk, DTU, 2026
    !> Original Matlab implementation (TetrahedralMeshAnalysis.m) by Andrea Roberto Insinga
    !> @brief
    !> TetrahedralUnstructuredMeshAnalysis analyzes a tetrahedral mesh
    !>
    !> The mesh is given by its nodes and by the node-element connectivity, i.e. exactly the
    !> information a mesh generator produces, and it is analyzed purely combinatorially: two
    !> tetrahedra are neighbours when they share three nodes, and the shared face is the face of
    !> each of them opposite the one node that is not shared. This is exact for a conforming mesh
    !> and needs no geometric search at all, which is what lets this routine run in O(N) where the
    !> corresponding routine for unstructured prisms has to fall back on lookup grids.
    !>
    !> GridInfo is filled with the same quantities that CartesianUnstructuredMeshAnalysis produces
    !> for prisms, so the exchange operator in computeDifferentialOperatorsFromMesh_DirectLap
    !> consumes it unchanged. DimsF has no meaning for a triangular face and is returned as zeros.
    !>
    !> Periodic boundary conditions.
    !> If exchPBC is set along a direction, the two boundary planes normal to that direction are
    !> linked together by one of two mechanisms, chosen per direction by whether the mesh itself
    !> is periodic along it.
    !>
    !> 1. Conforming - the surface triangulation on the two planes is a translated copy, which is
    !>    what a mesh generator produces when asked for a periodic mesh (gmsh: 'Periodic Surface',
    !>    COMSOL: a Copy Face mesh operation). Every boundary node then has a partner, the pairs
    !>    are merged into a single canonical node, and the whole analysis runs on the canonical
    !>    connectivity. A pair of tetrahedra on opposite sides of the domain shares three canonical
    !>    nodes and is picked up as an ordinary pair of neighbours, sharing one face with one area
    !>    and one normal. This is exact and involves no geometry beyond pairing the nodes.
    !>
    !> 2. Mortar - the two triangulations do not match, so no pairing exists. Each boundary
    !>    triangle on one plane is then clipped against the triangles it overlaps on the other, and
    !>    every intersection polygon becomes a sub-face shared by the two elements. Because the two
    !>    planes are exactly coplanar once the period is subtracted, the sub-faces tile both parent
    !>    triangles exactly: no gaps, no overlaps, and the flux balance stays exact. The
    !>    interpolation stencil of a sub-face is taken as the union of the stencils of its two
    !>    parent triangles, which is a superset of what a conforming face would give and, unlike a
    !>    geometric point-in-element search, needs no tolerance.
    !>
    !> The mechanism is chosen per direction, so a mesh that is periodic along one direction and
    !> not along another uses the exact path where it can. Downstream nothing distinguishes the
    !> two: computeDifferentialOperatorsFromMesh_DirectLap only needs the distance between a
    !> linked pair of elements to be measured through the boundary, which it already does using
    !> GridInfo%Lper.
    !>
    !> @param[in] nodes 3 x M array with the coordinates of the mesh nodes
    !> @param[in] elements (4 or 10) x N connectivity, 1-based. Only the four corner nodes are
    !>            used, so a quadratic mesh is analyzed as its linear counterpart
    !> @param[inout] GridInfo The grid information produced by the analysis
    !> @param[in] exchPBC Periodic boundary conditions along x, y and z for the exchange coupling
    !---------------------------------------------------------------------------
    subroutine TetrahedralUnstructuredMeshAnalysis( nodes, elements, GridInfo, exchPBC )
    real(dp), intent(in) :: nodes(:,:)
    integer, intent(in) :: elements(:,:)
    type(MicroMagGridInfo), intent(out) :: GridInfo
    integer, intent(in) :: exchPBC(3)

    integer :: N, M, K, i, j, p, q, t, s, loc, i2, totf, nsign, alloc_stat
    integer :: nOut, cnt, mask, nNew, nBndK, thatf, maxNdCnt, nT, nD, capT, capD
    integer :: nMP, nSub, idim, eA, eB, pA, pB, nfn
    integer :: fnodes(3), fnodesC(3), theseC(8)
    real(dp), allocatable :: Xel(:), Yel(:), Zel(:), Volumes(:)
    real(dp), allocatable :: Xf(:), Yf(:), Zf(:)
    real(dp), allocatable :: fNormX(:), fNormY(:), fNormZ(:), AreaFaces(:)
    real(dp), allocatable :: DimsF(:,:)
    integer, allocatable :: AllFaces(:,:), AllFacesC(:,:), nFaceNodes(:)
    integer, allocatable :: elemC(:,:), parent(:), canon(:)
    logical, allocatable :: used(:)
    integer, allocatable :: ndCnt(:), ndStart(:), ndItem(:)
    integer, allocatable :: nbElem(:,:), nbLoc(:,:), nNb(:), elemFace(:,:)
    integer, allocatable :: TheSigns_indices(:,:), TheTs_indices(:,:), TheDs_indices(:,:)
    integer, allocatable :: outElem(:), outMask(:), igrow(:,:)
    !--- mortar interfaces ---
    logical :: mortarDir(3)
    integer, allocatable :: mpElem(:), mpLoc(:), mpDir(:), mpSide(:)
    integer, allocatable :: sfA(:), sfB(:)
    real(dp), allocatable :: sfArea(:), sfCen(:,:)
    logical, allocatable :: isMortarParent(:,:)
    logical :: PBC(3), hasNb(4)
    real(dp) :: globMin(3), globMax(3), Lper(3), maxSpan(3), eMin(3), eMax(3)
    real(dp) :: minEdge, minEdge2, tol, d2, ax(3), bx(3), cx(3), dx(3), cr(3)
    character*(100) :: prog_str
    integer, save :: itimer=0

    call trace%begin( "TetrahedralUnstructuredMeshAnalysis", itimer=itimer )

    call displayGUIMessage( 'Starting tetrahedral mesh analysis' )

    N = size(elements, 2)
    M = size(nodes, 2)

    if ( N .lt. 1 ) then
        call displayGUIMessage( 'MagTense: the tetrahedral mesh holds no elements' )
        error stop 'TetrahedralUnstructuredMeshAnalysis: empty mesh'
    endif
    if ( size(elements,1) .lt. 4 ) then
        call displayGUIMessage( 'MagTense: the element connectivity needs at least four rows' )
        error stop 'TetrahedralUnstructuredMeshAnalysis: malformed connectivity'
    endif
    if ( (minval(elements(1:4,:)) .lt. 1) .or. (maxval(elements(1:4,:)) .gt. M) ) then
        call displayGUIMessage( 'MagTense: the element connectivity refers to a node outside the node array' )
        call displayGUIMessage( 'MagTense: the connectivity has to be 1-based' )
        error stop 'TetrahedralUnstructuredMeshAnalysis: node index out of range'
    endif

    PBC = ( exchPBC .ne. 0 )
    mortarDir = .false.

    !------------------------------------------------------------------------------------------
    ! Element centroids and volumes. The volume is written with abs(), so the orientation of the
    ! connectivity does not matter and an inverted tetrahedron still gets a positive volume.
    !------------------------------------------------------------------------------------------
    allocate(Xel(N), Yel(N), Zel(N), Volumes(N), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'Xel/Yel/Zel/Volumes' )

    do j = 1, N
        ax = nodes(:,elements(1,j))
        bx = nodes(:,elements(2,j))
        cx = nodes(:,elements(3,j))
        dx = nodes(:,elements(4,j))

        Xel(j) = ( ax(1) + bx(1) + cx(1) + dx(1) ) / 4.0_dp
        Yel(j) = ( ax(2) + bx(2) + cx(2) + dx(2) ) / 4.0_dp
        Zel(j) = ( ax(3) + bx(3) + cx(3) + dx(3) ) / 4.0_dp

        cr = crossProduct( bx - dx, cx - dx )
        Volumes(j) = abs( dot_product( ax - dx, cr ) ) / 6.0_dp
    end do

    if ( any( Volumes .le. 0.0_dp ) ) then
        call displayGUIMessage( 'MagTense: the mesh contains a tetrahedron of zero volume' )
        error stop 'TetrahedralUnstructuredMeshAnalysis: degenerate element'
    endif

    !------------------------------------------------------------------------------------------
    ! Bounding box, largest element extent and shortest edge. Only nodes that an element actually
    ! refers to are considered, so a stray unused node cannot inflate the period.
    !------------------------------------------------------------------------------------------
    allocate(used(M), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'used' )
    used = .false.

    globMin =  huge(1.0_dp)
    globMax = -huge(1.0_dp)
    maxSpan = 0.0_dp
    minEdge2 = huge(1.0_dp)
    do j = 1, N
        do p = 1, 4
            used(elements(p,j)) = .true.
        end do

        do i = 1, 3
            eMin(i) = minval( nodes(i,elements(1:4,j)) )
            eMax(i) = maxval( nodes(i,elements(1:4,j)) )
        end do
        globMin = min( globMin, eMin )
        globMax = max( globMax, eMax )
        maxSpan = max( maxSpan, eMax - eMin )

        do p = 1, 3
            do q = p+1, 4
                d2 = sum( ( nodes(:,elements(p,j)) - nodes(:,elements(q,j)) )**2 )
                !The comparison is done on the squared length, so that the initial huge() value
                !is never squared - that overflows, and the build traps on it with /fpe:0
                if ( (d2 .gt. 0.0_dp) .and. (d2 .lt. minEdge2) ) minEdge2 = d2
            end do
        end do
    end do
    Lper = globMax - globMin
    minEdge = sqrt(minEdge2)

    !The tolerance used to pair boundary nodes is tied to the mesh itself rather than to an
    !absolute length, since MagTense works in metres and a mesh is typically nanometre sized.
    tol = 1.0e-6_dp * minEdge

    !------------------------------------------------------------------------------------------
    ! Canonical node numbering. Without periodic boundary conditions every node is its own
    ! representative and everything below reduces to the plain non-periodic analysis.
    !------------------------------------------------------------------------------------------
    allocate(parent(M), canon(M), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'parent/canon' )
    do i = 1, M
        parent(i) = i
    end do

    if ( any(PBC) ) then
        do i = 1, 3
            !The domain has to be more than two elements thick along a periodic direction.
            !Otherwise an element reaches its own periodic image, and the pair of elements linked
            !through the boundary can no longer be identified unambiguously.
            if ( PBC(i) .and. ( Lper(i) .lt. (2.0_dp * maxSpan(i) + 1.0e-9_dp * Lper(i)) ) ) then
                call displayGUIMessage( 'MagTense: too few elements along a periodic direction' )
                call displayGUIMessage( 'MagTense: at least three elements are required for periodic exchange boundary conditions' )
                error stop 'TetrahedralUnstructuredMeshAnalysis: too few elements along a periodic direction'
            endif
        end do

        call buildPeriodicNodeMap( nodes, used, globMin, globMax, PBC, minEdge, tol, parent, mortarDir )
    endif

    do i = 1, M
        canon(i) = ufFind( parent, i )
    end do

    !------------------------------------------------------------------------------------------
    ! Canonical connectivity. Two corners of the same element sharing a canonical node means the
    ! element spans a whole period, which the thickness test above should already have caught.
    !------------------------------------------------------------------------------------------
    allocate(elemC(4,N), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'elemC' )
    do j = 1, N
        do p = 1, 4
            elemC(p,j) = canon( elements(p,j) )
        end do
        do p = 1, 3
            do q = p+1, 4
                if ( elemC(p,j) .eq. elemC(q,j) ) then
                    call displayGUIMessage( 'MagTense: an element touches its own periodic image' )
                    call displayGUIMessage( 'MagTense: the mesh is too thin along a periodic direction' )
                    error stop 'TetrahedralUnstructuredMeshAnalysis: element touches its own image'
                endif
            end do
        end do
    end do

    !------------------------------------------------------------------------------------------
    ! Node -> element lists, in compressed form. This single O(4N) pass replaces the repeated
    ! scans of the whole connectivity array that the Matlab implementation performs once per
    ! element and once per face, and which are what make it O(N^2).
    !
    ! The lists are filled with the elements in ascending order, so each one comes out sorted.
    ! The merges below rely on that: they walk the lists of the nodes of a face in step and
    ! therefore produce the candidate elements in ascending order without any sorting, which is
    ! also the order in which the Matlab implementation emits its entries.
    !------------------------------------------------------------------------------------------
    allocate(ndCnt(M), ndStart(M+1), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'ndCnt/ndStart' )
    ndCnt = 0
    do j = 1, N
        do p = 1, 4
            ndCnt(elemC(p,j)) = ndCnt(elemC(p,j)) + 1
        end do
    end do
    ndStart(1) = 1
    do i = 1, M
        ndStart(i+1) = ndStart(i) + ndCnt(i)
    end do
    allocate(ndItem(4*N), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'ndItem' )
    ndCnt = 0
    do j = 1, N
        do p = 1, 4
            i = elemC(p,j)
            ndItem(ndStart(i) + ndCnt(i)) = j
            ndCnt(i) = ndCnt(i) + 1
        end do
    end do

    !An element can appear at most once in a given node list, so the merge of at most eight node
    !lists - six of them for a mortar sub-face - can never produce more candidates than the total
    !length of those lists
    maxNdCnt = maxval(ndCnt)
    allocate(outElem(8*maxNdCnt), outMask(8*maxNdCnt), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'outElem/outMask' )

    !------------------------------------------------------------------------------------------
    ! Neighbours. Two tetrahedra that share exactly three canonical nodes share the face of each
    ! of them opposite the single node that is not shared, which is the one clear bit position in
    ! the mask returned by the merge. An element shares all four nodes with itself and is
    ! excluded by the same test.
    !------------------------------------------------------------------------------------------
    allocate(nbElem(4,N), nbLoc(4,N), nNb(N), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'nbElem/nbLoc/nNb' )
    nNb = 0

    do j = 1, N
        call mergeElementLists( 4, elemC(:,j), ndStart, ndItem, outElem, outMask, nOut )

        do t = 1, nOut
            mask = outMask(t)
            cnt = popcnt(mask)

            if ( cnt .eq. 4 ) then
                if ( outElem(t) .ne. j ) then
                    call displayGUIMessage( 'MagTense: two elements of the mesh share all four nodes' )
                    error stop 'TetrahedralUnstructuredMeshAnalysis: duplicate element'
                endif
                cycle
            endif
            if ( cnt .ne. 3 ) cycle

            !The local face index is the local node of j that the neighbour does not share
            loc = 0
            do p = 1, 4
                if ( .not. btest(mask, p-1) ) loc = p
            end do

            if ( nNb(j) .ge. 4 ) then
                call displayGUIMessage( 'MagTense: an element has more than four face neighbours' )
                call displayGUIMessage( 'MagTense: the mesh is not conforming - hanging nodes are not supported' )
                error stop 'TetrahedralUnstructuredMeshAnalysis: non-conforming mesh'
            endif
            do s = 1, nNb(j)
                if ( nbLoc(s,j) .eq. loc ) then
                    call displayGUIMessage( 'MagTense: two elements claim the same face of a third element' )
                    call displayGUIMessage( 'MagTense: the mesh is not conforming - hanging nodes are not supported' )
                    error stop 'TetrahedralUnstructuredMeshAnalysis: non-conforming mesh'
                endif
            end do

            nNb(j) = nNb(j) + 1
            nbElem(nNb(j),j) = outElem(t)
            nbLoc(nNb(j),j)  = loc
        end do
    end do

    !------------------------------------------------------------------------------------------
    ! Mortar interfaces, for the periodic directions along which the mesh is not periodic. The
    ! boundary faces on the two planes are clipped against each other and the sub-faces replace
    ! them: the parent faces themselves are skipped in the face construction below, and each
    ! sub-face is appended afterwards as an ordinary internal face with two elements and opposite
    ! signs.
    !------------------------------------------------------------------------------------------
    allocate(isMortarParent(4,N), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'isMortarParent' )
    isMortarParent = .false.
    nMP  = 0
    nSub = 0

    if ( any(mortarDir) ) then
        call buildMortarInterfaces( nodes, elements, nbLoc, nNb, globMin, globMax, mortarDir, &
            minEdge, tol, mpElem, mpLoc, mpDir, mpSide, nMP, sfA, sfB, sfArea, sfCen, nSub )
        do i = 1, nMP
            isMortarParent(mpLoc(i),mpElem(i)) = .true.
        end do
    endif

    !------------------------------------------------------------------------------------------
    ! The number of faces is known exactly once the neighbours and the sub-faces are, so the face
    ! arrays are allocated at their final size rather than at the 4*N upper bound that the Matlab
    ! version uses.
    !------------------------------------------------------------------------------------------
    nNew  = 0
    nBndK = 0
    do j = 1, N
        hasNb = .false.
        do t = 1, nNb(j)
            hasNb(nbLoc(t,j)) = .true.
            if ( nbElem(t,j) .gt. j ) nNew = nNew + 1
        end do
        do loc = 1, 4
            if ( (.not. hasNb(loc)) .and. (.not. isMortarParent(loc,j)) ) nBndK = nBndK + 1
        end do
    end do
    K = nNew + nBndK + nSub

    write(prog_str,'(A,I0,A,I0,A)') 'Tetrahedral mesh: ', N, ' elements, ', K, ' unique faces'
    call displayGUIMessage( trim(prog_str) )
    if ( nSub .gt. 0 ) then
        write(prog_str,'(A,I0,A)') 'Tetrahedral mesh: ', nSub, ' of the faces are mortar sub-faces'
        call displayGUIMessage( trim(prog_str) )
    endif

    allocate(AllFaces(K,3), AllFacesC(K,6), nFaceNodes(K), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'AllFaces/AllFacesC/nFaceNodes' )
    AllFacesC = 0
    allocate(Xf(K), Yf(K), Zf(K), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'Xf/Yf/Zf' )
    allocate(fNormX(K), fNormY(K), fNormZ(K), AreaFaces(K), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'fNormX/fNormY/fNormZ/AreaFaces' )
    allocate(elemFace(4,N), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'elemFace' )
    elemFace = 0

    !Every element contributes one incidence per local face, except that a mortar parent face is
    !replaced by its sub-faces, each of which carries two
    allocate(TheSigns_indices(4*N + 2*nSub,3), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'TheSigns_indices' )

    !------------------------------------------------------------------------------------------
    ! Build the unique faces and the sign matrix.
    !
    ! A face is created by the element with the lower index, which also fixes its orientation:
    ! the normal is flipped to point away from that element, so its own sign is +1 and the sign
    ! of the neighbour, whose normal then points inwards, is -1. A face with no neighbour is a
    ! boundary face and gets +1. This is the convention of the Matlab implementation and of
    ! CartesianUnstructuredMeshAnalysis, and it is what the flux terms downstream assume.
    !
    ! Faces are numbered in the order the Matlab implementation numbers them - elements
    ! ascending, within an element the neighbours ascending, then the boundary faces in ascending
    ! local order - so that the two can be compared entry by entry. The mortar sub-faces have no
    ! counterpart there and are appended at the end.
    !------------------------------------------------------------------------------------------
    totf = 0
    nsign = 0
    do j = 1, N
        hasNb = .false.

        do t = 1, nNb(j)
            i2  = nbElem(t,j)
            loc = nbLoc(t,j)
            hasNb(loc) = .true.

            if ( i2 .gt. j ) then
                totf = totf + 1
                elemFace(loc,j) = totf

                !The same face, seen from the neighbour. Recording it now makes the lookup of an
                !already created face, in the branch below, a single array read rather than the
                !scan of a growing list that the Matlab implementation performs.
                do s = 1, nNb(i2)
                    if ( nbElem(s,i2) .eq. j ) then
                        elemFace(nbLoc(s,i2),i2) = totf
                        exit
                    endif
                end do

                call faceNodesOfLocal( elements(1:4,j), elemC(:,j), loc, fnodes, fnodesC )
                AllFaces(totf,:)    = fnodes
                AllFacesC(totf,1:3) = fnodesC
                nFaceNodes(totf)    = 3
                call computeFaceInfoTriangle( nodes, fnodes, Xel(j), Yel(j), Zel(j), &
                    Xf(totf), Yf(totf), Zf(totf), fNormX(totf), fNormY(totf), fNormZ(totf), AreaFaces(totf) )

                nsign = nsign + 1
                TheSigns_indices(nsign,:) = [ j, totf, 1 ]
            else
                thatf = elemFace(loc,j)
                nsign = nsign + 1
                TheSigns_indices(nsign,:) = [ j, thatf, -1 ]
            endif
        end do

        do loc = 1, 4
            if ( hasNb(loc) ) cycle
            !A mortar parent is not a face of the mesh - its sub-faces are, and they are appended
            !below
            if ( isMortarParent(loc,j) ) cycle

            totf = totf + 1
            elemFace(loc,j) = totf

            call faceNodesOfLocal( elements(1:4,j), elemC(:,j), loc, fnodes, fnodesC )
            AllFaces(totf,:)    = fnodes
            AllFacesC(totf,1:3) = fnodesC
            nFaceNodes(totf)    = 3
            call computeFaceInfoTriangle( nodes, fnodes, Xel(j), Yel(j), Zel(j), &
                Xf(totf), Yf(totf), Zf(totf), fNormX(totf), fNormY(totf), fNormZ(totf), AreaFaces(totf) )

            nsign = nsign + 1
            TheSigns_indices(nsign,:) = [ j, totf, 1 ]
        end do
    end do

    !------------------------------------------------------------------------------------------
    ! Append the mortar sub-faces.
    !
    ! A sub-face is placed on the low plane, so its owner is the element on that side: the stored
    ! normal is the outward normal of that element, which gets +1, and the element on the high
    ! plane gets -1 and therefore sees the opposite normal - exactly the convention of an ordinary
    ! internal face. The sub-face sits a full period away from its second element, which is what
    ! GridInfo%Lper is for downstream.
    !
    ! The node set of a sub-face is the union of the nodes of its two parent triangles. It is used
    ! only to build the interpolation stencil below, where it gives the union of the stencils of
    ! the two parents.
    !------------------------------------------------------------------------------------------
    do i = 1, nSub
        pA = sfA(i)
        pB = sfB(i)
        eA = mpElem(pA)
        eB = mpElem(pB)
        idim = mpDir(pA)

        totf = totf + 1

        call faceNodesOfLocal( elements(1:4,eA), elemC(:,eA), mpLoc(pA), fnodes, fnodesC )
        AllFaces(totf,:)    = fnodes
        AllFacesC(totf,1:3) = fnodesC
        call faceNodesOfLocal( elements(1:4,eB), elemC(:,eB), mpLoc(pB), fnodes, fnodesC )
        AllFacesC(totf,4:6) = fnodesC
        nFaceNodes(totf)    = 6

        Xf(totf) = sfCen(1,i)
        Yf(totf) = sfCen(2,i)
        Zf(totf) = sfCen(3,i)

        !All sub-faces of one plane share the plane normal, pointing out of the low-side element
        fNormX(totf) = 0.0_dp
        fNormY(totf) = 0.0_dp
        fNormZ(totf) = 0.0_dp
        if ( idim .eq. 1 ) fNormX(totf) = -1.0_dp
        if ( idim .eq. 2 ) fNormY(totf) = -1.0_dp
        if ( idim .eq. 3 ) fNormZ(totf) = -1.0_dp

        AreaFaces(totf) = sfArea(i)

        nsign = nsign + 1
        TheSigns_indices(nsign,:) = [ eA, totf,  1 ]
        nsign = nsign + 1
        TheSigns_indices(nsign,:) = [ eB, totf, -1 ]
    end do

    if ( totf .ne. K ) then
        call displayGUIMessage( 'MagTense: internal error - the face count does not match' )
        error stop 'TetrahedralUnstructuredMeshAnalysis: face count mismatch'
    endif

    !------------------------------------------------------------------------------------------
    ! TheTs and TheDs.
    !
    ! TheTs(k,n) is set when element n shares at least two nodes with face k, i.e. an edge, and
    ! TheDs(k,n) when it shares at least one, i.e. a vertex. Both are read straight off the
    ! multiplicity that the merge of the node lists of the face returns, so the three full scans
    ! of the connectivity that the Matlab version does per face disappear entirely.
    !
    ! For a mortar sub-face the six parent nodes are merged instead of three, and the edge test is
    ! applied to each parent triple separately - an element sharing an edge with either parent
    ! shares an edge with the sub-face, whereas one node in each parent is not an edge.
    !
    ! The counts are not known ahead of time, so the buffers start modest and are doubled on
    ! demand. For a typical mesh TheDs holds of the order of fifty entries per face, which makes
    ! it by a wide margin the largest array this routine returns.
    !------------------------------------------------------------------------------------------
    capD = 16 * K
    capT = 8 * K
    allocate(TheDs_indices(capD,2), TheTs_indices(capT,2), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'TheDs_indices/TheTs_indices' )
    nD = 0
    nT = 0

    do i = 1, K
        nfn = nFaceNodes(i)
        theseC(1:nfn) = AllFacesC(i,1:nfn)
        call mergeElementLists( nfn, theseC, ndStart, ndItem, outElem, outMask, nOut )

        do t = 1, nOut
            mask = outMask(t)

            if ( nD .ge. capD ) then
                capD = 2 * capD
                allocate(igrow(capD,2), stat=alloc_stat)
                call checkAllocationTet( alloc_stat, 'TheDs_indices' )
                igrow(1:nD,:) = TheDs_indices(1:nD,:)
                call move_alloc(igrow, TheDs_indices)
            endif
            nD = nD + 1
            TheDs_indices(nD,1) = i
            TheDs_indices(nD,2) = outElem(t)

            !7 selects the bits of the first parent triple and 56 those of the second
            if ( ( popcnt(iand(mask,7)) .ge. 2 ) .or. &
                 ( (nfn .eq. 6) .and. (popcnt(iand(mask,56)) .ge. 2) ) ) then
                if ( nT .ge. capT ) then
                    capT = 2 * capT
                    allocate(igrow(capT,2), stat=alloc_stat)
                    call checkAllocationTet( alloc_stat, 'TheTs_indices' )
                    igrow(1:nT,:) = TheTs_indices(1:nT,:)
                    call move_alloc(igrow, TheTs_indices)
                endif
                nT = nT + 1
                TheTs_indices(nT,1) = i
                TheTs_indices(nT,2) = outElem(t)
            endif
        end do
    end do

    !------------------------------------------------------------------------------------------
    ! Fill GridInfo
    !------------------------------------------------------------------------------------------
    allocate(DimsF(K,3), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'DimsF' )
    !A triangular face has no side lengths. The field is part of the structure that is handed
    !back to Matlab, so it is zero filled rather than left unallocated.
    DimsF = 0.0_dp

    GridInfo%fNormX = fNormX
    GridInfo%fNormY = fNormY
    GridInfo%fNormZ = fNormZ
    GridInfo%AreaFaces = AreaFaces
    GridInfo%Volumes = Volumes
    GridInfo%Xel = Xel
    GridInfo%Yel = Yel
    GridInfo%Zel = Zel
    GridInfo%Xf = Xf
    GridInfo%Yf = Yf
    GridInfo%Zf = Zf
    GridInfo%DimsF = DimsF
    !For a mortar sub-face this holds the nodes of its parent triangle on the low plane, since the
    !sub-face itself is a polygon whose corners are not nodes of the mesh
    GridInfo%AllFaces = AllFaces
    GridInfo%TheSigns = TheSigns_indices(1:nsign,:)
    GridInfo%TheTs = TheTs_indices(1:nT,:)
    GridInfo%TheDs = TheDs_indices(1:nD,:)

    GridInfo%exchPBC = PBC
    GridInfo%Lper = Lper

    call displayGUIMessage( 'Tetrahedral mesh analysis done' )

    call trace%end( "TetrahedralUnstructuredMeshAnalysis", itimer=itimer )

    end subroutine TetrahedralUnstructuredMeshAnalysis


    !>-----------------------------------------
    !> @brief
    !> Walks the element lists of nl nodes in step and returns, in ascending element order, every
    !> element that appears in at least one of them, together with a bit mask of which of the
    !> lists it appeared in. Since the lists are themselves sorted this is a single pass over
    !> them, and it needs neither a sort nor a scratch array of length N.
    !> @param[in] nl The number of node lists to merge, at most 8
    !> @param[in] theseC The canonical node indices whose lists are merged
    !> @param[in] ndStart, ndItem The node -> element lists in compressed form
    !> @param[inout] outElem The elements found, ascending
    !> @param[inout] outMask Bit p-1 is set when the element appeared in the list of node p
    !> @param[out] nOut The number of elements found
    !---------------------------------------------------------------------------
    subroutine mergeElementLists( nl, theseC, ndStart, ndItem, outElem, outMask, nOut )
    integer, intent(in) :: nl
    integer, intent(in) :: theseC(:)
    integer, intent(in) :: ndStart(:), ndItem(:)
    integer, intent(inout) :: outElem(:), outMask(:)
    integer, intent(out) :: nOut

    integer :: ptr(8), fin(8)
    integer :: p, best, mask

    do p = 1, nl
        ptr(p) = ndStart(theseC(p))
        fin(p) = ndStart(theseC(p)+1) - 1
    end do

    nOut = 0
    do
        best = huge(1)
        do p = 1, nl
            if ( ptr(p) .le. fin(p) ) best = min( best, ndItem(ptr(p)) )
        end do
        if ( best .eq. huge(1) ) exit

        mask = 0
        do p = 1, nl
            if ( ptr(p) .le. fin(p) ) then
                if ( ndItem(ptr(p)) .eq. best ) then
                    mask = ibset(mask, p-1)
                    ptr(p) = ptr(p) + 1
                endif
            endif
        end do

        nOut = nOut + 1
        outElem(nOut) = best
        outMask(nOut) = mask
    end do

    end subroutine mergeElementLists


    !>-----------------------------------------
    !> @brief
    !> Returns the three corner nodes of the face of an element opposite its local node loc, both
    !> as the original node indices and as the canonical ones. They are listed in ascending local
    !> order, which is the order the Matlab implementation produces.
    !---------------------------------------------------------------------------
    subroutine faceNodesOfLocal( these, theseC, loc, fnodes, fnodesC )
    integer, intent(in) :: these(4), theseC(4), loc
    integer, intent(out) :: fnodes(3), fnodesC(3)

    integer :: p, n

    n = 0
    do p = 1, 4
        if ( p .eq. loc ) cycle
        n = n + 1
        fnodes(n)  = these(p)
        fnodesC(n) = theseC(p)
    end do

    end subroutine faceNodesOfLocal


    !>-----------------------------------------
    !> @brief
    !> Centroid, unit normal and area of a triangular face. The normal is flipped so that it
    !> points away from the centroid of the element that owns the face, which is what makes the
    !> +1 entry of TheSigns an outward flux for that element.
    !---------------------------------------------------------------------------
    subroutine computeFaceInfoTriangle( nodes, ThoseI, xC, yC, zC, xCf, yCf, zCf, NormX, NormY, NormZ, FaceArea )
    real(dp), intent(in) :: nodes(:,:)
    integer, intent(in) :: ThoseI(3)
    real(dp), intent(in) :: xC, yC, zC
    real(dp), intent(out) :: xCf, yCf, zCf, NormX, NormY, NormZ, FaceArea

    real(dp) :: v1(3), v2(3), v3(3), nrm(3), NormN, ThatProd, ThatSign

    v1 = nodes(:,ThoseI(1))
    v2 = nodes(:,ThoseI(2))
    v3 = nodes(:,ThoseI(3))

    xCf = ( v1(1) + v2(1) + v3(1) ) / 3.0_dp
    yCf = ( v1(2) + v2(2) + v3(2) ) / 3.0_dp
    zCf = ( v1(3) + v2(3) + v3(3) ) / 3.0_dp

    nrm = crossProduct( v2 - v1, v3 - v2 )

    !The length of the cross product is twice the area of the triangle
    NormN = sqrt( nrm(1)**2 + nrm(2)**2 + nrm(3)**2 )
    FaceArea = NormN / 2.0_dp

    ThatProd = (xCf - xC)*nrm(1) + (yCf - yC)*nrm(2) + (zCf - zC)*nrm(3)
    ThatSign = sign( 1.0_dp, ThatProd )

    NormX = nrm(1) * ThatSign / NormN
    NormY = nrm(2) * ThatSign / NormN
    NormZ = nrm(3) * ThatSign / NormN

    end subroutine computeFaceInfoTriangle


    !>-----------------------------------------
    !> @brief
    !> Identifies the nodes on the two boundary planes normal to each periodic direction and
    !> merges every matched pair into one canonical node, through the union-find structure in
    !> parent. The merge is transitive, which is what a node on an edge or a corner of the domain
    !> needs: with two periodic directions such a node has three partners and with three it has
    !> seven, and they all have to end up as a single canonical node rather than as a set of
    !> pairs.
    !>
    !> A direction is only merged when every node on the plane has a partner. If any node does not,
    !> nothing is merged along that direction and it is flagged for mortar coupling instead, so a
    !> half-merged plane can never occur - which is why the partners are collected first and the
    !> merging is done only afterwards.
    !>
    !> The search uses a uniform lookup grid whose cell is at least one element edge across, so
    !> the three by three by three neighbourhood of a query point always covers the tolerance ball
    !> around it and the candidate set is a superset of the true one.
    !> @param[in] nodes 3 x M coordinates
    !> @param[in] used Whether a node is referred to by an element
    !> @param[in] globMin, globMax The bounding box of the mesh
    !> @param[in] PBC The periodic directions
    !> @param[in] minEdge The shortest element edge in the mesh, used as the cell size
    !> @param[in] tol The distance below which two boundary nodes are taken to be the same node
    !> @param[inout] parent The union-find structure holding the canonical node numbering
    !> @param[inout] mortarDir Set for a periodic direction whose two planes do not match
    !---------------------------------------------------------------------------
    subroutine buildPeriodicNodeMap( nodes, used, globMin, globMax, PBC, minEdge, tol, parent, mortarDir )
    real(dp), intent(in) :: nodes(:,:)
    logical, intent(in) :: used(:)
    real(dp), intent(in) :: globMin(3), globMax(3)
    logical, intent(in) :: PBC(3)
    real(dp), intent(in) :: minEdge, tol
    integer, intent(inout) :: parent(:)
    logical, intent(inout) :: mortarDir(3)

    integer :: M, i, idim, alloc_stat
    integer :: nGrid(3), nGridTot, nReg, icell, gi(3), gx, gy, gz, ip, nn
    integer :: nLow, nHigh, nMatch, hit
    integer, allocatable :: cellCnt(:), cellStart(:), cellItem(:)
    integer, allocatable :: highList(:), partnerList(:)
    real(dp) :: gridH(3), gLo(3), gHi(3), pq(3), d2, best2
    real(dp) :: Lper(3)
    character*(100) :: prog_str
    character(1), parameter :: axisName(3) = [ 'x', 'y', 'z' ]

    M = size(nodes,2)
    Lper = globMax - globMin

    !--- the lookup grid over the nodes that an element actually uses -------------------------
    gLo = globMin
    gHi = globMax
    do i = 1, 3
        gridH(i) = minEdge
        if ( gHi(i) .gt. gLo(i) ) then
            nGrid(i) = max(1, int((gHi(i) - gLo(i)) / gridH(i)))
        else
            nGrid(i) = 1
        endif
    end do
    do while ( (real(nGrid(1),dp)*real(nGrid(2),dp)*real(nGrid(3),dp)) .gt. &
               (4.0_dp * real(M,dp) + 1024.0_dp) )
        if ( all(nGrid .le. 1) ) exit
        nGrid = max(1, nGrid / 2)
    end do
    !The cell can only have grown here, so it stays at least one element edge across and the
    !three by three by three neighbourhood keeps covering the tolerance ball
    do i = 1, 3
        if ( gHi(i) .gt. gLo(i) ) then
            gridH(i) = (gHi(i) - gLo(i)) / real(nGrid(i), dp)
        else
            gridH(i) = 1.0_dp
        endif
    end do
    nGridTot = nGrid(1) * nGrid(2) * nGrid(3)

    allocate(cellCnt(nGridTot), cellStart(nGridTot+1), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'cellCnt/cellStart' )
    cellCnt = 0
    do i = 1, M
        if ( .not. used(i) ) cycle
        call cellOf( nodes(:,i), gLo, gridH, nGrid, gi )
        icell = 1 + gi(1) + nGrid(1)*(gi(2) + nGrid(2)*gi(3))
        cellCnt(icell) = cellCnt(icell) + 1
    end do
    cellStart(1) = 1
    do i = 1, nGridTot
        cellStart(i+1) = cellStart(i) + cellCnt(i)
    end do
    nReg = cellStart(nGridTot+1) - 1
    allocate(cellItem(max(1,nReg)), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'cellItem' )
    cellCnt = 0
    do i = 1, M
        if ( .not. used(i) ) cycle
        call cellOf( nodes(:,i), gLo, gridH, nGrid, gi )
        icell = 1 + gi(1) + nGrid(1)*(gi(2) + nGrid(2)*gi(3))
        cellItem(cellStart(icell) + cellCnt(icell)) = i
        cellCnt(icell) = cellCnt(icell) + 1
    end do

    allocate(highList(M), partnerList(M), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'highList/partnerList' )

    !--- pair the two boundary planes of every periodic direction -----------------------------
    do idim = 1, 3
        if ( .not. PBC(idim) ) cycle

        nLow   = 0
        nHigh  = 0
        nMatch = 0

        do i = 1, M
            if ( .not. used(i) ) cycle
            if ( abs(nodes(idim,i) - globMin(idim)) .le. tol ) nLow  = nLow  + 1
            if ( abs(nodes(idim,i) - globMax(idim)) .le. tol ) nHigh = nHigh + 1
        end do

        do i = 1, M
            if ( .not. used(i) ) cycle
            if ( abs(nodes(idim,i) - globMax(idim)) .gt. tol ) cycle

            !The image of this node on the opposite boundary plane
            pq = nodes(:,i)
            pq(idim) = pq(idim) - Lper(idim)

            call cellOf( pq, gLo, gridH, nGrid, gi )

            hit = 0
            best2 = huge(1.0_dp)
            do gz = max(0,gi(3)-1), min(nGrid(3)-1,gi(3)+1)
              do gy = max(0,gi(2)-1), min(nGrid(2)-1,gi(2)+1)
                do gx = max(0,gi(1)-1), min(nGrid(1)-1,gi(1)+1)
                  icell = 1 + gx + nGrid(1)*(gy + nGrid(2)*gz)
                  do ip = cellStart(icell), cellStart(icell+1)-1
                      nn = cellItem(ip)
                      d2 = sum( (nodes(:,nn) - pq)**2 )
                      if ( d2 .le. tol*tol ) then
                          if ( (hit .ne. 0) .and. (d2 .ge. best2) ) cycle
                          hit = nn
                          best2 = d2
                      endif
                  end do
                end do
              end do
            end do

            if ( hit .ne. 0 ) then
                nMatch = nMatch + 1
                highList(nMatch)    = i
                partnerList(nMatch) = hit
            endif
        end do

        !Only merge when the two planes match completely. A partial merge would leave part of the
        !boundary linked and part of it open, which is worse than either mechanism on its own.
        if ( (nMatch .eq. nHigh) .and. (nLow .eq. nHigh) ) then
            do i = 1, nMatch
                call ufUnion( parent, highList(i), partnerList(i) )
            end do
            write(prog_str,'(A,A,A,I0,A)') 'Periodic exchange along ', axisName(idim), &
                ': mesh is periodic, ', nMatch, ' node pairs merged'
            call displayGUIMessage( trim(prog_str) )
        else
            !Mortar coupling is a fallback and not a free one, so it is announced as a warning
            !rather than as a status line: the user asked for a periodic boundary and is getting
            !a different, less accurate mechanism than the one they would get from a mesh that is
            !actually periodic, and they can only act on that if they are told.
            mortarDir(idim) = .true.
            call displayGUIMessage( '*** MagTense WARNING: the mesh is not periodic ***' )
            write(prog_str,'(A,A,A)') 'Periodic boundary conditions were requested along ', &
                axisName(idim), ', but the mesh'
            call displayGUIMessage( trim(prog_str) )
            write(prog_str,'(A,I0,A,I0,A)') 'does not match across it: the two planes hold ', &
                nLow, ' and ', nHigh, ' nodes,'
            call displayGUIMessage( trim(prog_str) )
            write(prog_str,'(A,I0,A)') 'and only ', nMatch, ' of them could be paired up.'
            call displayGUIMessage( trim(prog_str) )
            call displayGUIMessage( 'Falling back to MORTAR COUPLING of the boundary faces.' )
            call displayGUIMessage( 'The exchange operator stays conservative, but it is less' )
            call displayGUIMessage( 'accurate on that plane than in the interior of the mesh.' )
            write(prog_str,'(A,A,A)') 'For full accuracy, mesh the geometry periodically along ', &
                axisName(idim), ':'
            call displayGUIMessage( trim(prog_str) )
            call displayGUIMessage( 'gmsh ''Periodic Surface'', or a COMSOL ''Copy Face'' operation.' )
            call displayGUIMessage( '**************************************************' )
        endif
    end do

    deallocate(cellCnt, cellStart, cellItem, highList, partnerList)

    end subroutine buildPeriodicNodeMap


    !>-----------------------------------------
    !> @brief
    !> Builds the mortar interfaces for the periodic directions whose two boundary planes do not
    !> match node for node.
    !>
    !> The boundary faces on the two planes of such a direction are collected, and every face on
    !> the high plane is clipped against the faces it overlaps on the low plane. Each intersection
    !> polygon becomes a sub-face carrying its own area and centroid, shared by the element on
    !> either side. Because the two planes coincide exactly once the period is subtracted, the
    !> clipping is a plain two dimensional intersection of two triangles in the plane of the
    !> boundary and the sub-faces tile both parents exactly - the flux leaving one element through
    !> a parent face is the sum over its sub-faces and equals what enters on the other side.
    !>
    !> That exactness is checked: the summed sub-face area of every parent has to reproduce the
    !> area of the parent triangle. It does not when the two planes cover different regions, i.e.
    !> when the geometry rather than just the mesh is non-periodic, and that is an error - unlike
    !> a merely non-matching mesh, it cannot be coupled at all.
    !> @param[in] nodes, elements The mesh
    !> @param[in] nbLoc, nNb The face neighbours, used to find the boundary faces
    !> @param[in] globMin, globMax The bounding box of the mesh
    !> @param[in] mortarDir The directions to be mortar coupled
    !> @param[in] minEdge, tol Length scale and tolerance of the mesh
    !> @param[out] mpElem, mpLoc, mpDir, mpSide The parent faces: element, local face index,
    !>             direction, and 1 for the low plane or 2 for the high plane
    !> @param[out] nMP The number of parent faces
    !> @param[out] sfA, sfB The two parent faces of each sub-face, as indices into mpElem
    !> @param[out] sfArea, sfCen The area and the centroid of each sub-face
    !> @param[out] nSub The number of sub-faces
    !---------------------------------------------------------------------------
    subroutine buildMortarInterfaces( nodes, elements, nbLoc, nNb, globMin, globMax, mortarDir, &
        minEdge, tol, mpElem, mpLoc, mpDir, mpSide, nMP, sfA, sfB, sfArea, sfCen, nSub )
    real(dp), intent(in) :: nodes(:,:)
    integer, intent(in) :: elements(:,:)
    integer, intent(in) :: nbLoc(:,:), nNb(:)
    real(dp), intent(in) :: globMin(3), globMax(3)
    logical, intent(in) :: mortarDir(3)
    real(dp), intent(in) :: minEdge, tol
    integer, allocatable, intent(out) :: mpElem(:), mpLoc(:), mpDir(:), mpSide(:)
    integer, intent(out) :: nMP
    integer, allocatable, intent(out) :: sfA(:), sfB(:)
    real(dp), allocatable, intent(out) :: sfArea(:), sfCen(:,:)
    integer, intent(out) :: nSub

    integer :: N, j, loc, t, p, i, idim, iu, iv, alloc_stat, capS
    integer :: nA, nB, ia, ib, npoly, icell, gx, gy, gj0(2), gj1(2), ip
    integer :: nGrid2(2), nGridTot2, nReg
    integer, allocatable :: aList(:), bList(:)
    integer, allocatable :: cellCnt(:), cellStart(:), cellItem(:)
    integer, allocatable :: igrow(:), mpNode(:,:)
    real(dp), allocatable :: rgrow(:), rgrow2(:,:), parentArea(:), coveredArea(:)
    real(dp) :: triA(2,3), triB(2,3), poly(2,12), cen2(2)
    real(dp) :: gLo2(2), gHi2(2), gridH2(2), bb0(2), bb1(2)
    real(dp) :: area, tolLen, plane, relErr, worstErr
    integer :: fnodes(3), fnodesC(3)
    logical :: hasNb(4), onLow, onHigh
    character*(100) :: prog_str

    N = size(elements,2)
    tolLen = 1.0e-9_dp * minEdge

    !------------------------------------------------------------------------------------------
    ! Collect the parent faces: boundary faces whose three nodes all lie on one of the two planes
    ! of a mortar direction. A triangle lies in at most one of those planes, since any two of them
    ! are either parallel or perpendicular, so there is no ambiguity. The first pass counts and
    ! the second fills.
    !------------------------------------------------------------------------------------------
    nMP = 0
    do p = 1, 2
        if ( p .eq. 2 ) then
            allocate(mpElem(max(1,nMP)), mpLoc(max(1,nMP)), mpDir(max(1,nMP)), mpSide(max(1,nMP)), &
                mpNode(3,max(1,nMP)), stat=alloc_stat)
            call checkAllocationTet( alloc_stat, 'mpElem/mpLoc/mpDir/mpSide/mpNode' )
            nMP = 0
        endif

        do j = 1, N
            hasNb = .false.
            do t = 1, nNb(j)
                hasNb(nbLoc(t,j)) = .true.
            end do

            do loc = 1, 4
                if ( hasNb(loc) ) cycle

                call faceNodesOfLocal( elements(1:4,j), elements(1:4,j), loc, fnodes, fnodesC )

                do idim = 1, 3
                    if ( .not. mortarDir(idim) ) cycle

                    onLow  = all( abs(nodes(idim,fnodes) - globMin(idim)) .le. tol )
                    onHigh = all( abs(nodes(idim,fnodes) - globMax(idim)) .le. tol )
                    if ( .not. (onLow .or. onHigh) ) cycle

                    nMP = nMP + 1
                    if ( p .eq. 2 ) then
                        mpElem(nMP) = j
                        mpLoc(nMP)  = loc
                        mpDir(nMP)  = idim
                        if ( onLow ) then
                            mpSide(nMP) = 1
                        else
                            mpSide(nMP) = 2
                        endif
                        mpNode(:,nMP) = fnodes
                    endif
                    exit
                end do
            end do
        end do
    end do

    allocate(parentArea(max(1,nMP)), coveredArea(max(1,nMP)), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'parentArea/coveredArea' )
    parentArea  = 0.0_dp
    coveredArea = 0.0_dp

    capS = max(64, 4*nMP)
    allocate(sfA(capS), sfB(capS), sfArea(capS), sfCen(3,capS), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'sfA/sfB/sfArea/sfCen' )
    nSub = 0

    allocate(aList(max(1,nMP)), bList(max(1,nMP)), stat=alloc_stat)
    call checkAllocationTet( alloc_stat, 'aList/bList' )

    !------------------------------------------------------------------------------------------
    ! Clip, one direction at a time
    !------------------------------------------------------------------------------------------
    do idim = 1, 3
        if ( .not. mortarDir(idim) ) cycle

        !The two in-plane axes, taken cyclically
        iu = mod(idim,3) + 1
        iv = mod(idim+1,3) + 1

        nA = 0
        nB = 0
        do i = 1, nMP
            if ( mpDir(i) .ne. idim ) cycle
            call triOf( nodes, mpNode(:,i), iu, iv, triA )
            parentArea(i) = abs(triArea2D(triA))
            if ( parentArea(i) .le. 0.0_dp ) then
                call displayGUIMessage( 'MagTense: a boundary triangle on a periodic plane has zero area' )
                error stop 'TetrahedralUnstructuredMeshAnalysis: degenerate boundary face'
            endif
            if ( mpSide(i) .eq. 1 ) then
                nA = nA + 1
                aList(nA) = i
            else
                nB = nB + 1
                bList(nB) = i
            endif
        end do

        if ( (nA .eq. 0) .or. (nB .eq. 0) ) then
            call displayGUIMessage( 'MagTense: a periodic boundary plane carries no faces' )
            error stop 'TetrahedralUnstructuredMeshAnalysis: empty periodic boundary plane'
        endif

        !--- a two dimensional lookup grid over the low-plane triangles ------------------------
        gLo2 = [ globMin(iu), globMin(iv) ]
        gHi2 = [ globMax(iu), globMax(iv) ]
        do i = 1, 2
            gridH2(i) = minEdge
            if ( gHi2(i) .gt. gLo2(i) ) then
                nGrid2(i) = max(1, int((gHi2(i) - gLo2(i)) / gridH2(i)))
            else
                nGrid2(i) = 1
            endif
        end do
        do while ( (real(nGrid2(1),dp)*real(nGrid2(2),dp)) .gt. (4.0_dp*real(nA,dp) + 1024.0_dp) )
            if ( all(nGrid2 .le. 1) ) exit
            nGrid2 = max(1, nGrid2 / 2)
        end do
        do i = 1, 2
            if ( gHi2(i) .gt. gLo2(i) ) then
                gridH2(i) = (gHi2(i) - gLo2(i)) / real(nGrid2(i), dp)
            else
                gridH2(i) = 1.0_dp
            endif
        end do
        nGridTot2 = nGrid2(1) * nGrid2(2)

        allocate(cellCnt(nGridTot2), cellStart(nGridTot2+1), stat=alloc_stat)
        call checkAllocationTet( alloc_stat, 'cellCnt/cellStart' )
        cellCnt = 0
        do p = 1, 2
            if ( p .eq. 2 ) then
                cellStart(1) = 1
                do i = 1, nGridTot2
                    cellStart(i+1) = cellStart(i) + cellCnt(i)
                end do
                nReg = cellStart(nGridTot2+1) - 1
                allocate(cellItem(max(1,nReg)), stat=alloc_stat)
                call checkAllocationTet( alloc_stat, 'cellItem' )
                cellCnt = 0
            endif

            do ia = 1, nA
                call triOf( nodes, mpNode(:,aList(ia)), iu, iv, triA )
                bb0 = [ minval(triA(1,:)), minval(triA(2,:)) ] - tolLen
                bb1 = [ maxval(triA(1,:)), maxval(triA(2,:)) ] + tolLen
                call cellRange2D( bb0, bb1, gLo2, gridH2, nGrid2, gj0, gj1 )
                do gy = gj0(2), gj1(2)
                  do gx = gj0(1), gj1(1)
                    icell = 1 + gx + nGrid2(1)*gy
                    if ( p .eq. 1 ) then
                        cellCnt(icell) = cellCnt(icell) + 1
                    else
                        cellItem(cellStart(icell) + cellCnt(icell)) = aList(ia)
                        cellCnt(icell) = cellCnt(icell) + 1
                    endif
                  end do
                end do
            end do
        end do

        !--- clip every high-plane triangle against the low-plane triangles it overlaps --------
        plane = globMin(idim)
        do ib = 1, nB
            i = bList(ib)
            call triOf( nodes, mpNode(:,i), iu, iv, triB )
            bb0 = [ minval(triB(1,:)), minval(triB(2,:)) ] - tolLen
            bb1 = [ maxval(triB(1,:)), maxval(triB(2,:)) ] + tolLen
            call cellRange2D( bb0, bb1, gLo2, gridH2, nGrid2, gj0, gj1 )

            do gy = gj0(2), gj1(2)
              do gx = gj0(1), gj1(1)
                icell = 1 + gx + nGrid2(1)*gy
                do ip = cellStart(icell), cellStart(icell+1)-1
                    ia = cellItem(ip)

                    call triOf( nodes, mpNode(:,ia), iu, iv, triA )

                    !A triangle registers in every cell its bounding box touches, so the same pair
                    !can be reached from more than one cell. Counting it twice would double its
                    !flux, so only the cell holding the lower left corner of the overlap of the
                    !two bounding boxes is allowed to produce it.
                    if ( .not. ownsOverlap( triA, triB, gLo2, gridH2, nGrid2, gx, gy ) ) cycle

                    call clipTriangles( triA, triB, tolLen, poly, npoly )
                    if ( npoly .lt. 3 ) cycle

                    call polyAreaCentroid( poly, npoly, area, cen2 )
                    !Two triangles that only touch along an edge give a polygon of zero area,
                    !which is not a face. The threshold is far below any real sub-face and far
                    !above the rounding of the clip, so the tiling stays exact.
                    if ( area .le. 1.0e-12_dp * min(parentArea(ia), parentArea(i)) ) cycle

                    if ( nSub .ge. capS ) then
                        capS = 2 * capS
                        allocate(igrow(capS), stat=alloc_stat)
                        call checkAllocationTet( alloc_stat, 'sfA' )
                        igrow(1:nSub) = sfA(1:nSub)
                        call move_alloc(igrow, sfA)
                        allocate(igrow(capS), stat=alloc_stat)
                        call checkAllocationTet( alloc_stat, 'sfB' )
                        igrow(1:nSub) = sfB(1:nSub)
                        call move_alloc(igrow, sfB)
                        allocate(rgrow(capS), stat=alloc_stat)
                        call checkAllocationTet( alloc_stat, 'sfArea' )
                        rgrow(1:nSub) = sfArea(1:nSub)
                        call move_alloc(rgrow, sfArea)
                        allocate(rgrow2(3,capS), stat=alloc_stat)
                        call checkAllocationTet( alloc_stat, 'sfCen' )
                        rgrow2(:,1:nSub) = sfCen(:,1:nSub)
                        call move_alloc(rgrow2, sfCen)
                    endif

                    nSub = nSub + 1
                    sfA(nSub) = ia
                    sfB(nSub) = i
                    sfArea(nSub) = area
                    !The sub-face is placed on the low plane, i.e. on the side of its owner
                    sfCen(idim,nSub) = plane
                    sfCen(iu,nSub)   = cen2(1)
                    sfCen(iv,nSub)   = cen2(2)

                    coveredArea(ia) = coveredArea(ia) + area
                    coveredArea(i)  = coveredArea(i)  + area
                end do
              end do
            end do
        end do

        deallocate(cellCnt, cellStart, cellItem)
    end do

    !------------------------------------------------------------------------------------------
    ! Every parent has to be covered by its sub-faces. Anything else means the two planes do not
    ! describe the same cross section, i.e. the geometry itself is not periodic.
    !------------------------------------------------------------------------------------------
    worstErr = 0.0_dp
    do i = 1, nMP
        relErr = abs( coveredArea(i) - parentArea(i) ) / parentArea(i)
        worstErr = max( worstErr, relErr )
    end do

    if ( worstErr .gt. 1.0e-6_dp ) then
        write(prog_str,'(A, ES12.4)') 'MagTense: worst relative area mismatch on a periodic plane: ', worstErr
        call displayGUIMessage( trim(prog_str) )
        call displayGUIMessage( 'MagTense: the two boundary planes do not cover the same cross section' )
        call displayGUIMessage( 'MagTense: the geometry, not just the mesh, has to be periodic' )
        error stop 'TetrahedralUnstructuredMeshAnalysis: non-periodic geometry'
    endif

    write(prog_str,'(A,I0,A,I0,A)') 'Mortar coupling: ', nMP, ' boundary faces split into ', &
        nSub, ' sub-faces'
    call displayGUIMessage( trim(prog_str) )
    !The sub-faces tile their parents exactly, so this is a rounding level number. Anything
    !larger would have been rejected above, and it is reported so that it can be seen to be small.
    write(prog_str,'(A,ES10.2)') 'Mortar coupling: worst relative area mismatch ', worstErr
    call displayGUIMessage( trim(prog_str) )

    deallocate(aList, bList, parentArea, coveredArea, mpNode)

    end subroutine buildMortarInterfaces


    !>-----------------------------------------
    !> @brief
    !> The two in-plane coordinates of the three nodes of a boundary triangle
    !---------------------------------------------------------------------------
    subroutine triOf( nodes, ThoseI, iu, iv, tri )
    real(dp), intent(in) :: nodes(:,:)
    integer, intent(in) :: ThoseI(3), iu, iv
    real(dp), intent(out) :: tri(2,3)

    integer :: i

    do i = 1, 3
        tri(1,i) = nodes(iu,ThoseI(i))
        tri(2,i) = nodes(iv,ThoseI(i))
    end do

    end subroutine triOf


    !>-----------------------------------------
    !> @brief
    !> The signed area of a triangle in two dimensions, positive when it is counterclockwise
    !---------------------------------------------------------------------------
    real(dp) function triArea2D( tri ) result( a )
    real(dp), intent(in) :: tri(2,3)

    a = 0.5_dp * ( (tri(1,2)-tri(1,1))*(tri(2,3)-tri(2,1)) - (tri(1,3)-tri(1,1))*(tri(2,2)-tri(2,1)) )

    end function triArea2D


    !>-----------------------------------------
    !> @brief
    !> Whether the cell (gx,gy) is the one holding the lower left corner of the overlap of the
    !> bounding boxes of two triangles. A triangle is registered in every cell its bounding box
    !> touches, so a pair can be reached from several cells; letting only one of them produce the
    !> sub-face is what keeps each pair counted exactly once. That corner lies inside both bounding
    !> boxes whenever they overlap at all, so its cell is always one of the cells visited.
    !---------------------------------------------------------------------------
    logical function ownsOverlap( triA, triB, gLo2, gridH2, nGrid2, gx, gy ) result( owns )
    real(dp), intent(in) :: triA(2,3), triB(2,3), gLo2(2), gridH2(2)
    integer, intent(in) :: nGrid2(2), gx, gy

    real(dp) :: lo(2)
    integer :: gi2(2), i

    do i = 1, 2
        lo(i) = max( minval(triA(i,:)), minval(triB(i,:)) )
        gi2(i) = min( max( int(floor((lo(i) - gLo2(i)) / gridH2(i))), 0 ), nGrid2(i)-1 )
    end do

    owns = ( gi2(1) .eq. gx ) .and. ( gi2(2) .eq. gy )

    end function ownsOverlap


    !>-----------------------------------------
    !> @brief
    !> The range of cells a two dimensional bounding box touches, clamped to the grid
    !---------------------------------------------------------------------------
    subroutine cellRange2D( bb0, bb1, gLo2, gridH2, nGrid2, gj0, gj1 )
    real(dp), intent(in) :: bb0(2), bb1(2), gLo2(2), gridH2(2)
    integer, intent(in) :: nGrid2(2)
    integer, intent(out) :: gj0(2), gj1(2)

    integer :: i

    do i = 1, 2
        gj0(i) = min( max( int(floor((bb0(i) - gLo2(i)) / gridH2(i))), 0 ), nGrid2(i)-1 )
        gj1(i) = min( max( int(floor((bb1(i) - gLo2(i)) / gridH2(i))), 0 ), nGrid2(i)-1 )
    end do

    end subroutine cellRange2D


    !>-----------------------------------------
    !> @brief
    !> Clips triangle B against triangle A by Sutherland-Hodgman, which is exact here because the
    !> clip polygon is convex. The result is the convex intersection polygon, of at most six
    !> vertices, or fewer than three vertices when the two triangles do not overlap.
    !> @param[in] triA The clip triangle
    !> @param[in] triB The subject triangle
    !> @param[in] tolLen How far outside an edge a vertex may sit and still count as on it
    !> @param[out] poly The intersection polygon
    !> @param[out] npoly Its number of vertices
    !---------------------------------------------------------------------------
    subroutine clipTriangles( triA, triB, tolLen, poly, npoly )
    real(dp), intent(in) :: triA(2,3), triB(2,3)
    real(dp), intent(in) :: tolLen
    real(dp), intent(out) :: poly(2,12)
    integer, intent(out) :: npoly

    real(dp) :: clip(2,3), work(2,12)
    real(dp) :: e0(2), e1(2), ed(2), edLen, dCur, dPrev, tt
    real(dp) :: cur(2), prev(2)
    integer :: i, ie, nwork, ip

    !Sutherland-Hodgman needs the clip polygon counterclockwise, so that 'inside' is to the left
    !of every edge
    clip = triA
    if ( triArea2D(triA) .lt. 0.0_dp ) then
        clip(:,2) = triA(:,3)
        clip(:,3) = triA(:,2)
    endif

    npoly = 3
    poly(:,1:3) = triB

    do ie = 1, 3
        e0 = clip(:,ie)
        e1 = clip(:,mod(ie,3)+1)
        ed = e1 - e0
        edLen = sqrt( ed(1)**2 + ed(2)**2 )
        if ( edLen .le. 0.0_dp ) cycle

        nwork = 0
        do i = 1, npoly
            ip = i - 1
            if ( ip .eq. 0 ) ip = npoly

            cur  = poly(:,i)
            prev = poly(:,ip)

            !Signed distance to the edge, positive on the inside
            dCur  = ( ed(1)*(cur(2)  - e0(2)) - ed(2)*(cur(1)  - e0(1)) ) / edLen
            dPrev = ( ed(1)*(prev(2) - e0(2)) - ed(2)*(prev(1) - e0(1)) ) / edLen

            if ( dCur .ge. -tolLen ) then
                if ( dPrev .lt. -tolLen ) then
                    !dPrev < -tolLen <= dCur, so the denominator is strictly negative
                    tt = dPrev / (dPrev - dCur)
                    nwork = nwork + 1
                    work(:,nwork) = prev + tt*(cur - prev)
                endif
                nwork = nwork + 1
                work(:,nwork) = cur
            else
                if ( dPrev .ge. -tolLen ) then
                    tt = dPrev / (dPrev - dCur)
                    nwork = nwork + 1
                    work(:,nwork) = prev + tt*(cur - prev)
                endif
            endif

            !Clipping a convex polygon of at most six vertices by one line can add at most one
            !vertex, so this is only a guard against an unforeseen degeneracy
            if ( nwork .ge. 11 ) exit
        end do

        npoly = nwork
        if ( npoly .gt. 0 ) poly(:,1:npoly) = work(:,1:npoly)
        if ( npoly .lt. 3 ) exit
    end do

    end subroutine clipTriangles


    !>-----------------------------------------
    !> @brief
    !> The area and the area weighted centroid of a simple polygon, by the shoelace formula. The
    !> centroid, and not the average of the vertices, is what the gradient reconstruction
    !> downstream needs.
    !---------------------------------------------------------------------------
    subroutine polyAreaCentroid( poly, n, area, cen )
    real(dp), intent(in) :: poly(2,12)
    integer, intent(in) :: n
    real(dp), intent(out) :: area, cen(2)

    real(dp) :: a2, cr
    integer :: i, j

    a2 = 0.0_dp
    cen = 0.0_dp
    do i = 1, n
        j = mod(i,n) + 1
        cr = poly(1,i)*poly(2,j) - poly(1,j)*poly(2,i)
        a2 = a2 + cr
        cen(1) = cen(1) + ( poly(1,i) + poly(1,j) ) * cr
        cen(2) = cen(2) + ( poly(2,i) + poly(2,j) ) * cr
    end do
    a2 = 0.5_dp * a2

    if ( abs(a2) .gt. 0.0_dp ) then
        cen = cen / ( 6.0_dp * a2 )
    else
        !A degenerate polygon has no centroid of its own, and the caller discards it on the area
        !anyway
        cen(1) = sum( poly(1,1:n) ) / real(n,dp)
        cen(2) = sum( poly(2,1:n) ) / real(n,dp)
    endif

    area = abs(a2)

    end subroutine polyAreaCentroid


    !>-----------------------------------------
    !> @brief
    !> The zero based cell index of a point in the lookup grid, clamped to the grid
    !---------------------------------------------------------------------------
    subroutine cellOf( x, gLo, gridH, nGrid, gi )
    real(dp), intent(in) :: x(3), gLo(3), gridH(3)
    integer, intent(in) :: nGrid(3)
    integer, intent(out) :: gi(3)

    integer :: i

    do i = 1, 3
        gi(i) = min( max( int(floor((x(i) - gLo(i)) / gridH(i))), 0 ), nGrid(i)-1 )
    end do

    end subroutine cellOf


    !>-----------------------------------------
    !> @brief
    !> The representative of the class of node i, with path compression
    !---------------------------------------------------------------------------
    integer function ufFind( parent, i ) result( r )
    integer, intent(inout) :: parent(:)
    integer, intent(in) :: i

    integer :: cur, nxt

    r = i
    do while ( parent(r) .ne. r )
        r = parent(r)
    end do

    cur = i
    do while ( parent(cur) .ne. r )
        nxt = parent(cur)
        parent(cur) = r
        cur = nxt
    end do

    end function ufFind


    !>-----------------------------------------
    !> @brief
    !> Merges the classes of two nodes. The lower index always becomes the representative, so the
    !> canonical numbering comes out the same no matter in which order the pairs are merged.
    !---------------------------------------------------------------------------
    subroutine ufUnion( parent, a, b )
    integer, intent(inout) :: parent(:)
    integer, intent(in) :: a, b

    integer :: ra, rb

    ra = ufFind( parent, a )
    rb = ufFind( parent, b )
    if ( ra .eq. rb ) return

    if ( ra .lt. rb ) then
        parent(rb) = ra
    else
        parent(ra) = rb
    endif

    end subroutine ufUnion


    !>-----------------------------------------
    !> @brief
    !> The cross product of two vectors
    !---------------------------------------------------------------------------
    function crossProduct( a, b ) result( c )
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

    end function crossProduct


    !>-----------------------------------------
    !> @brief
    !> Checks the status value returned by an allocate statement and aborts with a meaningful
    !> message if the allocation failed. Without this an out-of-memory condition terminates the
    !> process (and thus Matlab) without any diagnostic.
    !> @param[in] stat The status value returned by the allocate statement
    !> @param[in] arrayname The name of the array that was being allocated
    !---------------------------------------------------------------------------
    subroutine checkAllocationTet( stat, arrayname )
    integer, intent(in) :: stat
    character(*), intent(in) :: arrayname

    if ( stat .ne. 0 ) then
        call displayGUIMessage( 'MagTense: out of memory in TetrahedralUnstructuredMeshAnalysis' )
        call displayGUIMessage( 'MagTense: failed to allocate '//trim(arrayname) )
        call displayGUIMessage( 'MagTense: reduce the number of mesh elements or free memory and retry' )
        error stop 'TetrahedralUnstructuredMeshAnalysis: allocation failure'
    endif

    end subroutine checkAllocationTet

end module TetrahedralMeshAnalysis
