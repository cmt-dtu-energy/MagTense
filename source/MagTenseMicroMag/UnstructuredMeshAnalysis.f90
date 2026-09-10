module UnstructuredMeshAnalysis
  use MKL_SPBLAS
  use BLAS95
  use MicroMagParameters
  use IO_GENERAL
  use trace_mod
  
  implicit none

    contains

    !>-----------------------------------------
    !> @author Rasmus Bjørk, rabj@dtu.dk, DTU, 2025
    !> Original Matlab implementation by Andrea Roberto Insinga
    !> @brief
    !> CartesianUnstructuredMeshAnalysis analyze a Cartesian Unstructured Mesh
    !> 
    !> GridInfo = CartesianUnstructuredMeshAnalysis(pos, dims)
    !> model is a structure:
    !> pos is an array with positions of all the mesh-elements
    !> dims is an array with dimensions of all the mesh-elements
    !> GridInfo is a structure with different information about the mesh 
    !> (e.g. normal to the faces, faces-elements connectivity matrix, etc.).
    !> It is used to compute the differential operators (and to plot the mesh)
    !> 
    !> start with all the elements centers ("pos") and cell sizes ("dims")
    !> round centers and sizes to minimum
    !> the elements are already the right number
    !> generate 6 faces for each of them and connect them (TheSigns matrix)
    !> there are too many faces
    !> cycle over each pair of faces.
    !> For each pair (A,B) these are the possible relations
    !>       A equals B
    !>       otherwise
    !>          A contains B
    !>          B contains A
    !>          otherwise
    !>               A shares an edge with B
    !>               A shares no edge with B
    !> if A contains {B_i}, then A needs to be removed
    !> the boundaries of the corresponding element now include all the {B_i}
    !> TheSigns is equal to -1 since the normals point inwards
    !>
    !> If periodic boundary conditions are requested through exchPBC, the elements at the two ends
    !> of a periodic direction are linked together, i.e. they are made to share a face and to enter
    !> each others interpolation stencils, exactly as an ordinary pair of neighbouring elements.
    !> The differential operators computed from the mesh therefore need no special treatment of the
    !> periodic boundaries, apart from measuring the distance between a linked pair of elements
    !> through the boundary rather than across the entire domain.
    !> @param[in] XXX
    !> @param[in] XXX
    !> @param[inout] XXX
    !---------------------------------------------------------------------------
    subroutine CartesianUnstructuredMeshAnalysis(pos, dims, GridInfo, exchPBC)
    real(dp), intent(in) :: pos(:,:)
    real(dp), intent(in) :: dims(:,:)
    type(MicroMagGridInfo), intent(out) :: GridInfo
    integer, intent(in) :: exchPBC(3)                  !> Periodic boundary conditions along x, y and z for the exchange coupling

    integer :: Nel, k, idim, ipm, j, n, kb, i, k_i, n_faces, k1, k2, Ncount
    real(dp), allocatable :: Xel(:), Yel(:), Zel(:)
    real(dp), allocatable :: Volumes(:)
    real(dp) :: DimsScales(3)
    real(dp), allocatable :: XXel(:,:)
    real(dp), allocatable :: Xf(:), Yf(:), Zf(:), XXF(:,:)
    real(dp), allocatable :: XXF_temp(:,:), DimsF_temp(:,:), AreaFaces_temp(:)
    real(dp), allocatable :: fNormX_temp(:),fNormY_temp(:),fNormZ_temp(:)
    real(dp), allocatable :: Xf_temp(:),Yf_temp(:),Zf_temp(:)
    real(dp), allocatable :: fNormX(:), fNormY(:), fNormZ(:)
    real(dp), allocatable :: AreaFaces(:)
    real(dp), allocatable :: DimsF(:,:), DimsF2(:,:)
    real(dp), allocatable :: UminA(:), UmaxA(:), VminA(:), VmaxA(:)
    real(dp), allocatable :: UminB(:), UmaxB(:), VminB(:), VmaxB(:)
    integer :: TheJ(3, 2), ThePM(2), TheEE(3, 3), TheNot(3, 2)
    integer, allocatable :: Aindex(:), Bindex(:), indexItContainsTrue(:,:), indexItContainsTrue_temp(:,:)
    real(dp), allocatable :: dimscopy(:,:)
    logical, allocatable :: SamePosAlongDim(:)
    logical, allocatable :: UAcontainsUB(:), VAcontainsVB(:), UBcontainsUA(:), VBcontainsVA(:)
    logical, allocatable :: AcontainsB(:), BcontainsA(:)
    integer, allocatable :: indxAB(:), indxBA(:)
    integer, allocatable :: firstindx(:), secondindx(:)
    integer, allocatable :: k1Mut_temp(:), k2Mut_temp(:), k1Mut(:), k2Mut(:)
    integer, allocatable :: kRmv(:), kSurv(:), nRmv(:)
    logical, allocatable :: iYes(:), mask1D(:)
    integer, allocatable :: cols(:)
    real(dp), allocatable :: xVertA(:,:), xVertB(:,:), xVertC(:,:), xVertD(:,:)
    real(dp), allocatable :: xxMinEl(:,:), xxMaxEl(:,:)
    real(dp), allocatable :: TheseMinXXel(:,:), TheseMaxXXel(:,:)
    logical, allocatable :: theseA(:), theseB(:), theseC(:), theseD(:)
    integer, allocatable :: iZero(:), count_temp(:,:), count1D(:)
    integer, allocatable :: theseA_int(:), theseB_int(:), theseC_int(:), theseD_int(:), theseNum(:)
    integer, allocatable :: kMut_F(:,:), kNonMut_F(:,:), kMut_F_temp(:,:), kNonMut_F_temp(:,:)
    integer, allocatable :: TheTs_indices_this(:), TheTs_temp(:,:), TheTs_indices(:,:)
    integer, allocatable :: TheDs_indices_this(:), TheDs_temp(:,:), TheDs_indices(:,:)
    integer, allocatable :: TheSigns_indices_pos(:,:), TheSigns_indices_neg(:,:), TheSigns_indices(:,:)
    integer, allocatable :: grow_temp(:,:), maskCount(:)
    integer :: Ncap, Ncap_T, Ncap_D, N_T, N_D, Nadd, alloc_stat
    logical :: PBC(3)
    real(dp) :: globMin(3), globMax(3), Lper(3), sVec(3)
    real(dp), allocatable :: shiftList(:,:)
    integer :: nShifts, ishift, ia, ib, ic
    integer, allocatable :: theseNumMax(:)
    !--- A-face lookup grid for the containment search ---
    integer, allocatable :: aCellCnt(:), aCellStart(:), aCellItem(:)
    integer, allocatable :: aStamp(:), candA(:), abList(:), baList(:)
    integer :: nAGrid(3), nAGridTot, aReg, nCandA, nAB, nBA
    integer :: ax1, ax2, ax3, ja0(3), ja1(3), jc, ju, jv, jj, aidx, itmpv, nsh, ish2
    real(dp) :: aH(3), aLo(3), aHi(3), aBox0(3), aBox1(3), cShift(3), dcoord
    logical :: okSame, okAB, okBA
    !--- first-index buckets for the mutual-pair test ---
    integer, allocatable :: pairCnt(:), pairStart(:), pairSecond(:)
    integer :: nPair, pa, pb, jp
    logical :: isMutual
    !--- element lookup grid for the T/D construction ---
    integer :: nFace, nGrid(3), nGridTot, nReg, icell, gi0(3), gi1(3)
    integer :: gx, gy, gz, iCell0, iv, ish, nn, ncand, idx, numMax
    integer :: nHit, capHit, capCand, iT, iD, ipos
    integer, allocatable :: cellCnt(:), cellStart(:), cellItem(:)
    integer, allocatable :: stampGen(:), stampIdx(:)
    integer, allocatable :: candElem(:), candCnt(:,:)
    integer, allocatable :: hitElem(:), hitFace(:), hitNum(:)
    integer, allocatable :: elemFill(:), elemStart(:), hitOrder(:)
    integer, allocatable :: igrow(:), igrow2(:,:)
    real(dp) :: gridH(3), pq(3), vtx(3,4), eMin(3), eMax(3)
    real(dp), parameter :: regPad = 1.0e-6_dp
    character*(40) :: prog_str
    integer, save :: itimer=0
    
    call trace%begin( "CartesianUnstructuredMeshAnalysis", itimer=itimer )
    
    call displayGUIMessage( 'Starting mesh analysis' )

    
    ! Initialize some variables
    Nel = size(pos, 1)
    allocate(Volumes(Nel), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'Volumes' )
    Volumes = dims(:,1) * dims(:,2) * dims(:,3)
    
    ! Rescaling
    DimsScales = minval(dims, dim=1) / 2.0_dp
    
    allocate(Xel(Nel), Yel(Nel), Zel(Nel), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'Xel/Yel/Zel' )

    Xel = anint(pos(:,1) / DimsScales(1))
    Yel = anint(pos(:,2) / DimsScales(2))
    Zel = anint(pos(:,3) / DimsScales(3))
    
   
    
    allocate(XXel(Nel,3), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'XXel' )
    XXel = reshape([Xel, Yel, Zel], [Nel, 3])

    allocate(dimscopy(Nel, 3), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'dimscopy' )
    do i = 1, 3
        dimscopy(:,i) = dims(:,i) / DimsScales(i)
    end do

    ! Periodic boundary conditions for the exchange interaction.
    ! The bounding box of the mesh defines the period along each of the three directions. Note that
    ! everything below is done in the rescaled coordinates used internally in this routine, so the
    ! period is converted back to the original scale when it is stored in GridInfo.
    PBC = ( exchPBC .ne. 0 )
    do i = 1, 3
        globMin(i) = minval( XXel(:,i) - dimscopy(:,i)/2.0 )
        globMax(i) = maxval( XXel(:,i) + dimscopy(:,i)/2.0 )
        Lper(i) = globMax(i) - globMin(i)
    end do

    do i = 1, 3
        ! The domain has to be thicker than two elements along a periodic direction. Otherwise an
        ! element becomes its own neighbour through the periodic image, and the distance between a
        ! linked pair of elements can no longer be identified unambiguously when the differential
        ! operators are computed.
        if ( PBC(i) .and. ( Lper(i) .lt. (2.0 * maxval(dimscopy(:,i)) + 1e-9) ) ) then
            call displayGUIMessage( 'MagTense: too few elements along a periodic direction' )
            call displayGUIMessage( 'MagTense: at least three elements are required for periodic exchange boundary conditions' )
            error stop 'CartesianUnstructuredMeshAnalysis: too few elements along a periodic direction'
        endif
    end do

    ! Construct all faces
    n_faces = 6 * Nel
    allocate(fNormX(n_faces), fNormY(n_faces), fNormZ(n_faces), AreaFaces(n_faces), DimsF(n_faces, 3), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'fNormX/fNormY/fNormZ/AreaFaces/DimsF' )
    allocate(Xf(n_faces), Yf(n_faces), Zf(n_faces), XXF(n_faces,3), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'Xf/Yf/Zf/XXF' )
    
    ThePM = [-1, 1]
    TheEE = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
    TheNot(:,1) = [2, 3, 1]
    TheNot(:,2) = [3, 1, 2]
        
    j = 0
    do idim = 1, 3
      do ipm = 1, 2
        j = j + 1
        TheJ(idim, ipm) = j
        fNormX((1 + (j-1) * Nel):(j * Nel)) = TheEE(idim, 1) * ThePM(ipm)
        fNormY((1 + (j-1) * Nel):(j * Nel)) = TheEE(idim, 2) * ThePM(ipm)
        fNormZ((1 + (j-1) * Nel):(j * Nel)) = TheEE(idim, 3) * ThePM(ipm)
        Xf((1 + (j-1) * Nel):(j * Nel)) = Xel + TheEE(idim, 1) * ThePM(ipm) * dimscopy(:, 1) / 2.0
        Yf((1 + (j-1) * Nel):(j * Nel)) = Yel + TheEE(idim, 2) * ThePM(ipm) * dimscopy(:, 2) / 2.0
        Zf((1 + (j-1) * Nel):(j * Nel)) = Zel + TheEE(idim, 3) * ThePM(ipm) * dimscopy(:, 3) / 2.0
        DimsF((1 + (j-1) * Nel):(j * Nel), :) = reshape([dimscopy(:, 1) * (1 - TheEE(idim, 1)), dimscopy(:, 2) * (1 - TheEE(idim, 2)), dimscopy(:, 3) * (1 - TheEE(idim, 3))], [Nel, 3])
        AreaFaces((1 + (j-1) * Nel):(j * Nel)) = (DimsScales(TheNot(idim, 1)) * dimscopy(:, TheNot(idim, 1))) * (DimsScales(TheNot(idim, 2)) * dimscopy(:, TheNot(idim, 2)))
      end do
    end do

    XXF(:,1) = Xf
    XXF(:,2) = Yf
    XXF(:,3) = Zf
    
    ! Check which faces are contained by other faces
    !The number of contained-face pairs is not known in advance. The theoretical worst case is
    !6*Nel*Nel (every A-face containing every B-face), but for an actual mesh the count is O(Nel):
    !for a uniform 20x20x20 mesh only 45600 of the 384000000 worst-case rows are used.
    !Allocating the worst case therefore costs 48*Nel^2 bytes (~3 GB at Nel = 8000) and, since the
    !expression is evaluated in default integer kind, it also overflows for Nel > 18918.
    !Instead the buffer starts at a modest size and is doubled on demand in the loop below.
    Ncap = 12 * Nel
    allocate(indexItContainsTrue_temp(Ncap, 2), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'indexItContainsTrue_temp' )
    allocate(UminA(Nel), UmaxA(Nel), VminA(Nel), VmaxA(Nel), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'UminA/UmaxA/VminA/VmaxA' )
    allocate(UminB(Nel), UmaxB(Nel), VminB(Nel), VmaxB(Nel), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'UminB/UmaxB/VminB/VmaxB' )
    allocate(Aindex(Nel), Bindex(Nel), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'Aindex/Bindex' )
    allocate(SamePosAlongDim(Nel), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'SamePosAlongDim' )
    allocate(UAcontainsUB(Nel), VAcontainsVB(Nel), UBcontainsUA(Nel), VBcontainsVA(Nel), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'UAcontainsUB/VAcontainsVB/UBcontainsUA/VBcontainsVA' )
    allocate(AcontainsB(Nel), BcontainsA(Nel), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'AcontainsB/BcontainsA' )
    
    !------------------------------------------------------------------------------------------
    ! Face containment.
    !
    ! Comparing every B face against every A face is O(3*Nel^2); measured, it was ~25% of this
    ! routine before the T/D loop was fixed and the largest remaining term afterwards. Both
    ! containment directions require the two (u,v) rectangles to overlap and the two faces to
    ! sit in the same plane, so the A faces are put in a uniform grid once per dimension and
    ! each B face then only examines the A faces sharing one of its cells.
    !
    ! As in the T/D construction, correctness does not depend on the cell size: two boxes that
    ! intersect have a point in common, and the cell holding that point is covered by both, so
    ! an overlapping A face is always among the candidates. The padding is far above rounding
    ! and far below a cell, making the candidate set a superset; the exact tests below then
    ! decide, written exactly as before against the same UminA/UmaxA/VminA/VmaxA arrays.
    !
    ! Output order is unchanged: dimensions outer, B faces ascending, and within each B face
    ! the "A contains B" rows in ascending A order followed by the "B contains A" rows, which
    ! is what the two pack() calls produced.
    !------------------------------------------------------------------------------------------
    allocate(aStamp(Nel), candA(Nel), abList(Nel), baList(Nel), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'aStamp/candA/abList/baList' )

   k = 1;
    do idim = 1, 3
      ! along a given dimension
      ! compare all the A faces with ipm = -1 with all the B faces with ipm=+1
      Aindex = [(1 + (TheJ(idim, 1) - 1) * Nel):(TheJ(idim, 1) * Nel)]
      Bindex = [(1 + (TheJ(idim, 2) - 1) * Nel):(TheJ(idim, 2) * Nel)]
      UminA = XXf(Aindex,TheNot(idim,1)) - 0.5_dp * DimsF(Aindex,TheNot(idim,1))
      UmaxA = XXf(Aindex,TheNot(idim,1)) + 0.5_dp * DimsF(Aindex,TheNot(idim,1))
      VminA = XXf(Aindex,TheNot(idim,2)) - 0.5_dp * DimsF(Aindex,TheNot(idim,2))
      VmaxA = XXf(Aindex,TheNot(idim,2)) + 0.5_dp * DimsF(Aindex,TheNot(idim,2))

      UminB = XXf(Bindex,TheNot(idim,1)) - 0.5_dp * DimsF(Bindex,TheNot(idim,1))
      UmaxB = XXf(Bindex,TheNot(idim,1)) + 0.5_dp * DimsF(Bindex,TheNot(idim,1))
      VminB = XXf(Bindex,TheNot(idim,2)) - 0.5_dp * DimsF(Bindex,TheNot(idim,2))
      VmaxB = XXf(Bindex,TheNot(idim,2)) + 0.5_dp * DimsF(Bindex,TheNot(idim,2))

      !--- grid axes: 1 = the plane normal, 2 and 3 = the two in-plane directions ------------
      ax1 = idim
      ax2 = TheNot(idim,1)
      ax3 = TheNot(idim,2)
      aLo = [globMin(ax1), globMin(ax2), globMin(ax3)]
      aHi = [globMax(ax1), globMax(ax2), globMax(ax3)]
      aH(1) = minval(dimscopy(:,ax1))
      aH(2) = minval(dimscopy(:,ax2))
      aH(3) = minval(dimscopy(:,ax3))
      do i = 1, 3
          if (aH(i) .le. 0.0_dp) aH(i) = 1.0_dp
          if (aHi(i) .gt. aLo(i)) then
              nAGrid(i) = max(1, int((aHi(i) - aLo(i)) / aH(i)))
          else
              nAGrid(i) = 1
          endif
      end do
      do while ( (real(nAGrid(1),dp)*real(nAGrid(2),dp)*real(nAGrid(3),dp)) .gt. &
                 (4.0_dp * real(Nel,dp) + 1024.0_dp) )
          if ( all(nAGrid .le. 1) ) exit
          nAGrid = max(1, nAGrid / 2)
      end do
      do i = 1, 3
          if (aHi(i) .gt. aLo(i)) then
              aH(i) = (aHi(i) - aLo(i)) / real(nAGrid(i), dp)
          else
              aH(i) = 1.0_dp
          endif
      end do
      nAGridTot = nAGrid(1) * nAGrid(2) * nAGrid(3)

      !--- register the A faces ---------------------------------------------------------------
      allocate(aCellCnt(nAGridTot), aCellStart(nAGridTot+1), stat=alloc_stat)
      call checkAllocation( alloc_stat, 'aCellCnt/aCellStart' )
      aCellCnt = 0
      do kb = 1, Nel
          aBox0 = [XXf(Aindex(kb),ax1), UminA(kb), VminA(kb)] - regPad
          aBox1 = [XXf(Aindex(kb),ax1), UmaxA(kb), VmaxA(kb)] + regPad
          do i = 1, 3
              ja0(i) = min(max(int(floor((aBox0(i) - aLo(i)) / aH(i))), 0), nAGrid(i)-1)
              ja1(i) = min(max(int(floor((aBox1(i) - aLo(i)) / aH(i))), 0), nAGrid(i)-1)
          end do
          do jv = ja0(3), ja1(3)
            do ju = ja0(2), ja1(2)
              do jc = ja0(1), ja1(1)
                jj = 1 + jc + nAGrid(1)*(ju + nAGrid(2)*jv)
                aCellCnt(jj) = aCellCnt(jj) + 1
              end do
            end do
          end do
      end do
      aCellStart(1) = 1
      do i = 1, nAGridTot
          aCellStart(i+1) = aCellStart(i) + aCellCnt(i)
      end do
      aReg = aCellStart(nAGridTot+1) - 1
      allocate(aCellItem(max(1,aReg)), stat=alloc_stat)
      call checkAllocation( alloc_stat, 'aCellItem' )
      aCellCnt = 0
      do kb = 1, Nel
          aBox0 = [XXf(Aindex(kb),ax1), UminA(kb), VminA(kb)] - regPad
          aBox1 = [XXf(Aindex(kb),ax1), UmaxA(kb), VmaxA(kb)] + regPad
          do i = 1, 3
              ja0(i) = min(max(int(floor((aBox0(i) - aLo(i)) / aH(i))), 0), nAGrid(i)-1)
              ja1(i) = min(max(int(floor((aBox1(i) - aLo(i)) / aH(i))), 0), nAGrid(i)-1)
          end do
          do jv = ja0(3), ja1(3)
            do ju = ja0(2), ja1(2)
              do jc = ja0(1), ja1(1)
                jj = 1 + jc + nAGrid(1)*(ju + nAGrid(2)*jv)
                aCellItem(aCellStart(jj) + aCellCnt(jj)) = kb
                aCellCnt(jj) = aCellCnt(jj) + 1
              end do
            end do
          end do
      end do

      !The periodic image of a B face along the plane normal. A face at one end of the domain is
      !also at "the same position" as a face at the other end, which is what the Lper term in the
      !original SamePosAlongDim test expressed; here it is one extra query per direction.
      nsh = 1
      cShift(1) = 0.0_dp
      if ( PBC(idim) ) then
          cShift(2) =  Lper(idim)
          cShift(3) = -Lper(idim)
          nsh = 3
      endif

      aStamp = 0

      do kb = 1, Nel
        !--- gather the A faces that could possibly overlap this B face -----------------------
        nCandA = 0
        do ish2 = 1, nsh
            aBox0 = [XXf(Bindex(kb),ax1) + cShift(ish2), UminB(kb), VminB(kb)] - regPad
            aBox1 = [XXf(Bindex(kb),ax1) + cShift(ish2), UmaxB(kb), VmaxB(kb)] + regPad
            do i = 1, 3
                ja0(i) = min(max(int(floor((aBox0(i) - aLo(i)) / aH(i))), 0), nAGrid(i)-1)
                ja1(i) = min(max(int(floor((aBox1(i) - aLo(i)) / aH(i))), 0), nAGrid(i)-1)
            end do
            do jv = ja0(3), ja1(3)
              do ju = ja0(2), ja1(2)
                do jc = ja0(1), ja1(1)
                  jj = 1 + jc + nAGrid(1)*(ju + nAGrid(2)*jv)
                  do i = aCellStart(jj), aCellStart(jj+1)-1
                      aidx = aCellItem(i)
                      if (aStamp(aidx) .ne. kb) then
                          aStamp(aidx) = kb
                          nCandA = nCandA + 1
                          candA(nCandA) = aidx
                      endif
                  end do
                end do
              end do
            end do
        end do

        !--- ascending A order, as pack() over Aindex produced -------------------------------
        do i = 2, nCandA
            itmpv = candA(i)
            jj = i - 1
            do while (jj .ge. 1)
                if (candA(jj) .le. itmpv) exit
                candA(jj+1) = candA(jj)
                jj = jj - 1
            end do
            candA(jj+1) = itmpv
        end do

        !--- exact tests, written as in the original ----------------------------------------
        nAB = 0
        nBA = 0
        do i = 1, nCandA
            aidx = candA(i)
            dcoord = XXf(Aindex(aidx),idim) - XXf(Bindex(kb),idim)

            okSame = (abs(dcoord) < 1e-9)
            if ( PBC(idim) ) then
                okSame = okSame .or. (abs( abs(dcoord) - Lper(idim) ) < 1e-9)
            endif
            if (.not. okSame) cycle

            okAB = (((UminA(aidx) - UminB(kb)) < +1e-9) .and. ((UmaxA(aidx) - UmaxB(kb)) > -1e-9)) &
             .and. (((VminA(aidx) - VminB(kb)) < +1e-9) .and. ((VmaxA(aidx) - VmaxB(kb)) > -1e-9))
            okBA = (((UminB(kb) - UminA(aidx)) < +1e-9) .and. ((UmaxB(kb) - UmaxA(aidx)) > -1e-9)) &
             .and. (((VminB(kb) - VminA(aidx)) < +1e-9) .and. ((VmaxB(kb) - VmaxA(aidx)) > -1e-9))

            if (okAB) then
                nAB = nAB + 1
                abList(nAB) = aidx
            endif
            if (okBA) then
                nBA = nBA + 1
                baList(nBA) = aidx
            endif
        end do

        !Grow the buffer (by doubling) if the entries about to be added do not fit
        if ((k + nAB + nBA - 1) .gt. Ncap) then
            do while ((k + nAB + nBA - 1) .gt. Ncap)
                Ncap = 2 * Ncap
            enddo
            allocate(grow_temp(Ncap,2), stat=alloc_stat)
            call checkAllocation( alloc_stat, 'indexItContainsTrue_temp' )
            grow_temp(1:(k-1),:) = indexItContainsTrue_temp(1:(k-1),:)
            call move_alloc(grow_temp, indexItContainsTrue_temp)
        endif

        do i = 1, nAB
            indexItContainsTrue_temp(k,1) = Aindex(abList(i))
            indexItContainsTrue_temp(k,2) = Bindex(kb)
            k = k+1
        end do
        do i = 1, nBA
            indexItContainsTrue_temp(k,1) = Bindex(kb)
            indexItContainsTrue_temp(k,2) = Aindex(baList(i))
            k = k+1
        end do
      end do

      deallocate(aCellCnt, aCellStart, aCellItem)
    end do

    deallocate(aStamp, candA, abList, baList)
    
    allocate(indexItContainsTrue(k-1,2), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'indexItContainsTrue' )
    indexItContainsTrue(:,:) = indexItContainsTrue_temp(1:(k-1),:)
    deallocate(indexItContainsTrue_temp)


    allocate(kMut_F_temp(size(indexItContainsTrue(:,1)),2), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'kMut_F_temp' )
    kMut_F_temp(:,:) = 0
    allocate(kNonMut_F_temp(size(indexItContainsTrue(:,1)),2), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'kNonMut_F_temp' )
    kNonMut_F_temp(:,:) = 0

    !------------------------------------------------------------------------------------------
    ! A pair (a,b) is mutual when the reversed pair (b,a) also appears in the list. Asking that
    ! by scanning the whole list for every row is O(P^2), and P grows as ~5.8*Nel, which made
    ! this loop ~32% of the routine at Nel = 32768.
    !
    ! Both columns hold face indices in [1, n_faces], so the rows are bucketed once by their
    ! first column and the query becomes "does bucket b contain the second value a?". The
    ! buckets hold well under one entry each on average (P < n_faces), so each lookup is O(1).
    ! Rows are still visited in ascending order, so kMut_F and kNonMut_F come out unchanged.
    !------------------------------------------------------------------------------------------
    nPair = size(indexItContainsTrue,1)
    allocate(pairCnt(n_faces), pairStart(n_faces+1), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'pairCnt/pairStart' )
    pairCnt = 0
    do i = 1, nPair
        pairCnt(indexItContainsTrue(i,1)) = pairCnt(indexItContainsTrue(i,1)) + 1
    end do
    pairStart(1) = 1
    do i = 1, n_faces
        pairStart(i+1) = pairStart(i) + pairCnt(i)
    end do
    allocate(pairSecond(max(1,nPair)), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'pairSecond' )
    pairCnt = 0
    do i = 1, nPair
        pa = indexItContainsTrue(i,1)
        pairSecond(pairStart(pa) + pairCnt(pa)) = indexItContainsTrue(i,2)
        pairCnt(pa) = pairCnt(pa) + 1
    end do

    k1 = 1
    k2 = 1
    do i = 1, nPair
        pa = indexItContainsTrue(i,1)
        pb = indexItContainsTrue(i,2)

        isMutual = .false.
        do jp = pairStart(pb), pairStart(pb+1)-1
            if (pairSecond(jp) .eq. pa) then
                isMutual = .true.
                exit
            endif
        end do

        if (isMutual) then
            kMut_F_temp(k1,:) = [pa, pb]
            k1 = k1+1
        else 
            kNonMut_F_temp(k2,:) = [pa, pb]
            k2 = k2+1
        end if
    end do
    deallocate(pairCnt, pairStart, pairSecond)
                
    !Allocate temporary arrays to reduce the size of the original arrays
    allocate(kMut_F(k1-1,2), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'kMut_F' )
    kMut_F(:,:) = 0
    kMut_F(1:(k1-1),:) = kMut_F_temp(1:(k1-1),:)
    deallocate(kMut_F_temp)

    allocate(kNonMut_F(k2-1,2), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'kNonMut_F' )
    kNonMut_F(:,:) = 0
    kNonMut_F(1:(k2-1),:) = kNonMut_F_temp(1:(k2-1),:)
    deallocate(kNonMut_F_temp)

    allocate(iYes(size(kMut_F,1)), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'iYes' )
    iYes = kMut_F(:,1) > kMut_F(:,2)

    allocate(k1Mut(count(iYes)), k2Mut(count(iYes)), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'k1Mut/k2Mut' )

    k1Mut = pack(kMut_F(:,1),iYes)  
    k2Mut = pack(kMut_F(:,2),iYes)
    
    
    allocate(kRmv(size(k1Mut)+size(kNonMut_F(:,1))), kSurv(size(k2Mut)+size(kNonMut_F(:,2))), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'kRmv/kSurv' )
    allocate(TheSigns_indices_pos(6*Nel,2), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'TheSigns_indices_pos' )
    allocate(TheSigns_indices_neg((size(k2Mut)+size(kNonMut_F(:,2))),2), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'TheSigns_indices_neg' )

    ! each of the faces being removed has:
    !   kSurv: one or more surviving contained (or equal) faces
    !   one element (nRmv) having the face as one of its 6 original boundaries
    kRmv  = [k1Mut, kNonMut_F(:,1)]
    kSurv = [k2Mut, kNonMut_F(:,2)]
    allocate(nRmv(size(kRmv)), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'nRmv' )
    nRmv  = MOD(kRmv-1,Nel)+1

    allocate(mask1D(n_faces), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'mask1D' )
    mask1D = .true.
    mask1D(kRmv(:)) = .false.

    !Ncount below is the number of surviving faces up to and including a given column. Evaluating
    !it as count(mask1D(1:...)) inside the two loops rescans the mask on every iteration, which
    !makes them O(n_faces^2). The running count is instead tabulated once here in a single pass,
    !so that each lookup becomes O(1). maskCount(i) == count(mask1D(1:i)) by construction.
    allocate(maskCount(n_faces), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'maskCount' )
    Ncount = 0
    do i = 1, n_faces
        if (mask1D(i)) Ncount = Ncount + 1
        maskCount(i) = Ncount
    end do

    !To construct TheSigns matrix, we need to remove all the columns indicated by kRmv
    !This is build into the mask1D array. To then filter the columns, we look at the column value
    !for each indices pair. For the negative values, the column value in the sparse matrix is given by kSurv(i)
    !If the mask1D entry for this is column is true (so the column is not removed), then we can add the entry to the TheSigns_indices_neg array
    !However, to account for the removed columns, its value is no longer kSurv(i), but rather the number of spaces it is moved to the "left" in the matrix
    !This is Ncount, which is simply the number of positive entries in the mask1D array up to the column kSurv(i)
    !The same applies for the positive values except here the y-values are 1:j+(i*Nel)
    k_i = 1
    do i=1,size(nRmv)        
        if (mask1D(kSurv(i))) then
            Ncount = maskCount(kSurv(i))
            TheSigns_indices_neg(k_i,1) = nRmv(i)
            TheSigns_indices_neg(k_i,2) = Ncount
            k_i = k_i + 1
        endif
    end do
    k1 = k_i - 1
        
    k_i = 1
    do i=0,5
        do j=1,Nel  
            if (mask1D(j+(i*Nel))) then
                Ncount = maskCount(j+(i*Nel))
                TheSigns_indices_pos(k_i,1) = j
                TheSigns_indices_pos(k_i,2) = Ncount
                k_i = k_i + 1
            endif
        end do
    end do
    k2 = k_i - 1

    deallocate(maskCount)

    !Combine the two arrays into the final TheSigns_indices array
    allocate(TheSigns_indices(k1+k2,3), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'TheSigns_indices' )
    TheSigns_indices(1:k1,1:2) = TheSigns_indices_neg(1:k1,1:2)
    TheSigns_indices(1:k1,3) = -1
    TheSigns_indices(k1+1:k1+k2,1:2) = TheSigns_indices_pos(1:k2,1:2)
    TheSigns_indices(k1+1:k1+k2,3) = 1
    
    
    !As there are duplicate entries in kRmv, this is the way to know the unique number of elements
    k = count(mask1D)
    
    !Allocate temporary arrays to reduce the size of the original arrays
    allocate(XXF_temp(k,3),DimsF_temp(k,3), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'XXF_temp/DimsF_temp' )
    allocate(fNormX_temp(k),fNormY_temp(k),fNormZ_temp(k), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'fNormX_temp/fNormY_temp/fNormZ_temp' )
    allocate(Xf_temp(k),Yf_temp(k),Zf_temp(k),AreaFaces_temp(k), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'Xf_temp/Yf_temp/Zf_temp/AreaFaces_temp' )
    do i=1,3
        XXF_temp(:,i) = pack(XXf(:,i),mask1D)
        DimsF_temp(:,i) = pack(DimsF(:,i),mask1D)
    end do
    fNormX_temp = pack(fNormX,mask1D)
    fNormY_temp = pack(fNormY,mask1D)
    fNormZ_temp = pack(fNormZ,mask1D)
    Xf_temp = pack(Xf,mask1D)
    Yf_temp = pack(Yf,mask1D)
    Zf_temp = pack(Zf,mask1D)
    AreaFaces_temp = pack(AreaFaces,mask1D)
        
    !Delete all the temporary arrays and put the values back into the original arrays
    !First deallocate the original arrays
    deallocate(XXF,DimsF)
    deallocate(fNormX,fNormY,fNormZ,Xf,Yf,Zf,AreaFaces)
    deallocate(mask1D)
    
    k = size(Xf_temp)
    allocate(XXF(k,3),DimsF(k,3),DimsF2(k,3), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'XXF/DimsF/DimsF2' )
    allocate(fNormX(k),fNormY(k),fNormZ(k),Xf(k),Yf(k),Zf(k),AreaFaces(k), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'fNormX/fNormY/fNormZ/Xf/Yf/Zf/AreaFaces' )
    !allocate(TheSigns(Nel,k))
    call move_alloc (XXF_temp,XXF)
    call move_alloc (DimsF_temp,DimsF)
    call move_alloc (fNormX_temp,fNormX)
    call move_alloc (fNormY_temp,fNormY)
    call move_alloc (fNormZ_temp,fNormZ)
    call move_alloc (Xf_temp,Xf)
    call move_alloc (Yf_temp,Yf)
    call move_alloc (Zf_temp,Zf)
    call move_alloc (AreaFaces_temp,AreaFaces)
    
    !! Construct T and D matrix
    ! TheTs is a k times Nel sparse matrix
    ! The entry TheTs(k,n) is 1 if the n-th element shares at least one edge
    ! (i.e. two points) with the k-th face
    ! The entry TheDs(k,n) is 1 if the n-th element shares at least one vertex
    ! with the k-th face
    ! each face has 4 (A,B,C,D) points with coordinates(xp,yp,zp)
    ! a point is on the boundary of an element if
    !    xp in [xel-dim(1)/2,xel+dim(1)/2]
    !    yp in [yel-dim(2)/2,yel+dim(2)/2]
    !    zp in [zel-dim(3)/2,zel+dim(3)/2]
    ! TheTs = sparse(k,Nel) ;
    ! TheDs = sparse(k,Nel) ;
    
    !k = size(Xf)
    allocate(xVertA(k,3), xVertB(k,3), xVertC(k,3), xVertD(k,3), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'xVertA/xVertB/xVertC/xVertD' )
    allocate(xxMinEl(Nel,3), xxMaxEl(Nel,3), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'xxMinEl/xxMaxEl' )
    allocate(iZero(k), count_temp(k,3), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'iZero/count_temp' )
    
    xVertA = XXf(:,:)+DimsF(:,:)/2.0
    xVertC = XXf(:,:)-DimsF(:,:)/2.0
    DimsF2(:,:) = DimsF(:,:)
    
    count_temp(:,1) = 1;
    count_temp(:,2) = 2;
    count_temp(:,3) = 3;
    
    allocate(mask1D(k), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'mask1D' )
    mask1D = .false.

    iZero(:) = 0
    do i = 1,3
        mask1D = (abs(DimsF(:,i)) < 1e-9)
        where (mask1D) iZero = iZero + count_temp(:,i)
    end do
    
    allocate(cols(size(iZero)), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'cols' )
    cols = mod(iZero+1,3)+1
    
    do i = 1,size(cols) 
        DimsF2(i,cols(i)) = -DimsF2(i,cols(i)) 
    end do
    
    xVertB(:,:) = XXf(:,:)+DimsF2(:,:)/2.0 
    xVertD(:,:) = XXf(:,:)-DimsF2(:,:)/2.0 

    xxMinEl(:,:) = XXel(:,:) - dimscopy(:,:)/2.0 
    xxMaxEl(:,:) = XXel(:,:) + dimscopy(:,:)/2.0 
    
    nFace = size(Xf)
    N_T = 0
    N_D = 0

    !The list of periodic images that each element is tested against below. Without periodic
    !boundary conditions this is just the single zero shift, and the loop below is then identical to
    !the non-periodic version. All combinations of the periodic directions are needed, as an element
    !can share a vertex with a face across two (or three) periodic directions at the same time.
    allocate(shiftList(27,3), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'shiftList' )
    nShifts = 0
    do ia = -1, 1
      do ib = -1, 1
        do ic = -1, 1
          if ( (ia .ne. 0 .and. .not. PBC(1)) .or. (ib .ne. 0 .and. .not. PBC(2)) &
                .or. (ic .ne. 0 .and. .not. PBC(3)) ) cycle
          nShifts = nShifts + 1
          shiftList(nShifts,:) = [ia*Lper(1), ib*Lper(2), ic*Lper(3)]
        end do
      end do
    end do

    !------------------------------------------------------------------------------------------
    ! T and D construction.
    !
    ! The direct formulation - for every element, test every face's four vertices - is
    ! O(Nel * nFace). Measured, it was ~60% of this routine and grew as Nel^2.1. The test is a
    ! point-in-box query, so it is inverted here: every element is registered once in a uniform
    ! lookup grid, and each face vertex then examines only the elements sharing its own cell.
    !
    ! Correctness does not depend on the cell size. If an element's padded box contains the
    ! query point p, then that box and the cell holding p have p in common, so the element was
    ! registered in that cell and is always among the candidates. The cell size therefore only
    ! trades memory against candidates per cell, and is picked to keep the grid near Nel cells.
    ! The padding used when registering (regPad) is far larger than any rounding difference
    ! between the two ways of writing the containment test, and far smaller than a cell, so the
    ! candidate set is always a superset of the true one; the exact test below then decides,
    ! written with the same arithmetic as the original so the accept/reject boundary is bit
    ! identical.
    !
    ! Output order is preserved exactly. Faces are visited in ascending order, and the
    ! (element, face) hits are then counting-sorted by element. That sort is stable, so faces
    ! stay ascending inside each element - reproducing the original "for each element n, faces
    ! in ascending order" layout row for row.
    !------------------------------------------------------------------------------------------

    !--- choose the grid ----------------------------------------------------------------------
    do i = 1, 3
        gridH(i) = minval(dimscopy(:,i))
        if (gridH(i) .le. 0.0_dp) gridH(i) = 1.0_dp
        if (globMax(i) .gt. globMin(i)) then
            nGrid(i) = max(1, int((globMax(i) - globMin(i)) / gridH(i)))
        else
            nGrid(i) = 1
        endif
    end do
    do while ( (real(nGrid(1),dp) * real(nGrid(2),dp) * real(nGrid(3),dp)) .gt. &
               (4.0_dp * real(Nel,dp) + 1024.0_dp) )
        if ( all(nGrid .le. 1) ) exit
        nGrid = max(1, nGrid / 2)
    end do
    do i = 1, 3
        if (globMax(i) .gt. globMin(i)) then
            gridH(i) = (globMax(i) - globMin(i)) / real(nGrid(i), dp)
        else
            gridH(i) = 1.0_dp
        endif
    end do
    nGridTot = nGrid(1) * nGrid(2) * nGrid(3)

    !--- register every element in the cells its padded box overlaps --------------------------
    allocate(cellCnt(nGridTot), cellStart(nGridTot+1), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'cellCnt/cellStart' )
    cellCnt = 0
    do n = 1, Nel
        do i = 1, 3
            gi0(i) = min(max(int(floor((xxMinEl(n,i) - regPad - globMin(i)) / gridH(i))), 0), nGrid(i)-1)
            gi1(i) = min(max(int(floor((xxMaxEl(n,i) + regPad - globMin(i)) / gridH(i))), 0), nGrid(i)-1)
        end do
        do gz = gi0(3), gi1(3)
          do gy = gi0(2), gi1(2)
            do gx = gi0(1), gi1(1)
              icell = 1 + gx + nGrid(1)*(gy + nGrid(2)*gz)
              cellCnt(icell) = cellCnt(icell) + 1
            end do
          end do
        end do
    end do

    cellStart(1) = 1
    do i = 1, nGridTot
        cellStart(i+1) = cellStart(i) + cellCnt(i)
    end do
    nReg = cellStart(nGridTot+1) - 1
    allocate(cellItem(max(1,nReg)), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'cellItem' )

    cellCnt = 0
    do n = 1, Nel
        do i = 1, 3
            gi0(i) = min(max(int(floor((xxMinEl(n,i) - regPad - globMin(i)) / gridH(i))), 0), nGrid(i)-1)
            gi1(i) = min(max(int(floor((xxMaxEl(n,i) + regPad - globMin(i)) / gridH(i))), 0), nGrid(i)-1)
        end do
        do gz = gi0(3), gi1(3)
          do gy = gi0(2), gi1(2)
            do gx = gi0(1), gi1(1)
              icell = 1 + gx + nGrid(1)*(gy + nGrid(2)*gz)
              cellItem(cellStart(icell) + cellCnt(icell)) = n
              cellCnt(icell) = cellCnt(icell) + 1
            end do
          end do
        end do
    end do

    !--- collect the (element, face) hits, faces in ascending order ---------------------------
    allocate(stampGen(Nel), stampIdx(Nel), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'stampGen/stampIdx' )
    stampGen = 0

    capCand = 64
    allocate(candElem(capCand), candCnt(capCand, max(1,nShifts)), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'candElem/candCnt' )

    capHit = 64 * Nel
    nHit = 0
    allocate(hitElem(capHit), hitFace(capHit), hitNum(capHit), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'hitElem/hitFace/hitNum' )

    do k = 1, nFace
        vtx(:,1) = xVertA(k,:)
        vtx(:,2) = xVertB(k,:)
        vtx(:,3) = xVertC(k,:)
        vtx(:,4) = xVertD(k,:)
        ncand = 0

        do ish = 1, nShifts
            sVec = shiftList(ish,:)

            do iv = 1, 4
                pq = vtx(:,iv) - sVec

                !A point outside the mesh bounding box cannot lie in any element, since globMin
                !and globMax are the extremes of the element boxes themselves. Skipping it also
                !keeps the clamped cell index below meaningful.
                if ( any(pq .lt. (globMin - 1e-9)) .or. any(pq .gt. (globMax + 1e-9)) ) cycle

                do i = 1, 3
                    gi0(i) = min(max(int(floor((pq(i) - globMin(i)) / gridH(i))), 0), nGrid(i)-1)
                end do
                iCell0 = 1 + gi0(1) + nGrid(1)*(gi0(2) + nGrid(2)*gi0(3))

                do ipos = cellStart(iCell0), cellStart(iCell0+1) - 1
                    nn = cellItem(ipos)

                    eMin = xxMinEl(nn,:) + sVec
                    eMax = xxMaxEl(nn,:) + sVec

                    !Skip the image if the shifted element falls entirely outside the mesh, exactly
                    !as the element-major version did. Only elements at a boundary have an image
                    !that can touch a face at the opposite boundary; the zero shift never skips.
                    if ( any( eMin .gt. (globMax + 1e-9) ) .or. &
                         any( eMax .lt. (globMin - 1e-9) ) ) cycle

                    if ( all( ((eMin - vtx(:,iv)) .lt. +1e-9) .and. &
                              ((eMax - vtx(:,iv)) .gt. -1e-9) ) ) then
                        if (stampGen(nn) .ne. k) then
                            ncand = ncand + 1
                            if (ncand .gt. capCand) then
                                allocate(igrow(2*capCand), stat=alloc_stat)
                                call checkAllocation( alloc_stat, 'candElem' )
                                igrow(1:capCand) = candElem(1:capCand)
                                call move_alloc(igrow, candElem)
                                allocate(igrow2(2*capCand, max(1,nShifts)), stat=alloc_stat)
                                call checkAllocation( alloc_stat, 'candCnt' )
                                igrow2(1:capCand,:) = candCnt(1:capCand,:)
                                call move_alloc(igrow2, candCnt)
                                capCand = 2 * capCand
                            endif
                            stampGen(nn) = k
                            stampIdx(nn) = ncand
                            candElem(ncand) = nn
                            candCnt(ncand,1:nShifts) = 0
                        endif
                        idx = stampIdx(nn)
                        candCnt(idx,ish) = candCnt(idx,ish) + 1
                    endif
                end do
            end do
        end do

        !The maximum, and not the sum, over the images ensures that an element is only entered
        !once in TheTs and TheDs, even if it touches the face both directly and through one of
        !its periodic images
        do idx = 1, ncand
            numMax = maxval(candCnt(idx,1:nShifts))
            if (numMax .ge. 1) then
                nHit = nHit + 1
                if (nHit .gt. capHit) then
                    capHit = 2 * capHit
                    allocate(igrow(capHit), stat=alloc_stat)
                    call checkAllocation( alloc_stat, 'hitElem' )
                    igrow(1:nHit-1) = hitElem(1:nHit-1)
                    call move_alloc(igrow, hitElem)
                    allocate(igrow(capHit), stat=alloc_stat)
                    call checkAllocation( alloc_stat, 'hitFace' )
                    igrow(1:nHit-1) = hitFace(1:nHit-1)
                    call move_alloc(igrow, hitFace)
                    allocate(igrow(capHit), stat=alloc_stat)
                    call checkAllocation( alloc_stat, 'hitNum' )
                    igrow(1:nHit-1) = hitNum(1:nHit-1)
                    call move_alloc(igrow, hitNum)
                endif
                hitElem(nHit) = candElem(idx)
                hitFace(nHit) = k
                hitNum(nHit)  = numMax
            endif
        end do
    end do

    !--- stable counting sort by element, so the rows come out element-major -------------------
    allocate(elemFill(Nel), elemStart(Nel+1), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'elemFill/elemStart' )
    elemFill = 0
    do i = 1, nHit
        elemFill(hitElem(i)) = elemFill(hitElem(i)) + 1
    end do
    elemStart(1) = 1
    do n = 1, Nel
        elemStart(n+1) = elemStart(n) + elemFill(n)
    end do
    allocate(hitOrder(max(1,nHit)), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'hitOrder' )
    elemFill = 0
    do i = 1, nHit
        n = hitElem(i)
        hitOrder(elemStart(n) + elemFill(n)) = i
        elemFill(n) = elemFill(n) + 1
    end do

    !--- emit TheTs and TheDs -----------------------------------------------------------------
    N_D = nHit
    N_T = 0
    do i = 1, nHit
        if (hitNum(i) .ge. 2) N_T = N_T + 1
    end do

    allocate(TheTs_indices(N_T,2), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'TheTs_indices' )
    allocate(TheDs_indices(N_D,2), stat=alloc_stat)
    call checkAllocation( alloc_stat, 'TheDs_indices' )

    iT = 0
    iD = 0
    do j = 1, nHit
        i = hitOrder(j)
        iD = iD + 1
        TheDs_indices(iD,1) = hitFace(i)
        TheDs_indices(iD,2) = hitElem(i)
        if (hitNum(i) .ge. 2) then
            iT = iT + 1
            TheTs_indices(iT,1) = hitFace(i)
            TheTs_indices(iT,2) = hitElem(i)
        endif
    end do

    deallocate(cellCnt, cellStart, cellItem, stampGen, stampIdx)
    deallocate(candElem, candCnt, hitElem, hitFace, hitNum)
    deallocate(elemFill, elemStart, hitOrder)

    !TheTs_indices and TheDs_indices are already allocated at exactly N_T and N_D rows
    !by the construction above, so no trimming pass is needed.

    ! Rescaling (inverse)
    Xel = Xel * DimsScales(1)
    Yel = Yel * DimsScales(2)
    Zel = Zel * DimsScales(3)
    Xf = Xf * DimsScales(1)
    Yf = Yf * DimsScales(2)
    Zf = Zf * DimsScales(3)
    DimsF(:, 1) = DimsF(:, 1) * DimsScales(1)
    DimsF(:, 2) = DimsF(:, 2) * DimsScales(2)
    DimsF(:, 3) = DimsF(:, 3) * DimsScales(3)
    ! areas and volumes are already in the original scale
    ! FaceNormals are already unit-normals

    ! Fill GridInfo structure
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
    GridInfo%TheTs = TheTs_indices
    GridInfo%TheDs = TheDs_indices
    GridInfo%TheSigns = TheSigns_indices

    !The periodic directions and the corresponding periods, converted back to the original scale.
    !These are needed when the differential operators are computed, as the distance between a pair
    !of elements linked across a periodic boundary has to be measured through that boundary.
    GridInfo%exchPBC = PBC
    GridInfo%Lper = Lper * DimsScales


    call displayGUIMessage( 'Mesh analysis done' )

    call trace%end( "CartesianUnstructuredMeshAnalysis", itimer=itimer )
    
  end subroutine CartesianUnstructuredMeshAnalysis

    !>-----------------------------------------
    !> @brief
    !> Checks the status value returned by an allocate statement and aborts with a
    !> meaningful message if the allocation failed. Without this an out-of-memory
    !> condition terminates the process (and thus Matlab) without any diagnostic.
    !> @param[in] stat The status value returned by the allocate statement
    !> @param[in] arrayname The name of the array that was being allocated
    !---------------------------------------------------------------------------
    subroutine checkAllocation( stat, arrayname )
    integer, intent(in) :: stat
    character(*), intent(in) :: arrayname

    if ( stat .ne. 0 ) then
        call displayGUIMessage( 'MagTense: out of memory in CartesianUnstructuredMeshAnalysis' )
        call displayGUIMessage( 'MagTense: failed to allocate '//trim(arrayname) )
        call displayGUIMessage( 'MagTense: reduce the number of mesh elements or free memory and retry' )
        error stop 'CartesianUnstructuredMeshAnalysis: allocation failure'
    endif

    end subroutine checkAllocation

end module UnstructuredMeshAnalysis
