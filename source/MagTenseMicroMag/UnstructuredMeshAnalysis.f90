module UnstructuredMeshAnalysis
  use MKL_SPBLAS
  use BLAS95
  use MicroMagParameters
  use IO_GENERAL
  
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
    !> generate 6 faces for each of them an connect them (TheSigns matrix)
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
    !> @param[in] XXX
    !> @param[in] XXX
    !> @param[inout] XXX
    !---------------------------------------------------------------------------
    subroutine CartesianUnstructuredMeshAnalysis(pos, dims, GridInfo)
    real(dp), intent(in) :: pos(:,:)
    real(dp), intent(in) :: dims(:,:)
    type(MicroMagGridInfo), intent(out) :: GridInfo

    integer :: Nel, K, idim, ipm, j, n, kb, i, nnz, k_i, indx, n_faces
    integer, allocatable :: TheJ(:,:)
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
    integer, allocatable :: ThePM(:), TheEE(:,:), TheNot(:,:)
    integer, allocatable :: Aindex(:), Bindex(:), indexItContainsTrue(:,:)
    integer, allocatable :: ia(:), ja(:), MutualContainsInd(:)
    real(dp), allocatable :: a(:)
    logical, allocatable :: ItContains(:,:), ItsContained(:,:), thisBoolean(:,:)
    real(dp), allocatable :: dimscopy(:,:)
    logical, allocatable :: SamePosAlongDim(:)
    logical, allocatable :: UAcontainsUB(:), VAcontainsVB(:), UBcontainsUA(:), VBcontainsVA(:)
    logical, allocatable :: AcontainsB(:), BcontainsA(:)
    integer, allocatable :: indxAB(:), indxBA(:)
    integer, allocatable :: firstindx(:), secondindx(:)
    integer, allocatable :: k1Mut_temp(:), k2Mut_temp(:), k1Mut(:), k2Mut(:), k1NonMut(:), k2NonMut(:)
    integer, allocatable :: kRmv(:), kSurv(:), nRmv(:), TheSigns(:,:), TheSigns_temp(:,:)
    logical, allocatable :: iYes(:), mask1D(:), mask2D1(:,:), mask2D2(:,:)
    integer, allocatable :: cols(:)
    real(dp), allocatable ::xVertA(:,:), xVertB(:,:), xVertC(:,:), xVertD(:,:)
    real(dp), allocatable ::xxMinEl(:,:), xxMaxEl(:,:)
    real(dp), allocatable :: TheseMinXXel(:,:), TheseMaxXXel(:,:)
    logical, allocatable :: theseA(:),theseB(:),theseC(:),theseD(:)
    integer, allocatable :: iZero(:), count_temp(:,:), count1D(:)
    integer, allocatable :: theseA_int(:),theseB_int(:),theseC_int(:),theseD_int(:),theseNum(:)
    logical, allocatable :: TheTs(:,:), TheDs(:,:)
    character*(40) :: prog_str

    ! Initialize some variables
    Nel = size(pos, 1)
    Volumes = dims(:,1) * dims(:,2) * dims(:,3)

    ! Rescaling
    DimsScales = minval(dims, dim=1) / 2.0_dp
       
    Xel = (pos(:,1) / DimsScales(1))
    Yel = (pos(:,2) / DimsScales(2))
    Zel = (pos(:,3) / DimsScales(3))
    XXel = reshape([Xel, Yel, Zel], [Nel, 3])
    
    allocate(dimscopy(Nel, 3))
    do i = 1, 3
        dimscopy(:,i) = dims(:,i) / DimsScales(i)
    end do
         
    ! Construct all faces
    n_faces = 6 * Nel
    allocate(fNormX(n_faces), fNormY(n_faces), fNormZ(n_faces), AreaFaces(n_faces), DimsF(n_faces, 3))
    allocate(TheJ(3, 2), ThePM(2), TheEE(3, 3), TheNot(3, 2))
    allocate(Xf(n_faces), Yf(n_faces), Zf(n_faces), XXF(n_faces,3))
    allocate(TheSigns(Nel,n_faces))

    do i=0,5
        do j=1,Nel
            TheSigns(j, j+(i*Nel)) = 1
        end do
    end do
    
    ThePM = [-1, 1]
    TheEE = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
    TheNot(:,1) = [2, 3, 1]
    TheNot(:,2) = [3, 1, 2]
        
    call displayGUIMessage( 'Test 3' )
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

    call displayGUIMessage( 'Test 7' )

    XXF(:,1) = Xf
    XXF(:,2) = Yf
    XXF(:,3) = Zf
    
    call displayGUIMessage( 'Test 4' )

    ! Check which faces are contained by other faces
    allocate(indexItContainsTrue(n_faces, 2))

    allocate(UminA(Nel), UmaxA(Nel), VminA(Nel), VmaxA(Nel))
    allocate(UminB(Nel), UmaxB(Nel), VminB(Nel), VmaxB(Nel))
    allocate(Aindex(Nel), Bindex(Nel))
    allocate(SamePosAlongDim(Nel))
    allocate(UAcontainsUB(Nel), VAcontainsVB(Nel), UBcontainsUA(Nel), VBcontainsVA(Nel))
    allocate(AcontainsB(Nel), BcontainsA(Nel))
    
    call displayGUIMessage( 'Test 5' )
    
    k = 1;
    do idim = 1, 3
      Aindex = [(1 + (TheJ(idim, 1) - 1) * Nel):(TheJ(idim, 1) * Nel)]
      Bindex = [(1 + (TheJ(idim, 2) - 1) * Nel):(TheJ(idim, 2) * Nel)]
      UminA = XXf(Aindex,TheNot(idim,1)) - 0.5_dp * DimsF(Aindex,TheNot(idim,1))
      UmaxA = XXf(Aindex,TheNot(idim,1)) + 0.5_dp * DimsF(Aindex,TheNot(idim,1))
      VminA = XXf(Aindex,TheNot(idim,2)) - 0.5_dp * DimsF(Aindex,TheNot(idim,2))
      VmaxA = XXf(Aindex,TheNot(idim,2)) + 0.5_dp * DimsF(Aindex,TheNot(idim,2))
    
      call displayGUIMessage( 'Test 6' )
      
      UminB = XXf(Bindex,TheNot(idim,1)) - 0.5_dp * DimsF(Bindex,TheNot(idim,1))
      UmaxB = XXf(Bindex,TheNot(idim,1)) + 0.5_dp * DimsF(Bindex,TheNot(idim,1))
      VminB = XXf(Bindex,TheNot(idim,2)) - 0.5_dp * DimsF(Bindex,TheNot(idim,2))
      VmaxB = XXf(Bindex,TheNot(idim,2)) + 0.5_dp * DimsF(Bindex,TheNot(idim,2))
      
      !Note to self: Code it such that indexItContainsTrue is defined to be an array that is (2*Nel, 2)
      !to ensure that there is enough room for the elements. Then just assign them sequentially, and afterwards make a copy of the array
      !to ensure that it is the right size.
      
      do kb = 1, Nel
        !SamePosAlongDim = (XXf(Aindex,idim) == XXf(Bindex(kb),idim)) ;
        SamePosAlongDim = (abs(XXf(Aindex,idim) - XXf(Bindex(kb),idim)) < 1e-9) ;
        !UAcontainsUB = ((UminA <= UminB(kb)) .and. (UmaxA >= UmaxB(kb)))
        !VAcontainsVB = ((VminA <= VminB(kb)) .and. (VmaxA >= VmaxB(kb)))
        UAcontainsUB = (((UminA - UminB(kb)) < +1e-9) .and. ((UmaxA - UmaxB(kb)) > -1e-9))
        VAcontainsVB = (((VminA - VminB(kb)) < +1e-9) .and. ((VmaxA - VmaxB(kb)) > -1e-9))
        
        !UBcontainsUA = ((UminB(kb) <= UminA) .and. (UmaxB(kb) >= UmaxA))
        !VBcontainsVA = ((VminB(kb) <= VminA) .and. (VmaxB(kb) >= VmaxA))
        UBcontainsUA = (((UminB(kb) - UminA) < +1e-9) .and. ((UmaxB(kb) - UmaxA) > -1e-9))
        VBcontainsVA = (((VminB(kb) - VminA) < +1e-9) .and. ((VmaxB(kb) - VmaxA) > -1e-9))
        
        AcontainsB = (SamePosAlongDim .and. UAcontainsUB .and. VAcontainsVB)
        BcontainsA = (SamePosAlongDim .and. UBcontainsUA .and. VBcontainsVA)
        
        indxAB = pack(Aindex, AcontainsB)
        indxBA = pack(Aindex, BcontainsA)
        
        allocate(firstindx(size(indxAB)+size(indxBA)))
        allocate(secondindx(size(indxAB)+size(indxBA)))
        
        k_i = 1
        do i = 1,size(indxAB)
            firstindx(i) = indxAB(i)
            k_i = k_i+1
        enddo
        do i = k_i,k_i+(size(indxBA)-1)
            firstindx(i) = Bindex(kb)
        enddo
        
        k_i = 1
        do i = 1,size(indxAB)
            secondindx(i) = Bindex(kb)
            k_i = k_i+1
        enddo
        do i = k_i,k_i+(size(indxBA)-1)
            secondindx(i) = indxBA(i-(k_i-1))
        enddo
        
        do i = 1,size(firstindx)
            indexItContainsTrue(k,1) = firstindx(i)
            indexItContainsTrue(k,2) = secondindx(i)
            k = k+1;
        enddo
                
        deallocate(firstindx, secondindx)
      end do               
    end do
    
    
    
    open(21,file='indexItContainsTrue1.txt',status='unknown',form='formatted',action='write')
    do i=1,n_faces
        write(21,*)  indexItContainsTrue(i,1)
    enddo
    close(21)
    
    open(21,file='indexItContainsTrue2.txt',status='unknown',form='formatted',action='write')
    do i=1,n_faces
        write(21,*)  indexItContainsTrue(i,2)
    enddo
    close(21)
    
    
    
    
    
    
    
    
    
    indx = findloc(indexItContainsTrue(:,1), VALUE = 0, DIM = 1)
    allocate(ItContains(n_faces, n_faces), ItsContained(n_faces, n_faces), thisBoolean(n_faces, n_faces))
    ItContains = .false.
    
    call displayGUIMessage( 'Test 109' )
    
    do i = 1,(indx-1)
        ItContains(indexItContainsTrue(i,1),indexItContainsTrue(i,2))  = .true.
    end do
    
    
    
    !open(21,file='ItContains.txt',status='unknown',form='formatted',action='write')
    !do i=1,n_faces 
    !    do j=1,n_faces  
    !        write(21,*)  ItContains(j,i)
    !    end do
    !end do
    !close(21)
    
    ItsContained = transpose(ItContains)
        
    thisBoolean = ItContains .and. ItsContained

    k = count(thisBoolean)
    allocate(k1Mut_temp(k))
    allocate(k2Mut_temp(k))
    k1Mut_temp(:) = 0
    k2Mut_temp(:) = 0
    
    k = 1;
    do i = 1,n_faces
        do j = 1,n_faces
            if (thisBoolean(j,i)) then
                k1Mut_temp(k) = j 
                k2Mut_temp(k) = i
                k = k+1
            end if
        end do
    end do
    !MutualContainsInd = pack(thisBoolean, .true.);
    call displayGUIMessage( 'Test 114' )
    
    allocate(iYes(k-1))
    
    iYes = k1Mut_temp > k2Mut_temp
    
    allocate(k1Mut(count(iYes)))
    allocate(k2Mut(count(iYes)))
    
    call displayGUIMessage( 'Test 115' )

    k1Mut = pack(k1Mut_temp,iYes)  
    k2Mut = pack(k2Mut_temp,iYes) ! rermove only one of the two
    
    call displayGUIMessage( 'Test 116' )
    
    thisBoolean(:,:) = ItContains(:,:)
    !jkC = ItContains == .true.
    do i = 1, n_faces
        do j = 1, n_faces
            if (ItContains(i, j)) then
                thisBoolean(i, j) = .not. (ItsContained(i, j))
            end if
        end do
    end do

    !open(21,file='thisBoolean_first.txt',status='unknown',form='formatted',action='write')
    !do i=1,n_faces 
    !    do j=1,n_faces  
    !        write(21,*)  thisBoolean(j,i)
    !    end do
    !end do
    !close(21)
    
    call displayGUIMessage( 'Test 117' )
    
    allocate(k1NonMut(COUNT(thisBoolean)))
    allocate(k2NonMut(COUNT(thisBoolean)))
    k1NonMut(:) = 0
    k2NonMut(:) = 0
    k = 1
    do i = 1,n_faces
        do j = 1,n_faces
            if (thisBoolean(j,i)) then
                k1NonMut(k) = j 
                k2NonMut(k) = i
                k = k+1
            end if
        end do
    end do
    
    call displayGUIMessage( 'Test 118' )
    
    allocate(kRmv(size(k1Mut)+size(k1NonMut)))
    allocate(kSurv(size(k2Mut)+size(k2NonMut)))
    
    kRmv = [k1Mut, k1NonMut]
    kSurv = [k2Mut, k2NonMut]
    nRmv = MOD(kRmv-1,Nel)+1
    
    call displayGUIMessage( 'Test 119' )

    TheSigns(nRmv,kSurv) = -1

    allocate(mask1D(n_faces))
    allocate(mask2D1(n_faces,3))
    allocate(mask2D2(Nel,n_faces))
    mask1D = .true.
    mask2D1 = .true.
    mask2D2 = .true.
    
    call displayGUIMessage( 'Test 120' )
    
    do i=1,size(kRmv)
        mask1D(kRmv(i)) = .false.
    end do
    
    !As there are duplicate entries in kRmv, this is the way to know the unique number of elements
    k = count(mask1D)   !n_faces-size(kRmv)
    
    write (prog_str,'(I10)') (k)
    call displayGUIMessage( prog_str )
    
    allocate(XXF_temp(k,3),DimsF_temp(k,3))
    allocate(fNormX_temp(k),fNormY_temp(k),fNormZ_temp(k))
    allocate(Xf_temp(k),Yf_temp(k),Zf_temp(k),AreaFaces_temp(k))
    allocate(TheSigns_temp(Nel,k))
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
    do i=1,Nel
        TheSigns_temp(i,:) = pack(TheSigns(i,:),mask1D)
    end do
    !TheSigns(:,unique(kRmv)) = [] ;
    
    call displayGUIMessage( 'Test 121' )
        
    !Delete all the temporary arrays and put the values back into the original arrays
    deallocate(XXF,DimsF,TheSigns)
    deallocate(fNormX,fNormY,fNormZ,Xf,Yf,Zf,AreaFaces)
    deallocate(mask1D)
    
    k = size(Xf_temp)
    allocate(XXF(k,3),DimsF(k,3))
    allocate(fNormX(k),fNormY(k),fNormZ(k))
    allocate(Xf(k),Yf(k),Zf(k),AreaFaces(k))
    allocate(TheSigns(Nel,k))
    allocate(DimsF2(k,3))
    call displayGUIMessage( 'Test 122' )
    XXF(:,:)    = XXF_temp(:,:)
    DimsF(:,:)  = DimsF_temp(:,:)
    call displayGUIMessage( 'Test 123' )
    fNormX(:) = fNormX_temp(:)
    fNormY(:) = fNormY_temp(:)
    fNormZ(:) = fNormZ_temp(:)
    call displayGUIMessage( 'Test 124' )
    Xf(:) = Xf_temp(:)
    Yf(:) = Yf_temp(:)
    Zf(:) = Zf_temp(:)
    call displayGUIMessage( 'Test 125' )
    AreaFaces(:) = AreaFaces_temp(:)
    TheSigns(:,:) = TheSigns_temp(:,:)
    call displayGUIMessage( 'Test 126' )
    deallocate(XXF_temp,DimsF_temp,TheSigns_temp)
    deallocate(fNormX_temp,fNormY_temp,fNormZ_temp,Xf_temp,Yf_temp,Zf_temp,AreaFaces_temp)
    
    call displayGUIMessage( 'Test 127' )
    
    
    K = size(Xf)
    
    allocate(xVertA(K,3), xVertB(K,3), xVertC(K,3), xVertD(K,3))
    allocate(xxMinEl(Nel,3), xxMaxEl(Nel,3))
    allocate(iZero(K), count_temp(K,3))
    
    
    call displayGUIMessage( 'Test 128' )
    
    xVertA = XXf(:,:)+DimsF(:,:)/2.0
    xVertC = XXf(:,:)-DimsF(:,:)/2.0
    DimsF2(:,:) = DimsF(:,:)
    
    call displayGUIMessage( 'Test 129' )
    
    count_temp(:,1) = 1;
    count_temp(:,2) = 2;
    count_temp(:,3) = 3;
    
    allocate(mask1D(K))
    mask1D = .false.
    
    call displayGUIMessage( 'Test 130' )
    
    iZero(:) = 0
    do i = 1,3
        mask1D = (abs(DimsF(:,i)) < 1e-9)
                
        where (mask1D) iZero = iZero + count_temp(:,i)
        !iZero(mask1D) = iZero(mask1D)+count_temp(mask1D,i)
    end do
    
    call displayGUIMessage( 'Test 131' )
    
    allocate(cols(size(iZero)))
    cols = mod(iZero+1,3)+1
    
    call displayGUIMessage( 'Test 131-1' )
    
    do i = 1,size(cols) 
        DimsF2(i,cols(i)) = -DimsF2(i,cols(i)) 
    end do
    
    call displayGUIMessage( 'Test 131-2' )
    
    call displayGUIMessage( 'Test 132' )
    
    xVertB(:,:) = XXf(:,:)+DimsF2(:,:)/2.0 
    xVertD(:,:) = XXf(:,:)-DimsF2(:,:)/2.0 

    call displayGUIMessage( 'Test 132-1' )
    
    xxMinEl(:,:) = XXel(:,:) - dimscopy(:,:)/2.0 
    xxMaxEl(:,:) = XXel(:,:) + dimscopy(:,:)/2.0 
    
    call displayGUIMessage( 'Test 133' )
    
    allocate(TheseMinXXel(size(Xf),3), TheseMaxXXel(size(Xf),3))
    allocate(theseA(size(Xf)),theseB(size(Xf)),theseC(size(Xf)),theseD(size(Xf)))
    allocate(theseA_int(size(Xf)),theseB_int(size(Xf)),theseC_int(size(Xf)),theseD_int(size(Xf)))
    allocate(theseNum(Nel))
    allocate(count1D(Nel))
    TheseMinXXel(:,:) = 0
    
    do n=1,Nel
        count1D(n) = n
    end do
    
    allocate(TheTs(size(Xf),Nel), TheDs(size(Xf),Nel))
    TheTs(:,:) = .false.
    TheDs(:,:) = .false.
    
    do n=1,Nel
        TheseMinXXel(:,1) = xxMinEl(n,1)
        TheseMinXXel(:,2) = xxMinEl(n,2)
        TheseMinXXel(:,3) = xxMinEl(n,3)
        TheseMaxXXel(:,1) = xxMaxEl(n,1)
        TheseMaxXXel(:,2) = xxMaxEl(n,2)
        TheseMaxXXel(:,3) = xxMaxEl(n,3)

        !theseA = all((TheseMinXXel <= xVertA ) .and. (TheseMaxXXel >= xVertA ),2)
        !theseB = all((TheseMinXXel <= xVertB ) .and. (TheseMaxXXel >= xVertB ),2)
        !theseC = all((TheseMinXXel <= xVertC ) .and. (TheseMaxXXel >= xVertC ),2)
        !theseD = all((TheseMinXXel <= xVertD ) .and. (TheseMaxXXel >= xVertD ),2)
        theseA = all((((TheseMinXXel - xVertA) < +1e-9 ) .and. ((TheseMaxXXel - xVertA ) > -1e-9)) ,2)
        theseB = all((((TheseMinXXel - xVertB) < +1e-9 ) .and. ((TheseMaxXXel - xVertB ) > -1e-9)) ,2)
        theseC = all((((TheseMinXXel - xVertC) < +1e-9 ) .and. ((TheseMaxXXel - xVertC ) > -1e-9)) ,2)
        theseD = all((((TheseMinXXel - xVertD) < +1e-9 ) .and. ((TheseMaxXXel - xVertD ) > -1e-9)) ,2)
                    
        theseA_int(:) = 0
        theseB_int(:) = 0
        theseC_int(:) = 0
        theseD_int(:) = 0
        
        where (theseA) theseA_int = 1
        where (theseB) theseB_int = 1
        where (theseC) theseC_int = 1
        where (theseD) theseD_int = 1
        
    !    theseA = all((TheseMinXXel <= xVertA ) & (TheseMaxXXel >= xVertA ),2) ;
    !    theseB = all((TheseMinXXel <= xVertB ) & (TheseMaxXXel >= xVertB ),2) ;
    !    theseC = all((TheseMinXXel <= xVertC ) & (TheseMaxXXel >= xVertC ),2) ;
    !    theseD = all((TheseMinXXel <= xVertD ) & (TheseMaxXXel >= xVertD ),2) ;
        theseNum = theseA_int + theseB_int + theseC_int + theseD_int
        
        
        !allocate(theseKs(count(mask1D)))
        !theseKs = pack(count1D,mask1D))
        !deallocate(theseKs)
        
        mask1D = (theseNum >= 2)
        where (mask1D) TheTs(:,n) = .true.
        
        if (n == 600) then
            write (prog_str,'(I10)') (xxMinEl(n,1))
            call displayGUIMessage( prog_str )
            write (prog_str,'(I10)') (TheseMinXXel(1,1))
            call displayGUIMessage( prog_str )
            write (prog_str,'(I10)') (TheseMinXXel(2,1))
            call displayGUIMessage( prog_str )
            write (prog_str,'(I10)') (TheseMinXXel(3,1))
            call displayGUIMessage( prog_str )
            
            open(21,file='TheseMinXXel1.txt',status='unknown',form='formatted',action='write')
                do i=1,Nel
                    write(21,*)  TheseMinXXel(i,1)
                enddo
            close(21)
            open(21,file='TheseMinXXel2.txt',status='unknown',form='formatted',action='write')
                do i=1,Nel
                    write(21,*)  TheseMinXXel(i,2)
                enddo
            close(21)
            open(21,file='TheseMinXXel3.txt',status='unknown',form='formatted',action='write')
                do i=1,Nel
                    write(21,*)  TheseMinXXel(i,3)
                enddo
            close(21)
    
            open(21,file='TheseMaxXXel1.txt',status='unknown',form='formatted',action='write')
                do i=1,Nel
                    write(21,*)  TheseMaxXXel(i,1)
                enddo
            close(21)
            
             open(21,file='theseA.txt',status='unknown',form='formatted',action='write')
                do j=1,size(theseA)
                    write(21,*)  theseA(j)
                enddo
             close(21)
             
             open(21,file='theseB.txt',status='unknown',form='formatted',action='write')
                do j=1,size(theseB)
                    write(21,*)  theseB(j)
                enddo
             close(21)
             
             open(21,file='theseC.txt',status='unknown',form='formatted',action='write')
                do j=1,size(theseC)
                    write(21,*)  theseC(j)
                enddo
             close(21)
             
             open(21,file='theseD.txt',status='unknown',form='formatted',action='write')
                do j=1,size(theseD)
                    write(21,*)  theseD(j)
                enddo
             close(21)
         
             open(21,file='theseNum.txt',status='unknown',form='formatted',action='write')
                do j=1,size(theseNum)
                    write(21,*)  theseNum(j)
                enddo
             close(21)
         
            
             open(21,file='mask1D1.txt',status='unknown',form='formatted',action='write')
                do j=1,size(mask1D)
                    write(21,*)  mask1D(j)
                enddo
             close(21)
        
             open(21,file='TheTs1.txt',status='unknown',form='formatted',action='write')
                do j=1,size(Xf)
                    write(21,*)  TheTs(j,n)
                enddo
             close(21)
        end if
        
        
        mask1D = (theseNum >= 1)
        where (mask1D) TheDs(:,n) = .true.
        
    !    theseNum = sum([theseA,theseB,theseC,theseD],2) ;
    !    theseKs = find(theseNum>=2) ;
    !    indexTheTs{n} = [ theseKs, repmat(n,size(theseKs,1),1) ];

     !   theseLs = find(theseNum>=1) ;
     !   indexTheDs{n} = [ theseLs, repmat(n,size(theseLs,1),1) ];
    end do

    open(21,file='k1Mut_temp.txt',status='unknown',form='formatted',action='write')
    do i=1,size(k1Mut_temp)
        write(21,*)  k1Mut_temp(i)
    enddo
    close(21)
    
    open(21,file='k2Mut_temp.txt',status='unknown',form='formatted',action='write')
    do i=1,size(k2Mut_temp)
        write(21,*)  k2Mut_temp(i)
    enddo
    close(21)
    
    open(21,file='iYes.txt',status='unknown',form='formatted',action='write')
    do i=1,size(iYes)
        write(21,*)  iYes(i)
    enddo
    close(21)
    
    open(21,file='k1Mut.txt',status='unknown',form='formatted',action='write')
    do i=1,size(k1Mut)
        write(21,*)  k1Mut(i)
    enddo
    close(21)
    
    open(21,file='k2Mut.txt',status='unknown',form='formatted',action='write')
    do i=1,size(k2Mut)
        write(21,*)  k2Mut(i)
    enddo
    close(21)
    
    open(21,file='k1NonMut.txt',status='unknown',form='formatted',action='write')
    do i=1,size(k1NonMut)
        write(21,*)  k1NonMut(i)
    enddo
    close(21)
    
    open(21,file='k2NonMut.txt',status='unknown',form='formatted',action='write')
    do i=1,size(k2NonMut)
        write(21,*)  k2NonMut(i)
    enddo
    close(21)
    
    write (prog_str,'(I10)') (size(Xf))
    call displayGUIMessage( prog_str )
            
    !open(21,file='Xf.txt',status='unknown',form='formatted',action='write')
    !do i=1,size(Xf)
    !    write(21,*)  Xf(i)
    !enddo
    !close(21)
    
    call displayGUIMessage( 'Test 140' )
    
    open(21,file='DimsF1.txt',status='unknown',form='formatted',action='write')
    do i=1,size(DimsF(:,1))
        write(21,*)  DimsF(i,1)
    enddo
    close(21)
    
    open(21,file='DimsF2.txt',status='unknown',form='formatted',action='write')
    do i=1,size(DimsF(:,2))
        write(21,*)  DimsF(i,2)
    enddo
    close(21)
    
    open(21,file='DimsF3.txt',status='unknown',form='formatted',action='write')
    do i=1,size(DimsF(:,3))
        write(21,*)  DimsF(i,3)
    enddo
    close(21)

    call displayGUIMessage( 'Test 140-1' )

    open(21,file='AreaFaces.txt',status='unknown',form='formatted',action='write')
    do i=1,size(AreaFaces)
        write(21,*)  AreaFaces(i)
    enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-2' )
    
    open(21,file='kRmv.txt',status='unknown',form='formatted',action='write')
        do i=1,size(kRmv)
            write(21,*)  kRmv(i)
        enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-3' )
        
    open(21,file='iZero.txt',status='unknown',form='formatted',action='write')
        do i=1,size(iZero)
            write(21,*)  iZero(i)
        enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-7' )
    
    open(21,file='cols.txt',status='unknown',form='formatted',action='write')
        do i=1,size(cols)
            write(21,*)  cols(i)
        enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-8' )
    
    open(21,file='XXel1.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  XXel(i,1)
        enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-9' )
    
    open(21,file='dims1.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  dims(i,1)
        enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-10' )
    
    open(21,file='DimsF21.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  DimsF2(i,1)
        enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-11' )
    
    open(21,file='xxMinEl1.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  xxMinEl(i,1)
        enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-12' )
    
    open(21,file='xxMinEl2.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  xxMinEl(i,2)
        enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-13' )
    
    open(21,file='xxMinEl3.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  xxMinEl(i,3)
        enddo
    close(21)
    
    open(21,file='xxMaxEl1.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  xxMaxEl(i,1)
        enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-12' )
    
    open(21,file='xxMaxEl2.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  xxMaxEl(i,2)
        enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-13' )
    
    open(21,file='xxMaxEl3.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  xxMaxEl(i,3)
        enddo
    close(21)
    
    call displayGUIMessage( 'Test 140-14' )
    
    open(21,file='xVertA1.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertA(i,1)
        enddo
    close(21)
    
    open(21,file='xVertA2.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertA(i,2)
        enddo
    close(21)
    
    open(21,file='xVertA3.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertA(i,3)
        enddo
    close(21)
    
    open(21,file='xVertB1.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertB(i,1)
        enddo
    close(21)
    
    open(21,file='xVertB2.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertB(i,2)
        enddo
    close(21)
    
    open(21,file='xVertB3.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertB(i,3)
        enddo
    close(21)
    
    open(21,file='xVertC1.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertC(i,1)
        enddo
    close(21)
    
    open(21,file='xVertC2.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertC(i,2)
        enddo
    close(21)
    
    open(21,file='xVertC3.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertC(i,3)
        enddo
    close(21)
    
    open(21,file='xVertD1.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertD(i,1)
        enddo
    close(21)
    
    open(21,file='xVertD2.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertD(i,2)
        enddo
    close(21)
    
    open(21,file='xVertD3.txt',status='unknown',form='formatted',action='write')
        do i=1,size(Xf)
            write(21,*)  xVertD(i,3)
        enddo
    close(21)
    call displayGUIMessage( 'Test 140-15' )
    
    open(21,file='TheTs.txt',status='unknown',form='formatted',action='write')
    do i=1,Nel 
        do j=1,size(Xf)        
            write(21,*)  TheTs(j,i)
        end do
    end do
    close(21)
    
    open(21,file='TheDs.txt',status='unknown',form='formatted',action='write')
    do i=1,Nel 
        do j=1,size(Xf)        
            write(21,*)  TheDs(j,i)
        end do
    end do
    close(21)
    
    call displayGUIMessage( 'Test 141' )
    
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

    ! Use Intel MKL to create sparse matrices (TheSigns, TheTs, TheDs)
    ! Further implementation required to create and fill these sparse matrices using MKL
  end subroutine CartesianUnstructuredMeshAnalysis

end module UnstructuredMeshAnalysis