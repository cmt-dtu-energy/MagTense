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

    integer :: Nel, K, idim, ipm, j, n, kb, i, k_i, indx, n_faces, i_end, k1, k2
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
    integer, allocatable :: k1Mut_temp(:), k2Mut_temp(:), k1Mut(:), k2Mut(:), k1NonMut(:), k2NonMut(:)
    integer, allocatable :: kRmv(:), kSurv(:), nRmv(:), TheSigns(:,:), TheSigns_temp(:,:)
    logical, allocatable :: iYes(:), mask1D(:)
    integer, allocatable :: cols(:)
    real(dp), allocatable ::xVertA(:,:), xVertB(:,:), xVertC(:,:), xVertD(:,:)
    real(dp), allocatable ::xxMinEl(:,:), xxMaxEl(:,:)
    real(dp), allocatable :: TheseMinXXel(:,:), TheseMaxXXel(:,:)
    logical, allocatable :: theseA(:),theseB(:),theseC(:),theseD(:)
    integer, allocatable :: iZero(:), count_temp(:,:), count1D(:)
    integer, allocatable :: theseA_int(:),theseB_int(:),theseC_int(:),theseD_int(:),theseNum(:)
    logical, allocatable :: TheTs(:,:), TheDs(:,:)
    integer, allocatable :: thisBoolean_arr(:,:), kMut_F(:,:), kNonMut_F(:,:), kMut_F_temp(:,:), kNonMut_F_temp(:,:)
    character*(40) :: prog_str

    call displayGUIMessage( 'Starting mesh analysis' )
    
    ! Initialize some variables
    Nel = size(pos, 1)
    Volumes = dims(:,1) * dims(:,2) * dims(:,3)

    ! Rescaling
    DimsScales = minval(dims, dim=1) / 2.0_dp
    
    allocate(Xel(Nel), Yel(Nel), Zel(Nel))
    
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
    allocate(indexItContainsTrue_temp(n_faces, 2))

    allocate(UminA(Nel), UmaxA(Nel), VminA(Nel), VmaxA(Nel))
    allocate(UminB(Nel), UmaxB(Nel), VminB(Nel), VmaxB(Nel))
    allocate(Aindex(Nel), Bindex(Nel))
    allocate(SamePosAlongDim(Nel))
    allocate(UAcontainsUB(Nel), VAcontainsVB(Nel), UBcontainsUA(Nel), VBcontainsVA(Nel))
    allocate(AcontainsB(Nel), BcontainsA(Nel))
    
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
      
      !Note: indexItContainsTrue is defined to be an array that is (2*Nel, 2) to ensure that there is enough room for the elements. 
      !These are assigned sequentially, and afterwards make a copy of the array to ensure that it is the right size.
      
      do kb = 1, Nel
        SamePosAlongDim = (abs(XXf(Aindex,idim) - XXf(Bindex(kb),idim)) < 1e-9) ;
        UAcontainsUB = (((UminA - UminB(kb)) < +1e-9) .and. ((UmaxA - UmaxB(kb)) > -1e-9))
        VAcontainsVB = (((VminA - VminB(kb)) < +1e-9) .and. ((VmaxA - VmaxB(kb)) > -1e-9))
       
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
            indexItContainsTrue_temp(k,1) = firstindx(i)
            indexItContainsTrue_temp(k,2) = secondindx(i)
            k = k+1;
        enddo
                
        deallocate(firstindx, secondindx)
      end do               
    end do
    
    allocate(indexItContainsTrue(k-1,2))
    indexItContainsTrue(:,:) = indexItContainsTrue_temp(1:(k-1),:)
    deallocate(indexItContainsTrue_temp)
    
    
    allocate(kMut_F_temp(size(indexItContainsTrue(:,1)),2))
    kMut_F_temp(:,:) = 0
    allocate(kNonMut_F_temp(size(indexItContainsTrue(:,1)),2))
    kNonMut_F_temp(:,:) = 0
    
    allocate(mask1D(size(indexItContainsTrue(:,1))))
    mask1D = .false.
    
    k1 = 1
    k2 = 1
    do i = 1,size(indexItContainsTrue(:,1))
        
        mask1D = ((indexItContainsTrue(i,1) .eq. indexItContainsTrue(:,2)) .and. (indexItContainsTrue(i,2) .eq. indexItContainsTrue(:,1)))
        
        if (any(mask1D, DIM = 1)) then
            write (prog_str,'(I10)') (k1)
            call displayGUIMessage( prog_str )
            
            kMut_F_temp(k1,:) = [indexItContainsTrue(i,1), indexItContainsTrue(i,2)]
            k1 = k1+1
        else 
            kNonMut_F_temp(k2,:) = [indexItContainsTrue(i,1), indexItContainsTrue(i,2)]
            k2 = k2+1
        end if
        
    end do
    deallocate(mask1D)
                
    !Allocate temporary arrays to reduce the size of the original arrays
    allocate(kMut_F(k1-1,2))
    kMut_F(:,:) = 0
    kMut_F(1:(k1-1),:) = kMut_F_temp(1:(k1-1),:)
    deallocate(kMut_F_temp)   
    
    allocate(kNonMut_F(k2-1,2))
    kNonMut_F(:,:) = 0
    kNonMut_F(1:(k2-1),:) = kNonMut_F_temp(1:(k2-1),:)
    deallocate(kNonMut_F_temp)   
    
    iYes = kMut_F(:,1) > kMut_F(:,2)
        
    allocate(k1Mut(count(iYes)))
    allocate(k2Mut(count(iYes)))

    k1Mut = pack(kMut_F(:,1),iYes)  
    k2Mut = pack(kMut_F(:,2),iYes)
    
    
    allocate(kRmv(size(k1Mut)+size(kNonMut_F(:,1))))
    allocate(kSurv(size(k2Mut)+size(kNonMut_F(:,2))))
    allocate(nRmv(size(k2Mut)+size(kNonMut_F(:,2))))
    
    ! each of the faces being removed has:
    !   kSurv: one or more surviving contained (or equal) faces
    !   one element (nRmv) having the face as one of its 6 original boundaries
    kRmv  = [k1Mut, kNonMut_F(:,1)]
    kSurv = [k2Mut, kNonMut_F(:,2)]
    nRmv  = MOD(kRmv-1,Nel)+1
    
    TheSigns(nRmv,kSurv) = -1

    allocate(mask1D(n_faces))
    mask1D = .true.
    mask1D(kRmv(:)) = .false.
    
    !As there are duplicate entries in kRmv, this is the way to know the unique number of elements
    k = count(mask1D)
    
    !Allocate temporary arrays to reduce the size of the original arrays
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
        
    !Delete all the temporary arrays and put the values back into the original arrays
    !First deallocate the original arrays
    deallocate(XXF,DimsF,TheSigns)
    deallocate(fNormX,fNormY,fNormZ,Xf,Yf,Zf,AreaFaces)
    deallocate(mask1D)
    
    k = size(Xf_temp)
    allocate(XXF(k,3),DimsF(k,3),DimsF2(k,3))
    allocate(fNormX(k),fNormY(k),fNormZ(k),Xf(k),Yf(k),Zf(k),AreaFaces(k))
    allocate(TheSigns(Nel,k))
    XXF(:,:)    = XXF_temp(:,:)
    DimsF(:,:)  = DimsF_temp(:,:)
    fNormX(:) = fNormX_temp(:)
    fNormY(:) = fNormY_temp(:)
    fNormZ(:) = fNormZ_temp(:)
    Xf(:) = Xf_temp(:)
    Yf(:) = Yf_temp(:)
    Zf(:) = Zf_temp(:)
    AreaFaces(:) = AreaFaces_temp(:)
    TheSigns(:,:) = TheSigns_temp(:,:)
    deallocate(XXF_temp,DimsF_temp,TheSigns_temp)
    deallocate(fNormX_temp,fNormY_temp,fNormZ_temp,Xf_temp,Yf_temp,Zf_temp,AreaFaces_temp)
    
    !! Construct T and D matrix
    ! TheTs is a K times Nel sparse matrix
    ! The entry TheTs(k,n) is 1 if the n-th element shares at least one edge
    ! (i.e. two points) with the k-th face
    ! The entry TheDs(k,n) is 1 if the n-th element shares at least one vertex
    ! with the k-th face
    ! each face has 4 (A,B,C,D) points with coordinates(xp,yp,zp)
    ! a point is on the boundary of an element if
    !    xp in [xel-dim(1)/2,xel+dim(1)/2]
    !    yp in [yel-dim(2)/2,yel+dim(2)/2]
    !    zp in [zel-dim(3)/2,zel+dim(3)/2]
    ! TheTs = sparse(K,Nel) ;
    ! TheDs = sparse(K,Nel) ;
    K = size(Xf)
    allocate(xVertA(K,3), xVertB(K,3), xVertC(K,3), xVertD(K,3))
    allocate(xxMinEl(Nel,3), xxMaxEl(Nel,3))
    allocate(iZero(K), count_temp(K,3))
    
    xVertA = XXf(:,:)+DimsF(:,:)/2.0
    xVertC = XXf(:,:)-DimsF(:,:)/2.0
    DimsF2(:,:) = DimsF(:,:)
    
    count_temp(:,1) = 1;
    count_temp(:,2) = 2;
    count_temp(:,3) = 3;
    
    allocate(mask1D(K))
    mask1D = .false.
    
    iZero(:) = 0
    do i = 1,3
        mask1D = (abs(DimsF(:,i)) < 1e-9)
        where (mask1D) iZero = iZero + count_temp(:,i)
    end do
    
    allocate(cols(size(iZero)))
    cols = mod(iZero+1,3)+1
    
    do i = 1,size(cols) 
        DimsF2(i,cols(i)) = -DimsF2(i,cols(i)) 
    end do
    
    xVertB(:,:) = XXf(:,:)+DimsF2(:,:)/2.0 
    xVertD(:,:) = XXf(:,:)-DimsF2(:,:)/2.0 

    xxMinEl(:,:) = XXel(:,:) - dimscopy(:,:)/2.0 
    xxMaxEl(:,:) = XXel(:,:) + dimscopy(:,:)/2.0 
    
    k = size(Xf)
    allocate(TheseMinXXel(k,3), TheseMaxXXel(k,3))
    allocate(theseA(k),theseB(K),theseC(k),theseD(k),theseA_int(k),theseB_int(k),theseC_int(k),theseD_int(k))
    allocate(theseNum(Nel),count1D(Nel))
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
    
        theseNum = theseA_int + theseB_int + theseC_int + theseD_int
        
        mask1D = (theseNum >= 2)
        where (mask1D) TheTs(:,n) = .true.
        
        mask1D = (theseNum >= 1)
        where (mask1D) TheDs(:,n) = .true.
    end do
    
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

    call displayGUIMessage( 'Saving GridInfo to disk' )
    
    open(21,file='GridInfo_fNormX.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%fNormX)
        write(21,*)  GridInfo%fNormX(i)
    enddo
    close(21)
    
    open(21,file='GridInfo_fNormY.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%fNormY)
        write(21,*)  GridInfo%fNormY(i)
    enddo
    close(21)
    
    open(21,file='GridInfo_fNormZ.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%fNormZ)
        write(21,*)  GridInfo%fNormZ(i)
    enddo
    close(21)
    
    open(21,file='GridInfo_AreaFaces.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%AreaFaces)
        write(21,*)  GridInfo%AreaFaces(i)
    enddo
    close(21)
    
    open(21,file='GridInfo_Volumes.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%Volumes)
        write(21,*)  GridInfo%Volumes(i)
    enddo
    close(21)
    
    open(21,file='GridInfo_Xel.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%Xel)
        write(21,*)  GridInfo%Xel(i)
    enddo
    close(21)
    
    open(21,file='GridInfo_Yel.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%Yel)
        write(21,*)  GridInfo%Yel(i)
    enddo
    close(21)
    
    open(21,file='GridInfo_Zel.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%Zel)
        write(21,*)  GridInfo%Zel(i)
    enddo
    close(21)
    
    open(21,file='GridInfo_Xf.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%Xf)
        write(21,*)  GridInfo%Xf(i)
    enddo
    close(21)
    
    open(21,file='GridInfo_Yf.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%Yf)
        write(21,*)  GridInfo%Yf(i)
    enddo
    close(21)
    
    open(21,file='GridInfo_Zf.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%Zf)
        write(21,*)  GridInfo%Zf(i)
    enddo
    close(21)
    
    open(21,file='GridInfo_DimsF1.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%DimsF(:,1))
        write(21,*)  GridInfo%DimsF(i,1)
    enddo
    close(21)
    
    open(21,file='GridInfo_DimsF2.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%DimsF(:,2))
        write(21,*)  GridInfo%DimsF(i,2)
    enddo
    close(21)
    
    open(21,file='GridInfo_DimsF3.txt',status='unknown',form='formatted',action='write')
    do i=1,size(GridInfo%DimsF(:,3))
        write(21,*)  GridInfo%DimsF(i,3)
    enddo
    close(21)
        
    open(21,file='GridInfo_TheTs.txt',status='unknown',form='formatted',action='write')
    do i=1,Nel 
        do j=1,size(Xf)        
            write(21,*)  TheTs(j,i)
        end do
    end do
    close(21)
    
    open(21,file='GridInfo_TheDs.txt',status='unknown',form='formatted',action='write')
    do i=1,Nel 
        do j=1,size(Xf)        
            write(21,*)  TheDs(j,i)
        end do
    end do
    close(21)
    
    open(21,file='GridInfo_TheSigns.txt',status='unknown',form='formatted',action='write')
    do i=1,size(Xf) 
        do j=1,Nel    
            write(21,*)  TheSigns(j,i)
        end do
    end do
    close(21)
    
    call displayGUIMessage( 'Mesh analysis done' )
    
    stop
    
    ! Use Intel MKL to create sparse matrices (TheSigns, TheTs, TheDs)
    ! Further implementation required to create and fill these sparse matrices using MKL
  end subroutine CartesianUnstructuredMeshAnalysis

end module UnstructuredMeshAnalysis