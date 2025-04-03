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
    !TheNot = reshape([2, 3, 3, 1, 1, 2], [3, 2])
    TheNot(:,1) = [2, 3, 1]
    TheNot(:,2) = [3, 1, 2]

    !call displayGUIMessage( 'TheNot' )
    !write (prog_str,'(I10)') (TheNot(1, 1))
    !call displayGUIMessage( prog_str )
    !write (prog_str,'(I10)') (TheNot(2, 1))
    !call displayGUIMessage( prog_str )
    !write (prog_str,'(I10)') (TheNot(3, 1))
    !call displayGUIMessage( prog_str )
    !write (prog_str,'(I10)') (TheNot(1, 2))
    !call displayGUIMessage( prog_str )
    !write (prog_str,'(I10)') (TheNot(2, 2))
    !call displayGUIMessage( prog_str )
    !write (prog_str,'(I10)') (TheNot(3, 2))
    !call displayGUIMessage( prog_str )
            
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
        !call displayGUIMessage( 'AreaFaces' )
        !write (prog_str,'(ES18.8)') AreaFaces(1)
        !call displayGUIMessage( prog_str )
        call displayGUIMessage( 'Test 7' )
      end do
    end do
    XXF(:,1) = Xf
    XXF(:,2) = Yf
    XXF(:,3) = Zf
    
    call displayGUIMessage( 'Test 4' )

    open(21,file='Xf.txt',status='unknown',form='formatted',action='write')
    do i=1,n_faces
        write(21,*)  Xf(i)
    enddo
    close(21)
    
    open(21,file='DimsF.txt',status='unknown',form='formatted',action='write')
    do i=1,n_faces
        write(21,*)  DimsF(i,1)
    enddo
    close(21)
    
    open(21,file='AreaFaces.txt',status='unknown',form='formatted',action='write')
    do i=1,n_faces
        write(21,*)  AreaFaces(i)
    enddo
    close(21)
    
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
        
        !open(21,file='BcontainsA.txt',status='unknown',form='formatted',action='write')
        !    do i=1,Nel
        !        write(21,*)  BcontainsA(i)
        !enddo
        !close(21) 
            
        
        
        indxAB = pack(Aindex, AcontainsB)
        indxBA = pack(Aindex, BcontainsA)
        
        !open(21,file='indxBA.txt',status='unknown',form='formatted',action='write')
        !    do i=1,size(indxBA)
        !        write(21,*)  indxBA(i)
        !enddo
        !close(21) 
            
        
        
        !call displayGUIMessage( 'Test 6.1' )
        
        !write (prog_str,'(I10)') (size(indxAB))
        !call displayGUIMessage( prog_str )
        !write (prog_str,'(I10)') (size(indxBA))
        !call displayGUIMessage( prog_str )
        
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
        
        !open(21,file='firstindx.txt',status='unknown',form='formatted',action='write')
        !    do i=1,size(firstindx)
        !        write(21,*)  firstindx(i)
        !enddo
        !close(21) 
        
        !if (k == 474) then
        !    call displayGUIMessage( 'Test 474' )
            
        !    open(21,file='SamePosAlongDim.txt',status='unknown',form='formatted',action='write')
        !    do i=1,Nel
        !        write(21,*)  SamePosAlongDim(i)
        !    enddo
        !    close(21)
            
        !    open(21,file='UBcontainsUA.txt',status='unknown',form='formatted',action='write')
        !    do i=1,Nel
        !        write(21,*)  UBcontainsUA(i)
        !    enddo
        !    close(21)
            
        !    open(21,file='VBcontainsVA.txt',status='unknown',form='formatted',action='write')
        !    do i=1,Nel
        !        write(21,*)  VBcontainsVA(i)
        !    enddo
        !    close(21)  
            
           
        !    call displayGUIMessage( 'Test 474 end' )
        !endif
        
        !do i = 1,size(firstindx)
        !    write (prog_str,'(I10)') (firstindx(i))
        !    call displayGUIMessage( prog_str )
        !enddo
        
        k_i = 1
        do i = 1,size(indxAB)
            secondindx(i) = Bindex(kb)
            k_i = k_i+1
        enddo
        do i = k_i,k_i+(size(indxBA)-1)
            secondindx(i) = indxBA(i-(k_i-1))
        enddo
        
        !do i = 1,size(secondindx)
        !    write (prog_str,'(I10)') (secondindx(i))
        !    call displayGUIMessage( prog_str )
        !enddo
      
        !call displayGUIMessage( 'Test 6.2' )
            
        !firstindx = [indxAB repmat(Bindex(kb),1,length(indxBA))]
        !secondindx = [repmat(Bindex(kb),1,length(indxAB)) indxBA]

        do i = 1,size(firstindx)
            indexItContainsTrue(k,1) = firstindx(i)
            indexItContainsTrue(k,2) = secondindx(i)
            k = k+1;
        enddo
                
        deallocate(firstindx, secondindx)
      end do               
    end do
    

          open(21,file='indexItContainsTrue1.txt',status='unknown',form='formatted',action='write')
          do i=1,Nel*2*3
              write(21,*)  indexItContainsTrue(i,1)
          enddo
          close(21)
          
          open(21,file='indexItContainsTrue2.txt',status='unknown',form='formatted',action='write')
          do i=1,Nel*2*3
              write(21,*)  indexItContainsTrue(i,2)
          enddo
          close(21)


    indx = findloc(indexItContainsTrue(:,1), VALUE = 0, DIM = 1)
    allocate(ItContains(n_faces, n_faces), ItsContained(n_faces, n_faces), thisBoolean(n_faces, n_faces))
    ItContains = .false.
    
    call displayGUIMessage( 'Test 109' )
    write (prog_str,'(I10)') (indx)
    call displayGUIMessage( prog_str )
    
    call displayGUIMessage( 'Test 110' )
    
    !ItContains(indexItContainsTrue(1:(indx-1),1), indexItContainsTrue(1:(indx-1),2)) = .true.
    do i = 1,(indx-1)
        ItContains(indexItContainsTrue(i,1),indexItContainsTrue(i,2))  = .true.
    end do
    
    !ItContains(indexItContainsTrue(1:(indx-1),1),indexItContainsTrue(1:(indx-1),2)) = .true.
    
    k = COUNT(ItContains)
    write (prog_str,'(I10)') (k)
    call displayGUIMessage( prog_str )
    
    call displayGUIMessage( 'Test 111' )
    ItsContained = transpose(ItContains)
    
    k = COUNT(ItsContained)
    write (prog_str,'(I10)') (k)
    call displayGUIMessage( prog_str )
    
    
    !open(21,file='ItContains.txt',status='unknown',form='formatted',action='write')
    !do i = 1,n_faces
    !    do j = 1,n_faces
    !          write(21,*)  ItContains(i,j)
    !    enddo
    !enddo
    !close(21)
          
    !      open(21,file='ItsContained.txt',status='unknown',form='formatted',action='write')
    !      do i = 1,n_faces
    !    do j = 1,n_faces
    !          write(21,*)  ItsContained(i,2)
    !    enddo
    !    enddo
    !      close(21)
    
    call displayGUIMessage( 'Test 112' )
    thisBoolean = ItContains .and. ItsContained
    call displayGUIMessage( 'Test 113' )
    
    k = COUNT(thisBoolean)
    write (prog_str,'(I10)') (k)
    call displayGUIMessage( prog_str )
    
    allocate(k1Mut_temp(k))
    allocate(k2Mut_temp(k))
    k = 1;
    do i = 1,n_faces
        do j = 1,n_faces
            if (thisBoolean(i,j)) then
                k1Mut_temp(k) = i 
                k2Mut_temp(k) = j
                k = k+1
            end if
        end do
    end do
    !MutualContainsInd = pack(thisBoolean, .true.);
    call displayGUIMessage( 'Test 114' )
    
    allocate(iYes(k-1))
    
    iYes = k1Mut_temp > k2Mut_temp
    
    k = COUNT(iYes)
    write (prog_str,'(I10)') (k)
    call displayGUIMessage( prog_str )
    
    
    allocate(k1Mut(k))
    allocate(k2Mut(k))
    
    call displayGUIMessage( 'Test 115' )

    k1Mut = pack(k1Mut_temp,iYes)  
    k2Mut = pack(k2Mut_temp,iYes) ! rermove only one of the two
    
    call displayGUIMessage( 'Test 116' )
    
    thisBoolean = ItContains
    !jkC = ItContains == .true.
    do i = 1, n_faces
        do j = 1, n_faces
            if (ItContains(i, j)) then
                thisBoolean(i, j) = .not. (ItsContained(i, j))
            end if
        end do
    end do

    call displayGUIMessage( 'Test 117' )
    
    allocate(k1NonMut(COUNT(thisBoolean)))
    allocate(k2NonMut(COUNT(thisBoolean)))
    k = 1
    do i = 1,n_faces
        do j = 1,n_faces
            if (thisBoolean(i,j)) then
                k1NonMut(k) = i 
                k2NonMut(k) = j
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

    open(21,file='kRmv.txt',status='unknown',form='formatted',action='write')
        do i=1,size(kRmv)
            write(21,*)  kRmv(i)
        enddo
    close(21)
    
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
    
    open(21,file='DimsF_temp1.txt',status='unknown',form='formatted',action='write')
        do i=1,size(DimsF_temp(:,1))
            write(21,*)  DimsF_temp(i,1)
        enddo
    close(21)
    
    open(21,file='DimsF_temp2.txt',status='unknown',form='formatted',action='write')
        do i=1,size(DimsF_temp(:,2))
            write(21,*)  DimsF_temp(i,2)
        enddo
    close(21)
    
    open(21,file='DimsF_temp3.txt',status='unknown',form='formatted',action='write')
        do i=1,size(DimsF_temp(:,3))
            write(21,*)  DimsF_temp(i,3)
        enddo
    close(21)
    
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
    call displayGUIMessage( 'Test 125-1' )
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
        
        if (i == 1) then
             open(21,file='mask1D1.txt',status='unknown',form='formatted',action='write')
                do j=1,size(mask1D)
                    write(21,*)  mask1D(j)
                enddo
             close(21)
             
             where (mask1D) iZero = count_temp(:,i)
             
             !iZero(mask1D) = count_temp(mask1D,i)
             
             open(21,file='iZero.txt',status='unknown',form='formatted',action='write')
                do j=1,size(iZero)
                    write(21,*)  iZero(j)
                enddo
             close(21)
             
             iZero(:) = 0
             
        end if
        
        where (mask1D) iZero = iZero + count_temp(:,i)
        !iZero(mask1D) = iZero(mask1D)+count_temp(mask1D,i)
    end do
   
    open(21,file='iZero2.txt',status='unknown',form='formatted',action='write')
        do i=1,size(iZero)
            write(21,*)  iZero(i)
        enddo
    close(21)
    
    
    call displayGUIMessage( 'Test 131' )
    
    allocate(cols(size(iZero)))
    cols = mod(iZero+1,3)+1
    
    open(21,file='cols.txt',status='unknown',form='formatted',action='write')
        do i=1,size(cols)
            write(21,*)  cols(i)
        enddo
    close(21)
    
    open(21,file='XXel.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  XXel(i,1)
        enddo
    close(21)
    
    open(21,file='dims.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  dims(i,1)
        enddo
    close(21)
    
    open(21,file='DimsF2.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  DimsF2(i,1)
        enddo
    close(21)
    
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
    
    allocate(TheseMinXXel(Nel,3), TheseMaxXXel(Nel,3))
    allocate(theseA(size(Xf)),theseB(size(Xf)),theseC(size(Xf)),theseD(size(Xf)))
    allocate(theseA_int(size(Xf)),theseB_int(size(Xf)),theseC_int(size(Xf)),theseD_int(size(Xf)))
    allocate(theseNum(Nel))
    allocate(count1D(Nel))
    TheseMinXXel(:,:) = 0
    
    do n=1,Nel
        count1D(n) = n
    end do
    
    open(21,file='xxMinEl1.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  xxMinEl(i,1)
        enddo
    close(21)
    
    open(21,file='xxMinEl2.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  xxMinEl(i,2)
        enddo
    close(21)
    
    open(21,file='xxMinEl3.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  xxMinEl(i,3)
        enddo
    close(21)
    
    open(21,file='xVertA.txt',status='unknown',form='formatted',action='write')
        do i=1,Nel
            write(21,*)  xVertA(i,1)
        enddo
    close(21)
    
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
        
        if (n == 1) then
            
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
    endif
    
    !    TheseMinXXel = repmat(xxMinEl(n,:),K,1) ;
    !    TheseMaxXXel = repmat(xxMaxEl(n,:),K,1) ;
        theseA = all((TheseMinXXel <= xVertA ) .and. (TheseMaxXXel >= xVertA ),2)
        theseB = all((TheseMinXXel <= xVertB ) .and. (TheseMaxXXel >= xVertB ),2) ;
        theseC = all((TheseMinXXel <= xVertC ) .and. (TheseMaxXXel >= xVertC ),2) ;
        theseD = all((TheseMinXXel <= xVertD ) .and. (TheseMaxXXel >= xVertD ),2) ;
    
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
        
         if (n == 1) then
             open(21,file='theseA.txt',status='unknown',form='formatted',action='write')
                do j=1,size(theseA)
                    write(21,*)  theseA(j)
                enddo
             close(21)
         end if
         
         if (n == 1) then
             open(21,file='theseNum.txt',status='unknown',form='formatted',action='write')
                do j=1,size(theseNum)
                    write(21,*)  theseNum(j)
                enddo
             close(21)
         end if
                
        mask1D = (theseNum >= 2)
        where (mask1D) TheTs(:,n) = .true.
        
        if (n == 1) then
             open(21,file='mask1D1.txt',status='unknown',form='formatted',action='write')
                do j=1,size(mask1D)
                    write(21,*)  mask1D(j)
                enddo
             close(21)
        end if
        
        if (n == 1) then
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

    open(21,file='TheTs.txt',status='unknown',form='formatted',action='write')
    do j=1,Nel 
        do i=1,size(Xf)        
            write(21,*)  TheTs(i,j)
        end do
    end do
    close(21)
    
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