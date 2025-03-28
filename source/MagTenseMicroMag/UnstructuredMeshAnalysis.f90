module UnstructuredMeshAnalysis
  use MKL_SPBLAS
  use BLAS95
  use MicroMagParameters
  use IO_GENERAL
  
  implicit none

    contains

        !>-----------------------------------------
    !> @author Rasmus Bjørk, rabj@dtu.dk, DTU, 2025
    !> Original Matlab implementation by Emil Poulsen and Andrea Roberto Insinga
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

    integer :: Nel, K, idim, ipm, j, n, kb, i, nnz
    integer, allocatable :: TheJ(:,:)
    real(dp), allocatable :: Xel(:), Yel(:), Zel(:)
    real(dp), allocatable :: Volumes(:)
    real(dp) :: DimsScales(3)
    real(dp), allocatable :: XXel(:,:)
    real(dp), allocatable :: Xf(:), Yf(:), Zf(:)
    real(dp), allocatable :: fNormX(:), fNormY(:), fNormZ(:)
    real(dp), allocatable :: AreaFaces(:)
    real(dp), allocatable :: DimsF(:,:)
    integer, allocatable :: ThePM(:), TheEE(:,:), TheNot(:,:)
    integer, allocatable :: Aindex(:), Bindex(:), indexItContainsTrue(:,:)
    integer, allocatable :: ia(:), ja(:)
    real(dp), allocatable :: a(:)
    logical, allocatable :: ItContains(:,:)
    real(dp), allocatable :: dimscopy(:,:)

    ! Initialize some variables
    Nel = size(pos, 1)
    Volumes = dims(:,1) * dims(:,2) * dims(:,3)

    ! Rescaling
    DimsScales = minval(dims, dim=1) / 2.0_dp
    Xel = (pos(:,1) / DimsScales(1))
    Yel = (pos(:,2) / DimsScales(2))
    Zel = (pos(:,3) / DimsScales(3))
    XXel = reshape([Xel, Yel, Zel], [Nel, 3])
    dimscopy = dims / reshape(DimsScales, [1, 3])

     call displayGUIMessage( 'Test 2' )
    ! Construct all faces
    K = 6 * Nel
    allocate(fNormX(K), fNormY(K), fNormZ(K), AreaFaces(K), DimsF(K, 3))
    allocate(TheJ(3, 2), ThePM(2), TheEE(3, 3), TheNot(3, 2))
    ThePM = [-1, 1]
    TheEE = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
    TheNot = reshape([2, 3, 3, 1, 1, 2], [3, 2])

     call displayGUIMessage( 'Test 3' )
    j = 0
    do idim = 1, 3
      do ipm = 1, 2
        j = j + 1
        TheJ(idim, ipm) = j
        fNormX((1 + (j-1) * Nel):(j * Nel)) = TheEE(idim, 1) * ThePM(ipm)
        fNormY((1 + (j-1) * Nel):(j * Nel)) = TheEE(idim, 2) * ThePM(ipm)
        fNormZ((1 + (j-1) * Nel):(j * Nel)) = TheEE(idim, 3) * ThePM(ipm)
        call displayGUIMessage( 'Test 5' )
        Xf((1 + (j-1) * Nel):(j * Nel)) = Xel + TheEE(idim, 1) * ThePM(ipm) * dimscopy(:, 1) / 2.0
        call displayGUIMessage( 'Test 5' )
        Yf((1 + (j-1) * Nel):(j * Nel)) = Yel + TheEE(idim, 2) * ThePM(ipm) * dimscopy(:, 2) / 2.0
        Zf((1 + (j-1) * Nel):(j * Nel)) = Zel + TheEE(idim, 3) * ThePM(ipm) * dimscopy(:, 3) / 2.0
        call displayGUIMessage( 'Test 6' )
        DimsF((1 + (j-1) * Nel):(j * Nel), :) = reshape([dimscopy(:, 1) * (1 - TheEE(idim, 1)), dimscopy(:, 2) * (1 - TheEE(idim, 2)), dimscopy(:, 3) * (1 - TheEE(idim, 3))], [Nel, 3])
        AreaFaces((1 + (j-1) * Nel):(j * Nel)) = (DimsScales(TheNot(idim, 1)) * dimscopy(:, TheNot(idim, 1))) * (DimsScales(TheNot(idim, 2)) * dimscopy(:, TheNot(idim, 2)))
        call displayGUIMessage( 'Test 7' )
      end do
    end do
    
    call displayGUIMessage( 'Test 4' )

    open(21,file='AreaFaces.txt',status='unknown',form='formatted',action='write')
    do i=1,j
        write(21,*)  AreaFaces(i)
    enddo
    close(21)
    
    ! Check which faces are contained by other faces
    allocate(ItContains(K, K))
    ItContains = .false.
    allocate(indexItContainsTrue(0, 2))

    !do idim = 1, 3
    !  Aindex = (1 + (TheJ(idim, 1) - 1) * Nel):(TheJ(idim, 1) * Nel)
    !  Bindex = (1 + (TheJ(idim, 2) - 1) * Nel):(TheJ(idim, 2) * Nel)
    !  do kb = 1, Nel
    !    do i = 1, size(Aindex)
    !      if (XXf(Aindex(i), idim) == XXf(Bindex(kb), idim)) then
    !        if ((XXf(Aindex(i), TheNot(idim, 1)) - 0.5_dp * DimsF(Aindex(i), TheNot(idim, 1)) <= XXf(Bindex(kb), TheNot(idim, 1)) - 0.5_dp * DimsF(Bindex(kb), TheNot(idim, 1))) .and. &
    !            (XXf(Aindex(i), TheNot(idim, 1)) + 0.5_dp * DimsF(Aindex(i), TheNot(idim, 1)) >= XXf(Bindex(kb), TheNot(idim, 1)) + 0.5_dp * DimsF(Bindex(kb), TheNot(idim, 1))) .and. &
    !            (XXf(Aindex(i), TheNot(idim, 2)) - 0.5_dp * DimsF(Aindex(i), TheNot(idim, 2)) <= XXf(Bindex(kb), TheNot(idim, 2)) - 0.5_dp * DimsF(Bindex(kb), TheNot(idim, 2))) .and. &
    !            (XXf(Aindex(i), TheNot(idim, 2)) + 0.5_dp * DimsF(Aindex(i), TheNot(idim, 2)) >= XXf(Bindex(kb), TheNot(idim, 2)))) then
    !          allocate(indexItContainsTrue(size(indexItContainsTrue, 1) + 1, 2))
    !          indexItContainsTrue(size(indexItContainsTrue, 1), :) = [Aindex(i), Bindex(kb)]
    !        end if
    !      end if
    !    end do
    !  end do
    !end do

    !nnz = size(indexItContainsTrue, 1)
    !allocate(ia(K + 1), ja(nnz), a(nnz))
    !ia = 0
    !ja = indexItContainsTrue(:, 2)
    !a = 1.0_dp

    !call mkl_sparse_d_create_csr(GridInfo%TheSigns, SPARSE_INDEX_BASE_ZERO, K, K, ia, ia(2:), ja, a)

    ! Remove faces containing other faces
    ! Your complete implementation goes here

    ! Construct T and D matrices
    ! Your complete implementation goes here

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