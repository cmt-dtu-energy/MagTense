MODULE MagPrism_CALL
use PARAMETERS_CALL
use UTIL_CALL
use qsort_c_module
implicit none
    contains
    
!::Needs a better name
!::pos is a (n,3) array with the x,y,z coordinates for each prism
!::
!::dims is a (n,3) array with the dimensions (a,b,c) for each prism
!::
!::dir is a (n,3) array with the direction unit vector for each prism
!::
!::magType is a (n,1) array with indication of the type of magnet 
!::(defined in the constants so that 0 is a soft ferromagnet and 1
!::is a hard ferromagnet). Hard magnets are dominated by the direction
!:: while the direction of soft magnets is determined by the local field
!::
!:: Hext has size(n,1) and contains an externally applied field
!::
!::T has size(n,1) and contains the temperature at each prism
!::
!::n is the number of prisms
!::
!:: stateFcn is an array of custom type stateFunction. Each element in this array thus
!:: contains a state function (M(T,H))
!::
!:: stateFcnIndices is a (n,1) size array where each element gives the index to the
!:: state function for the corresponding prism. Thus, the maximum value of stateFcnIndices
!:: should not exceed m
!::
!:: m is the number of state functions available
!::
!:: solPts is an (l,3) array specifying the points at which the solution is to be returned
!::
!:: l is the number of points in the points at which the solution is required.
!::
!:: spacedim is the spatial dimension of the solution. It can be 2 for 2D or 3 for 3D
!::
!:: Hout is an (n+l,3) array with the magnetic field returned in each prism (the first m points)
!:: and the required points (from n+1 to n+l)
!::
!:: Mout is similar to Hout but containing the magnetization and only for the n prisms, i.e. Mout 
!:: has the size (n,3)
!::
!:: nIte is an integer defining the max number of iterations
!::
!:: rotAngles is an (n,3) array with the rotation angles (yaw, pitch and roll) for each prism
!::    
!:: hyst_map_init is an(n,1) array with zeros and ones indicating whether the state of the i'th unit
!:: is to use the cooling curve (zero) or the heating curve (one)
!::
!:: formType is a (n,1) array with indication of the type of shape
!:: that the individual object has. It is defined in the constants 
!:: so that 0 is a prism and 1 is an ellipsoid
!:: NOTE THAT IN 2D THE CODE ONLY WORKS FOR PRISMS
!::
!:: 13-07-2017 rabj:
!:: Added ellipsoids (spheroids) in 3D
!:: 3-11-2016 kaki:
!:: hyst_map returns the hystereis map as found by the iterated solution and has the same size as hyst_map_init
!:: 4-11-2016 kaki
!:: added support for defining an initial guess
!:: Hint_init is an optional initial field (starting guess). If present it is used in stead of the applied field (Hext) as
!:: the initial value for Hint
subroutine calcH( pos, dims, dir, magType, formType, Hext, T, n, stateFcn, stateFcnIndices, m, solPts, l, spacedim, Hout, Mout, nIte, rotAngles, hyst_map_init, hyst_map, Hint_init )
real,intent(in),dimension(n,3) :: pos, dims,dir,rotAngles
real,intent(in),dimension(n+l,3) :: Hext
integer,intent(in),dimension(n) :: magType, formType
real,intent(in),dimension(n) :: T
type(stateFunction),intent(in),dimension(2*m) :: stateFcn
integer,intent(in),dimension(n) :: stateFcnIndices
integer,intent(in) :: n,m,l,spacedim
real,intent(in),dimension(l,3) :: solPts
real,intent(out),dimension(:,:),allocatable :: Hout,Mout
integer,intent(inout) :: nIte
integer,intent(in),dimension(n) :: hyst_map_init
integer,intent(inout),dimension(n) :: hyst_map
real,intent(in),dimension(n,3),optional :: Hint_init

real :: lambda,maxRelDiff,lambda_old
real,dimension(4) :: maxRelDiffArr
real,dimension(:,:),allocatable :: Hdem,H_old
integer :: i,j,convCnt,lambdaCnt
integer,parameter :: maxIte=1000
real,parameter :: tol = 0.01
logical :: lCh
real,dimension(:),allocatable :: tmp,Hnorm,Hnorm_old
real,dimension(:,:,:),allocatable :: rotMat,rotMatInv
real,dimension(:,:,:,:),allocatable :: N_out
integer :: debugRetval,endInd
integer,parameter :: recl=2
real,parameter :: conv_fractile=0.85
real,dimension(2) :: debArr
lambda = 1.0
lambda_old = lambda
convCnt = 1
lambdaCnt = 1

call writeVersionString()

!::check that maxval(stateFcnIndices) <= m
if ( maxval(stateFcnIndices) .le. m ) then
    !::Allocate the solution
    allocate(Hout(n+l,3),Mout(n,3),Hdem(n,3),H_old(n,3))
    allocate(tmp(n),Hnorm(n),Hnorm_old(n))
    Hout(:,:) = 0
    if ( present( Hint_init ) ) then
        Hout(1:n,:) = Hint_init
    else        
        Hout(1:n,:) = Hext(1:n,:)
    endif
    
    Mout(:,:) = 0
    H_old(:,:) = 0
    hyst_map = hyst_map_init
    
    !::Get the rotation matrix and its inverse
    call findRotationMatrices( rotAngles, n, rotMat, rotMatInv)
    
    call calcNTensor( pos, dims, N_out, n, spacedim, rotMat, formType )
    
    !::First we find the self-consistent solution, i.e. the magnetization (magnitude and direction) in each prism    
    do i=0,maxIte
        nIte = i
        H_old = Hout(1:n,:)
        !::Get the magnetization in each node
        call getM(Hout(1:n,:), T, Mout, dir, magType, n, stateFcn, stateFcnIndices, m, hyst_map, rotMat)
        
        !::Get the demagnetizing field
        call getHdem( Mout, n, Hdem, N_out, rotMat, rotMatInv)
        
        !::Update the solution
        Hout(1:n,:) = Hout(1:n,:) + lambda * ( (Hext(1:n,:) + Hdem) - Hout(1:n,:) )
        
        !::update the hysteresis tracking
        call updateHysteresisTracking( T, Hout, hyst_map_init, hyst_map, n, stateFcn, stateFcnIndices, m )
        
        !::Find the relative change in the solution
        Hnorm = sqrt( Hout(1:n,1)**2 + Hout(1:n,2)**2 + Hout(1:n,3)**2 )
        Hnorm_old = sqrt( H_old(1:n,1)**2 + H_old(1:n,2)**2 + H_old(1:n,3)**2 );
        tmp = -1
         where( Hnorm_old .ne. 0 )
            tmp = abs( ( Hnorm - Hnorm_old ) / Hnorm_old )
        endwhere
        
        call QsortC(tmp)
        
        endInd = ceiling( conv_fractile * n )
        
        
        maxRelDiff = maxval(tmp)
        maxRelDiffArr = cshift( maxRelDiffArr, 1 )
        maxRelDiffArr(1) = maxRelDiff
        !call writeDebugString( 'maxreldiff', i, maxRelDiff )
        lambda_old = lambda
        if ( i .ge. 4 .AND. convCnt .eq. 1 .AND. lambdaCnt .ge. 4 ) then
            call writeDebugStringArr1D( 'UPD LAMBDA', maxRelDiffArr )
            call updateLambda( lambda, maxRelDiffArr, lCh )
            if ( lambda .ne. lambda_old ) then
                lambdaCnt = 1
            endif
            
        endif
        call writeDebugString( 'lambda', i, lambda  )
        lambdaCnt = lambdaCnt + 1
        !write(*,*) 'maxRelDiff,lambda,tol',maxRelDiff,lambda,tol*lambda
        !call writeDebugStringArr1D( 'Mout ', sqrt(Mout(:,1)**2+Mout(:,2)**2+Mout(:,3)**2) )
        !call writeDebugStringArr1D( 'Hout ', sqrt(Hout(:,1)**2+Hout(:,2)**2+Hout(:,3)**2) )
        debArr(1) = tol * lambda
        debArr(2) = tmp(endInd)
        call writeDebugStringArr1D( 'lambda', debArr )
        
        if (  tmp(endInd) .lt. tol * lambda ) then
            call writeDebugString('SolConv')            
            exit    
        endif        
        
        if ( maxRelDiff .lt. tol * lambda  ) then
            if ( lambda .ne. 1 .AND. convCnt .lt. 5 ) then
                convCnt = convCnt + 1
            elseif ( i .gt. 2 ) then         
                call writeDebugString('SolConv')
                call writeDebugString('maxRelDiff',i,maxRelDiff)
                
                
                
                exit
            endif
        else
                convCnt = 1
        endif
       
        
        !open (14, file='H_M.dat',	&
			     !      status='unknown', form='unformatted',action='write',	&
			     !      access='direct', recl=recl*n*3)
        !write(14,rec=2*i+1) Hout(1:n,:)
        !write(14,rec=2*i+2) Mout
        !close(14)
        
    enddo
    
    
    !::Having obtained the solution we compute the magnetic field at each of the required points
    call getSolution( solPts, l, Mout, dims, pos, n, spacedim, Hout(n+1:n+l,:), Hext(n+1:n+l,:), rotMat, rotMatInv, formType )
    
endif

end subroutine calcH

!::Returns the solution at the required points
subroutine getSolution( solPts, l, M, dims, pos, n, spacedim, Hout, Hext, rotMat, rotMatInv, formType )
real,intent(in),dimension(l,3) :: solPts
real,intent(in),dimension(l,3) :: Hext
real,intent(inout),dimension(l,3) :: Hout
real,intent(in),dimension(n,3) :: M, dims, pos
real,intent(in),dimension(n,3,3) :: rotMat,rotMatInv
integer,intent(in),dimension(n) :: formType
integer,intent(in) :: l,n, spacedim
real,dimension(3) :: diffPos, dotProd
real,dimension(3,3) :: N_out
integer :: i,j


!::Loop over each point where the solution is required
do i=1,l
    !::Loop over each prism
    do j=1,n
        diffPos = pos(j,:) - solPts(i,:)
        !::Rotate the position vector to accomodate rotation of the j'th prism
        diffPos = matmul( rotMat(j,:,:), diffPos )
        if ( spacedim .eq. 3 ) then
           call getN_3D(dims(j,1),dims(j,2),dims(j,3), diffPos(1), diffPos(2), diffPos(3), formType(j), N_out )
        endif
        if ( spacedim .eq. 2 ) then
           call getN_2D(dims(j,1),dims(j,2),dims(j,3), diffPos(1), diffPos(2), diffPos(3), N_out )
        endif
        
        call getDotProd( N_out, matmul( rotMat(j,:,:), M(j,:) ), dotProd )
                
        !::Rotate the solution back. 
        dotProd = matmul( rotMatInv(j,:,:), dotProd )
        
        Hout(i,:) = Hout(i,:) - dotProd       
    enddo
    
    Hout = Hout + Hext
    
enddo


end subroutine getSolution


!:: divides lambda by 2 and sets lCh to true, if maxDiff oscillates.
subroutine updateLambda( lambda, maxDiff, lCh )
real,intent(inout) :: lambda
real,intent(in),dimension(4) :: maxDiff
logical,intent(inout) :: lCh

lCh = .false.

if ( (maxDiff(1) .gt. maxDiff(2) .AND. maxDiff(2) .lt. maxDiff(3) .AND. maxDiff(3) .gt. maxDiff(4)) .OR. &
     (maxDiff(1) .lt. maxDiff(2) .AND. maxDiff(2) .gt. maxDiff(3) .AND. maxDiff(3) .lt. maxDiff(4))   ) then

    lambda = lambda * 0.95
    lCh = .true.
endif


end subroutine updateLambda

!::Returns the magnitude of the magnetization in each prism
!::3 Nov 2016. kaki updated this function to include a hyst_map,
!::which has the size of the number of units (n). If hyst_map(i) is
!::zero then the cooling curve is used while if hyst_map(i) is one then
!::the heating curve is used. It is assumed that the stateFcn array has
!::two state functions for each material modeled, i.e. 2*m elements. The first
!:: m elements are cooling curves while the elements from m+1 to 2m are heating curves
subroutine getM( H, T, Mout, dir, magType, n, stateFcn, stateFcnIndices, m, hyst_map, rotMat)
real,intent(in),dimension(n,3) :: H,dir
real,intent(in),dimension(n) :: T
real,intent(inout),dimension(n,3) :: Mout
real,intent(in),dimension(n,3,3) :: rotMat
integer,intent(in),dimension(n) :: magType
integer,intent(in) :: n,m
type(stateFunction),intent(in),dimension(2*m) :: stateFcn
integer,intent(in),dimension(n) :: stateFcnIndices
integer,intent(in),dimension(n) :: hyst_map

real,dimension(:),allocatable :: Hnorm

real :: Mtmp
integer :: i

allocate(Hnorm(n))

Hnorm = sqrt( H(:,1)**2 + H(:,2)**2 + H(:,3)**2 )


!::For each prism
!$OMP PARALLEL DO PRIVATE(i,Mtmp)
do i=1,n
    call getBilinInterp( stateFcn(stateFcnIndices(i)+hyst_map(i)*m)%M, stateFcn(stateFcnIndices(i)+hyst_map(i)*m)%T, &
                         stateFcn(stateFcnIndices(i)+hyst_map(i)*m)%H, stateFcn(stateFcnIndices(i)+hyst_map(i)*m)%nT, &
                         stateFcn(stateFcnIndices(i)+hyst_map(i)*m)%nH, T(i), Hnorm(i), Mtmp )    
    if ( magType(i) .eq. magTypeSoft ) then
        !::Then the magnetization is along the direction of the local field
        if ( Hnorm(i) .ne. 0 ) then
            Mout(i,:) = Mtmp * H(i,:) / Hnorm(i)                  
        else
            Mout(i,1) = Mtmp
            Mout(i,2:3) = 0        

        endif        
    elseif ( magType(i) .eq. magTypeHard ) then
        !::Then the magnetization is along the direction pre-specified
        Mout(i,:) = Mtmp * dir(i,:)
    endif
    
enddo
!$OMP END PARALLEL DO


end subroutine getM

!::Updates the hysteresis map based on the temperature, magnetic field norm and phase-diagram
!:: T is an (n,1) array with the temperatures
!:: H is an (n,3) array with the magnetic vector field
!:: hyst_map_init is an (n,1) array with the initial hysteresis map
!:: hyst_map is an (n,1) array with the current hysteresis map
!:: n is an integer indicating the number of units
!:: stateFcn is an array of the custom type stateFcn with 2*m elements containing the state functions used
!:: stateFncIndices is an (n,1) integer array with the indices mapping from the n units into the m state functions
!:: m is an integer 
subroutine updateHysteresisTracking( T, H, hyst_map_init, hyst_map, n, stateFcn, stateFcnIndices, m )
real,dimension(n),intent(in) :: T
real,dimension(n,3),intent(in) :: H
integer,dimension(n),intent(in) :: hyst_map_init
integer,dimension(n),intent(inout) :: hyst_map
integer,intent(in) :: n,m
type(stateFunction),intent(in),dimension(2*m) :: stateFcn
integer,intent(in),dimension(n) :: stateFcnIndices

integer :: i
real :: Hcrit_cool,Hcrit_heat
real,dimension(:),allocatable :: Hnorm

allocate(Hnorm(n))


Hnorm = sqrt( H(:,1)**2+H(:,2)**2+H(:,3)**2 )

!:: for each element
do i=1,n
    call interp1( stateFcn(stateFcnIndices(i))%Tcrit_cool, stateFcn(stateFcnIndices(i))%Hcrit_cool, T(i), stateFcn(stateFcnIndices(i))%nCritCool, Hcrit_cool )
    
    call interp1( stateFcn(stateFcnIndices(i))%Tcrit_heat, stateFcn(stateFcnIndices(i))%Hcrit_heat, T(i), stateFcn(stateFcnIndices(i))%nCritHeat, Hcrit_heat )
    
    
    if ( Hnorm(i) .ge. Hcrit_cool ) then
        !::Then the heating curve should be used
        hyst_map(i) = 1
    elseif ( Hnorm(i) .le. Hcrit_heat ) then
        !::Then the cooling curve should be used
        hyst_map(i) = 0
    else
        !::then we are in the mixed region and the state should be the initial state
        hyst_map(i) = hyst_map_init(i)
    endif
    
    
    
enddo

!call writeDebugStringArr1DInt('Hystm', hyst_map )

deallocate(Hnorm)
end subroutine updateHysteresisTracking


!::Find the demagnetizing field
subroutine getHdem( M, n, Hdem, N_out, rotMat, rotMatInv)
real,intent(in),dimension(n,3) :: M
integer,intent(in) :: n
real,intent(in),dimension(n,3,3) :: rotMat, rotMatInv
real,intent(inout),dimension(n,3) :: Hdem
real,intent(in),dimension(n,n,3,3) :: N_out
real,dimension(3) :: dotProd
integer :: i,j

Hdem(:,:) = 0

!$OMP PARALLEL DO PRIVATE(i,j,dotProd)
do i=1,n
    !::For each prism find the contribution from all other prisms
    do j=1,n        

        !::Compute the dot product between N and M
        !call getDotProd( N_out(i,j,:,:), M(j,:), dotProd )
        
        !::Rotate the solution back.
        dotProd = matmul( rotMatInv(j,:,:), dotProd )        
        
        !::Add to the solution
        Hdem(i,:) = Hdem(i,:) - dotProd
                
    enddo        
enddo
!$OMP END PARALLEL DO

end subroutine getHdem

subroutine calcNTensor( pos, dims, N_out, n, spacedim, rotMat, formType )
real,intent(in),dimension(n,3) :: dims,pos
real,intent(inout),allocatable,dimension(:,:,:,:) :: N_out
real,intent(in),dimension(n,3,3) :: rotMat
integer,intent(in),dimension(n) :: formType
integer,intent(in) :: n, spacedim
integer :: i,j
real,dimension(3) :: diffPos

allocate( N_out(n,n,3,3) )

!$OMP PARALLEL DO PRIVATE(i,j,diffPos)
do i=1,n
    !::For each prism find the contribution from all other prisms
    do j=1,n        
        diffPos(1) = pos(i,1)-pos(j,1)
        diffPos(2) = pos(i,2)-pos(j,2)
        diffPos(3) = pos(i,3)-pos(j,3)
        
        !::Rotate the postion vector
        diffPos = matmul( rotMat(j,:,:), diffPos )
        
        !::Find the demag. tensor
        if ( spacedim .eq. 3 ) then
            call getN_3D(dims(j,1),dims(j,2),dims(j,3),diffPos(1),diffPos(2),diffPos(3), formType(j), N_out(i,j,:,:) )
        elseif ( spacedim .eq. 2 ) then
            call getN_2D(dims(j,1),dims(j,2),dims(j,3),diffPos(1),diffPos(2),diffPos(3), N_out(i,j,:,:) )
        endif
    enddo
enddo
!$OMP END PARALLEL DO

end subroutine calcNTensor

!::
!::Calculates the rotation matrix and its inverse for each prism
!::
subroutine findRotationMatrices( rotAngles, n, rotMat, rotMatInv)
real,intent(in),dimension(n,3) :: rotAngles
integer,intent(in) :: n
real,intent(inout),dimension(:,:,:),allocatable :: rotMat,rotMatInv

integer :: i
real,dimension(3,3) :: RotX,RotY,RotZ

allocate( rotMat(n,3,3), rotMatInv(n,3,3) )

do i=1,n
    !::The minus sign is important since a rotated prism can be represented with a rotation about the given axis in the opposite direction
    call getRotX( -rotAngles(i,1), RotX )
    call getRotY( -rotAngles(i,2), RotY )
    call getRotZ( -rotAngles(i,3), RotZ )

    !::Find the rotation matrix as defined by yaw, pitch and roll
    !::Rot = RotZ * RotY * RotX
    rotMat(i,:,:) = matmul( matmul( RotZ, RotY ), RotX )
    !::As the rotation matrices are orthogonal we have that înv(R) = transp(R)
    !:: And thus inv(R1R2) = transp(R2)transp(R1) (note the change in multiplication order)
    !::and that transp( R1R2R3 ) = transp(R3)transp(R1R2) = transp(R3) transp(R2) transp(R1)
    !:: inv(Rot) = inv(RotZ * RotY * RotX ) =>
    !:: inv(Rot) = transp(Rot) = transp( RotZ * RotY * RotX )
    !::                        = transp( RotX ) * transp( RotZ * RotY )
    !::                        = transp( RotX ) * transp( RotY ) * transp( RotZ )
    rotMatInv(i,:,:) = matmul( matmul( transpose(RotX), transpose(RotY) ), transpose( RotZ ) )
enddo

end subroutine findRotationMatrices

!::
!::Returns the rotation matrix for a rotation of angle radians about the x-axis
!::
subroutine getRotX( angle, rotMat )
real,intent(in) :: angle
real,intent(inout),dimension(3,3) :: rotMat

!::fortran matrices are (row,col)
rotMat(1,1) = 1
!::Top row
rotMat(1,2:3) = 0
!::left column
rotMat(2:3,1) = 0
rotMat(2,2) = cos(angle)
rotMat(3,3) = cos(angle)
rotMat(2,3) = -sin(angle)
rotMat(3,2) = sin(angle)

end subroutine getRotX

!::
!::Returns the rotation matrix for a rotation of angle radians about the y-axis
!::
subroutine getRotY( angle, rotMat )
real,intent(in) :: angle
real,intent(inout),dimension(3,3) :: rotMat

!::fortran matrices are (row,col)
!::top row
rotMat(1,1) = cos(angle)
rotMat(1,2) = 0
rotMat(1,3) = sin(angle)
!::middle row
rotMat(2,1) = 0
rotMat(2,2) = 1
rotMat(2,3) = 0
!::bottom row
rotMat(3,1) = -sin(angle)
rotMat(3,2) = 0
rotMat(3,3) = cos(angle)

end subroutine getRotY

!::
!::Returns the rotation matrix for a rotation of angle radians about the z-axis
!::
subroutine getRotZ( angle, rotMat )
real,intent(in) :: angle
real,intent(inout),dimension(3,3) :: rotMat

!::top row
rotMat(1,1) = cos(angle)
rotMat(1,2) = -sin(angle)
rotMat(1,3) = 0
!::middle row
rotMat(2,1) = sin(angle)
rotMat(2,2) = cos(angle)
rotMat(2,3) = 0
!::bottom row
rotMat(3,1) = 0
rotMat(3,2) = 0
rotMat(3,3) = 1

end subroutine getRotZ


!:: Calculates the dot product between (3x3)-matrix and (3)-vector
subroutine getDotProd( N, M, dot_prod )
real,intent(in),dimension(3,3) :: N
real,intent(in),dimension(3) :: M
real,intent(inout),dimension(3) :: dot_prod

dot_prod(1) = sum(N(1,:) * M(:))
dot_prod(2) = sum(N(2,:) * M(:))
dot_prod(3) = sum(N(3,:) * M(:))


end subroutine getDotProd


!::Calculates N from the analytical expression in 3D
!::Given the dimensions of the prism (a,b,c) and the distance vector
!::to it
subroutine getN_3D( a, b, c, x, y, z, form, N_out )
real,intent(in) :: a,b,c,x,y,z
integer,intent(in) :: form
real,intent(out),dimension(3,3) :: N_out
real :: nom,denom,nom_l,nom_h,denom_l,denom_h
!real,parameter :: zero_limit=1.006258
real,parameter :: lim_scl_h=1.01,lim_scl_l=0.98
real :: xl,xh,yl,yh,zl,zh,al,ah,bl,bh,cl,ch,lim
real :: a_temp,b_temp,c_temp,a_o,b_o,c_o,x_o,y_o,z_o,f,w,A_c,alpha
real,dimension(6) :: tmp
integer,dimension(3) :: indx_semi
real,dimension(3) :: coor,semi_axis
real,dimension(3,3) :: N_demag, N_temp

!--- The demagnetization factor of a prism
if ( form .eq. formPrisme ) then
    !::Diagonal elements
    N_out(1,1) = 1./(4.*pi) * ( atan(f_3D(a,b,c,x,y,z))   + atan(f_3D(a,b,c,-x,y,z))  + atan(f_3D(a,b,c,x,-y,z)) + &
                                atan(f_3D(a,b,c,x,y,-z))  + atan(f_3D(a,b,c,-x,-y,z)) + atan(f_3D(a,b,c,x,-y,-z)) + &
                                atan(f_3D(a,b,c,-x,y,-z)) + atan(f_3D(a,b,c,-x,-y,-z)) )



    N_out(2,2) = 1./(4.*pi) * ( atan(g_3D(a,b,c,x,y,z))   + atan(g_3D(a,b,c,-x,y,z))  + atan(g_3D(a,b,c,x,-y,z)) + &
                                atan(g_3D(a,b,c,x,y,-z))  + atan(g_3D(a,b,c,-x,-y,z)) + atan(g_3D(a,b,c,x,-y,-z)) + &
                                atan(g_3D(a,b,c,-x,y,-z)) + atan(g_3D(a,b,c,-x,-y,-z)) )
                            
                            

    N_out(3,3) = 1./(4.*pi) * ( atan(h_3D(a,b,c,x,y,z))   + atan(h_3D(a,b,c,-x,y,z))  + atan(h_3D(a,b,c,x,-y,z)) + &
                                atan(h_3D(a,b,c,x,y,-z))  + atan(h_3D(a,b,c,-x,-y,z)) + atan(h_3D(a,b,c,x,-y,-z)) + &
                                atan(h_3D(a,b,c,-x,y,-z)) + atan(h_3D(a,b,c,-x,-y,-z)) )                            
                            

    !::Off-diagonal elements
    nom = FF_3D(a,b,c,x,y,z)  * FF_3D(-a,-b,c,x,y,z) * FF_3D(a,-b,-c,x,y,z) * FF_3D(-a,b,-c,x,y,z)
    denom = FF_3D(a,-b,c,x,y,z) * FF_3D(-a,b,c,x,y,z)  * FF_3D(a,b,-c,x,y,z)  * FF_3D(-a,-b,-c,x,y,z)

    if ( denom .eq. 0 .or. nom .eq. 0 ) then
        !Find the limit
        lim = getF_limit(a,b,c,x,y,z,FF_3D)
        N_out(1,2) = -1./(4.*pi) * log( lim )        
    else
        N_out(1,2) = -1./(4.*pi) * log( nom / denom )
    endif



    N_out(2,1) = N_out(1,2)

    nom = GG_3D(a,b,c,x,y,z)  * GG_3D(-a,-b,c,x,y,z) * GG_3D(a,-b,-c,x,y,z) * GG_3D(-a,b,-c,x,y,z)
    denom = GG_3D(a,-b,c,x,y,z) * GG_3D(-a,b,c,x,y,z)  * GG_3D(a,b,-c,x,y,z)  * GG_3D(-a,-b,-c,x,y,z)

    if ( denom .eq. 0 .or. nom .eq. 0 ) then
    !Find the limit
      
        lim = getF_limit(a,b,c,x,y,z,GG_3D)       
    
        N_out(2,3) = -1./(4.*pi) * log( lim )
    else
        N_out(2,3) = -1./(4.*pi) * log( nom/denom )
    endif


    N_out(3,2) = N_out(2,3)

    nom = HH_3D(a,b,c,x,y,z)  * HH_3D(-a,-b,c,x,y,z) * HH_3D(a,-b,-c,x,y,z) * HH_3D(-a,b,-c,x,y,z)
    denom = HH_3D(a,-b,c,x,y,z) * HH_3D(-a,b,c,x,y,z)  * HH_3D(a,b,-c,x,y,z)  * HH_3D(-a,-b,-c,x,y,z)

    if ( denom .eq. 0 .or. nom .eq. 0 ) then
    
        lim = getF_limit(a,b,c,x,y,z,HH_3D)           
    
        N_out(1,3) = -1./(4.*pi) * log( lim )
    else
        N_out(1,3) = -1./(4.*pi) * log( nom / denom )
    endif


    N_out(3,1) = N_out(1,3)
endif

!--- The demagnetization factor of ellipsoids
!--- At present only spheroids are implemented, i.e. 
!--- ellipsoids where two of the semi-axis are equal
!---
!--- The general ellipsoid with a .neq. b .neq. c
!--- can also be implemented, but involved computing
!--- incomplete elliptical integrals
if ( form .eq. formEllipse ) then
    semi_axis(1:3) = (/ a, b, c /)
    coor(1:3) = (/ x, y, z /)

    !--- Manual sort
    !--- We can have the following cases:
    !--- a = b > c  %Do nothing
    !--- a = b < c  %Switch a and c
    !--- a < b = c  %Switch a and c
    !--- a > b = c  %Do nothing
    !--- a > b < c  %Switch b and c
    !--- a < b > c  %Switch b and a

    indx_semi(1:3) = (/ 1, 2, 3 /)

    !--- a = b < c  %Switch a and c
    if ((a == b .AND. b < c) .OR. (a < b .AND. b == c)) then
        indx_semi(1) = 3
        indx_semi(3) = 1
    end if

    !--- a > b < c  %Switch b and c
    if (a > b .AND. b < c) then
        indx_semi(2) = 3
        indx_semi(3) = 2
    end if

    !--- a < b > c  %Switch b and a
    if (a < b .AND. b > c) then
        indx_semi(1) = 2
        indx_semi(2) = 1
    end if

    a_o = semi_axis(indx_semi(1))
    b_o = semi_axis(indx_semi(2))
    c_o = semi_axis(indx_semi(3))
    
    x_o = coor(indx_semi(1))
    y_o = coor(indx_semi(2))
    z_o = coor(indx_semi(3))

    !--- Internal field
    if (x == 0 .AND. y == 0 .AND. z == 0) then
        !--- Oblate
        if (a_o == b_o) then
            alpha = c_o/a_o
        
            N_demag(3,3) = 1./(1.-alpha**2)*(1.-alpha/(sqrt(1.-alpha**2))*acos(alpha))
            N_demag(1,1) = (1.-N_demag(3,3))/2.
            N_demag(2,2) = N_demag(1,1)
        end if
    
        !--- Prolate
        if (b_o == c_o) then
            alpha = a_o/c_o;
        
            N_demag(1,1) = 1./(alpha**2-1.)*(alpha/(sqrt(alpha**2-1.))*acosh(alpha)-1)
            N_demag(2,2) = (1.-N_demag(1,1))/2.
            N_demag(3,3) = N_demag(2,2)
        end if
    
        !--- Sphere
        if (a_o == b_o .AND. b_o == c_o) then
            N_demag(1,1) = 1./3.
            N_demag(2,2) = 1./3.
            N_demag(3,3) = 1./3.
        end if
    
        N_demag(1,2) = 0
        N_demag(1,3) = 0
        N_demag(2,1) = 0
        N_demag(2,3) = 0
        N_demag(3,1) = 0
        N_demag(3,2) = 0
        
    !--- External field    
    else
    
        !--- Oblate
        if (a_o == b_o) then
            f   = sqrt(a_o**2-c_o**2);
            A_c = sqrt((x_o**2+y_o**2+z_o**2-f**2)**2+4*z_o**2*f**2);
            w   = sqrt(1./(2*f**2)*(x_o**2+y_o**2+z_o**2+f**2+A_c));

            N_demag(1,1:3) = (/ -(1./2*(asin(1./w)-sqrt(w**2-1)/w**2)-x_o**2*sqrt(w**2-1)/(w**4*A_c)) , x_o*y_o*sqrt(w**2-1)/(w**4*A_c) , z_o*x_o/(w**2*sqrt(w**2-1)*A_c) /);
            N_demag(2,1:3) = (/ x_o*y_o*sqrt(w**2-1)/(w**4*A_c) , -(1./2*(asin(1./w)-sqrt(w**2-1)/w**2)-y_o**2*sqrt(w**2-1)/(w**4*A_c)) , z_o*y_o/(w**2*sqrt(w**2-1)*A_c) /);
            N_demag(3,1:3) = (/ x_o*z_o/(w**2*sqrt(w**2-1)*A_c) , y_o*z_o/(w**2*sqrt(w**2-1)*A_c) , -(1./sqrt(w**2-1)-asin(1./w)-z_o**2./((w**2-1)**(3./2.)*A_c)) /);
            N_demag = N_demag*c_o*a_o**2*f**(-3) 
        end if

        !--- Prolate
        if (b_o == c_o) then
            f   = sqrt(a_o**2-c_o**2);
            A_c = sqrt((x_o**2+y_o**2+z_o**2+f**2)**2-4*x_o**2*f**2);
            w   = sqrt(1./(2*f**2)*(x_o**2+y_o**2+z_o**2+f**2+A_c));

            N_demag(1,1:3) = (/ -(1./2*log((w+1)/(w-1))-1./w-x_o**2./(w**3*A_c)) , y_o*x_o/(w*(w**2-1)*A_c) , z_o*x_o/(w*(w**2-1)*A_c) /)
            N_demag(2,1:3) = (/x_o*y_o/(w*(w**2-1)*A_c) , -(1./2*(w/(w**2-1)+1./2*log((w-1)/(w+1)))-y_o**2*w/((w**2-1)**2*A_c)) , z_o*y_o*w/((w**2-1)**2*A_c) /)
            N_demag(3,1:3) = (/x_o*z_o/(w*(w**2-1)*A_c) , z_o*y_o*w/((w**2-1)**2*A_c) , -(1./2*(w/(w**2-1)+1./2*log((w-1)/(w+1)))-z_o**2*w/((w**2-1)**2*A_c)) /)
            N_demag = N_demag*a_o*c_o**2*f**(-3)
        end if
    
        !--- Sphere
        if (a_o == b_o .AND. b_o == c_o) then
            N_demag(1,1:3) = (/ x_o**2-sqrt(x_o**2+y_o**2+z_o**2)**2/3., x_o*y_o, x_o*z_o /)
            N_demag(2,1:3) = (/ x_o*y_o, y_o**2-sqrt(x_o**2+y_o**2+z_o**2)**2/3., y_o*z_o /)
            N_demag(3,1:3) = (/ x_o*z_o, y_o*z_o, z_o**2-sqrt(x_o**2+y_o**2+z_o**2)**2/3. /)
            N_demag = 3./(4.*pi*sqrt(x_o**2+y_o**2+z_o**2)**5)*N_demag
            N_demag = 4./3.*pi*a**3*N_demag;
        end if
    
        N_demag = -N_demag
    end if

    !--- Reorder N based on indx_semi
    !--- Switch rows and columns
    N_temp = N_demag
    N_demag(1:3,1) = N_temp(1:3,indx_semi(1));
    N_demag(1:3,2) = N_temp(1:3,indx_semi(2));
    N_demag(1:3,3) = N_temp(1:3,indx_semi(3));
    N_temp = N_demag;
    N_demag(1,1:3) = N_temp(indx_semi(1),1:3);
    N_demag(2,1:3) = N_temp(indx_semi(2),1:3);
    N_demag(3,1:3) = N_temp(indx_semi(3),1:3);
    
    N_out = N_demag    
endif

!open(11,file='debug.txt',status='unknown',access='sequential',form='formatted',position='append',action='write')
!write(11,*) "N_out"
!write(11,*) N_out
!write(11,*) a
!write(11,*) b
!write(11,*) c
!write(11,*) x
!write(11,*) y
!write(11,*) z
!close(11)
    
!call writeDebugStringArrD2( 'N_out ', N_out)

end subroutine getN_3D

function getF_limit(a,b,c,x,y,z,func)
real,intent(in) :: x,y,z,a,b,c
real :: getF_limit
real :: xh,yh,zh,xl,yl,zl,nom_l,nom_h,denom_l,denom_h,nom,denom
real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999
real :: func
external func
   
    !Find the limit
    xh = x * lim_scl_h
    yh = y * lim_scl_h
    zh = z * lim_scl_h
    
    xl = x * lim_scl_l
    yl = y * lim_scl_l
    zl = z * lim_scl_l
    
        
    nom_l = func(a,b,c,xl,yl,zl)  * func(-a,-b,c,xl,yl,zl) * func(a,-b,-c,xl,yl,zl) * func(-a,b,-c,xl,yl,zl)
    nom_h = func(a,b,c,xh,yh,zh)  * func(-a,-b,c,xh,yh,zh) * func(a,-b,-c,xh,yh,zh) * func(-a,b,-c,xh,yh,zh)
    
    denom_l = func(a,-b,c,xl,yl,zl) * func(-a,b,c,xl,yl,zl)  * func(a,b,-c,xl,yl,zl)  * func(-a,-b,-c,xl,yl,zl)
    denom_h = func(a,-b,c,xh,yh,zh) * func(-a,b,c,xh,yh,zh)  * func(a,b,-c,xh,yh,zh)  * func(-a,-b,-c,xh,yh,zh)
    
    nom = 0.5 * ( nom_l + nom_h )
    
    denom = 0.5 * ( denom_l + denom_h )
    
    
    getF_limit = nom / denom

end function getF_limit

!:: Dims is defined as the dimensions (a/2.,b/2.,c/2.) for each prism
!:: As the code below assumes that dims is center +/- dims, we have to divide dims by two


!::function f from eq. 11
function f_3D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: f_3D, f_3D_l, f_3D_h, xh,xl
real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999

!::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
!::So we need to check for this and return a large number to ensure numerical stability

if ( a/2. - x .eq. 0 ) then
    !f_3D = sign(1e10, (b/2. - y) * (c/2. - z))
    xh = x * lim_scl_h
    xl = x * lim_scl_l
    
    f_3D_h = (b/2. - y) * (c/2. - z) / ( (a/2. - xh) * sqrt( (a/2. - xh)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
    f_3D_l = (b/2. - y) * (c/2. - z) / ( (a/2. - xl) * sqrt( (a/2. - xl)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
    
    f_3D = 0.5 * ( f_3D_h + f_3D_l )
    
else    
    f_3D = (b/2. - y) * (c/2. - z) / ( (a/2. - x) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
endif



return
end function f_3D


!::function g from eq. 11
function g_3D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: g_3D, g_3D_l, g_3D_h, yh,yl
real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999

!::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
!::So we need to check for this and return a large number to ensure numerical stability
if ( b/2. - y .eq. 0 ) then
    yh = y * lim_scl_h
    yl = y * lim_scl_l
    
    g_3D_h = (a/2. - x) * (c/2. - z) / ( (b/2. - yh) * sqrt( (a/2. - x)**2 + (b/2. - yh)**2 + (c/2. - z)**2 ) )
    g_3D_l = (a/2. - x) * (c/2. - z) / ( (b/2. - yl) * sqrt( (a/2. - x)**2 + (b/2. - yl)**2 + (c/2. - z)**2 ) )
    
    g_3D = 0.5 * ( g_3D_h + g_3D_l )
else    
    g_3D = (a/2. - x) * (c/2. - z) / ( (b/2. - y) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
endif
return
end function g_3D


!::function h from eq. 11
function h_3D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: h_3D,h_3D_l,h_3D_h,zh,zl
real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999
!::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
!::So we need to check for this and return a large number to ensure numerical stability
if ( c/2. - z .eq. 0 ) then
    zh = z * lim_scl_h
    zl = z * lim_scl_l
    
    h_3D_h = (a/2. - x) * (b/2. - y) / ( (c/2. - zh) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - zh)**2 ) )
    h_3D_l = (a/2. - x) * (b/2. - y) / ( (c/2. - zl) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - zl)**2 ) )
    
    h_3D = 0.5 * ( h_3D_l + h_3D_h )
    
else    
    h_3D = (a/2. - x) * (b/2. - y) / ( (c/2. - z) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
endif
return
end function h_3D


!::function F from eq. 11
function FF_3D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: FF_3D

FF_3D = (c/2. - z) + sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 )

return
end function FF_3D


!::function G from eq. 11
function GG_3D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: GG_3D

GG_3D = (a/2. - x) + sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 )

return
end function GG_3D


!::function H from eq. 11
function HH_3D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: HH_3D

HH_3D = (b/2. - y) + sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 )

return
end function HH_3D



!::Calculates N from the analytical expression in 2D
!::Given the dimensions of the prism (a,b,c) and the distance vector
!::to it
!:: In 2D, the off-diagonal elements are given by 
!:: Matlab: syms a b c x y z; pretty(limit(((-c-z)+sqrt((a-x)^2+(-b-y)^2+(-c-z)^2))*((-c-z)+sqrt((-a-x)^2+(b-y)^2+(-c-z)^2))/(((-c-z)+sqrt((a-x)^2+(b-y)^2+(-c-z)^2))*((-c-z)+sqrt((-a-x)^2+(-b-y)^2+(-c-z)^2))),c,inf))

subroutine getN_2D( a, b, c, x, y, z, N_out )
real,intent(in) :: a,b,c,x,y,z
real,intent(out),dimension(3,3) :: N_out
real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999
real :: xl,xh,yl,yh,nom_l,nom_h,denom_l,denom_h,lim,nom,denom
!::Diagonal elements
N_out(1,1) = 1./(4.*pi) * ( atan(f_2D(a,b,c,x,y,z))   + atan(f_2D(a,b,c,-x,y,z))  + atan(f_2D(a,b,c,x,-y,z)) + &
                            atan(f_2D(a,b,c,x,y,-z))  + atan(f_2D(a,b,c,-x,-y,z)) + atan(f_2D(a,b,c,x,-y,-z)) + &
                            atan(f_2D(a,b,c,-x,y,-z)) + atan(f_2D(a,b,c,-x,-y,-z)) )



N_out(2,2) = 1./(4.*pi) * ( atan(g_2D(a,b,c,x,y,z))   + atan(g_2D(a,b,c,-x,y,z))  + atan(g_2D(a,b,c,x,-y,z)) + &
                            atan(g_2D(a,b,c,x,y,-z))  + atan(g_2D(a,b,c,-x,-y,z)) + atan(g_2D(a,b,c,x,-y,-z)) + &
                            atan(g_2D(a,b,c,-x,y,-z)) + atan(g_2D(a,b,c,-x,-y,-z)) )
                            
                            

N_out(3,3) = 0                           


nom = (((a/2. + x)**2)/2. + ((b/2. - y)**2)/2.) * (((a/2. - x)**2)/2. + ((b/2. + y)**2)/2.)
denom = ((((a/2. - x)**2)/2. + ((b/2. - y)**2)/2.) * (((a/2. + x)**2)/2. + ((b/2. + y)**2)/2.))

if ( denom .eq. 0 .or. nom .eq. 0 ) then
    !Find the limit
    xl = x * lim_scl_l
    xh = x * lim_scl_h
    
    yl = y * lim_scl_l
    yh = y * lim_scl_h
    
    nom_l = (((a/2. + xl)**2)/2. + ((b/2. - yl)**2)/2.) * (((a/2. - xl)**2)/2. + ((b/2. + yl)**2)/2.)
    nom_h = (((a/2. + xh)**2)/2. + ((b/2. - yh)**2)/2.) * (((a/2. - xh)**2)/2. + ((b/2. + yh)**2)/2.)
    
    
    denom_l = ((((a/2. - xl)**2)/2. + ((b/2. - yl)**2)/2.) * (((a/2. + xl)**2)/2. + ((b/2. + yl)**2)/2.))
    denom_h = ((((a/2. - xh)**2)/2. + ((b/2. - yh)**2)/2.) * (((a/2. + xh)**2)/2. + ((b/2. + yh)**2)/2.))
    
    nom = 0.5 * ( nom_l + nom_h )
    denom = 0.5 * ( denom_l + denom_h )
    
    lim = nom / denom
    
    N_out(1,2) = -1./(4.*pi) * log( lim )        
else
    N_out(1,2) = -1./(4.*pi) * log( nom / denom )
endif


!::Off-diagonal elements
!N_out(1,2) = -1./(4.*pi) * log( (((a/2. + x)**2)/2. + ((b/2. - y)**2)/2.) * (((a/2. - x)**2)/2. + ((b/2. + y)**2)/2.) / &
!                               ((((a/2. - x)**2)/2. + ((b/2. - y)**2)/2.) * (((a/2. + x)**2)/2. + ((b/2. + y)**2)/2.)) )
                                
N_out(2,1) = N_out(1,2)

N_out(2,3) = 0
                                
N_out(3,2) = N_out(2,3)

N_out(1,3) = 0

N_out(3,1) = N_out(1,3)


end subroutine getN_2D

!:: Dims is defined as the dimensions (a/2.,b/2.,c/2.) for each prism
!:: As the code below assumes that dims is center +/- dims, we have to divide dims by two

!::function f from eq. 11
function f_2D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: f_2D,f_2D_l,f_2D_h,xl,xh
real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999
!::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
!::So we need to check for this and return a large number to ensure numerical stability
if ( a/2. - x .eq. 0 ) then
    xl = x * lim_scl_l
    xh = x * lim_scl_h
    
    f_2D_l = (b/2. - y) / (a/2. - xl)
    f_2D_h = (b/2. - y) / (a/2. - xh)
    
    f_2D = 0.5 * ( f_2D_l + f_2D_h )
else    
    f_2D = (b/2. - y) / (a/2. - x)
endif
return
end function f_2D


!::function g from eq. 11
function g_2D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: g_2D,g_2D_l,g_2D_h,yl,yh
real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999
!::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
!::So we need to check for this and return a large number to ensure numerical stability
if ( b/2. - y .eq. 0 ) then
    yl = y * lim_scl_l
    yh = y * lim_scl_h
    
    g_2D_l = (a/2. - x) / (b/2. - yl)
    g_2D_h = (a/2. - x) / (b/2. - yh)
    
    g_2D = 0.5 * ( g_2D_l + g_2D_h )
else    
    g_2D = (a/2. - x) / (b/2. - y)
endif
return
end function g_2D

   
END MODULE MagPrism_CALL
