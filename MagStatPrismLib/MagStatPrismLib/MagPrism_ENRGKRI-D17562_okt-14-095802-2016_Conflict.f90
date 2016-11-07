MODULE MagPrism_CALL
use PARAMETERS_CALL
use UTIL_CALL
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
subroutine calcH( pos, dims, dir, magType, Hext,T, n, stateFcn, stateFcnIndices, m, solPts, l, spacedim, Hout, Mout, nIte, rotAngles )
real,intent(in),dimension(n,3) :: pos, dims,dir,Hext,rotAngles
integer,intent(in),dimension(n) :: magType
real,intent(in),dimension(n) :: T
type(stateFunction),intent(in),dimension(m) :: stateFcn
integer,intent(in),dimension(n) :: stateFcnIndices
integer,intent(in) :: n,m,l,spacedim
real,intent(in),dimension(l,3) :: solPts
real,intent(out),dimension(:,:),allocatable :: Hout,Mout
integer,intent(inout) :: nIte

real :: lambda,maxRelDiff,lambda_old
real,dimension(4) :: maxRelDiffArr
real,dimension(:,:),allocatable :: Hdem,H_old
integer :: i,j,convCnt,lambdaCnt
integer,parameter :: maxIte=1000
real,parameter :: tol = 0.001
logical :: lCh
real,dimension(:,:),allocatable :: tmp
real,dimension(:,:,:),allocatable :: rotMat,rotMatInv

lambda = 1.0
lambda_old = lambda
convCnt = 1
lambdaCnt = 1
call writeDebugString( 'start' )


!::check that maxval(stateFcnIndices) <= m
if ( maxval(stateFcnIndices) .le. m ) then
    !::Allocate the solution
    allocate(Hout(n+l,3),Mout(n,3),Hdem(n,3),H_old(n,3),tmp(n,3))
    Hout(:,:) = 0
    Hout(1:n,:) = Hext
    Mout(:,:) = 0
    H_old(:,:) = 0
    
    !::Get the rotation matrix and its inverse
    call findRotationMatrices( rotAngles, n, rotMat, rotMatInv)
    
    !::First we find the self-consistent solution, i.e. the magnetization (magnitude and direction) in each prism    
    do i=0,maxIte
        nIte = i
        H_old = Hout(1:n,:)
        !::Get the magnetization in each node
        call getM(Hout(1:n,:), T, Mout, dir, magType, n, stateFcn, stateFcnIndices, m)
        !call writeDebugString( 'index', i )        
        !call writeDebugString( 'minM ',i, minval(Mout) )
        !call writeDebugString( 'maxM ',i, maxval(Mout) )
        !::Get the demagnetizing field
        call getHdem(dims,pos,Mout,n,spacedim,Hdem, rotMat, rotMatInv)
        !call writeDebugStringArr2D('Hdem ', Hdem )
        !::Update the solution
        Hout(1:n,:) = Hout(1:n,:) + lambda * ( (Hext + Hdem) - Hout(1:n,:) )
        !call writeDebugString( 'index', i )        
        !call writeDebugString( 'minH ',i, minval(Hout(1:n,:)) )
        !call writeDebugString( 'maxH ',i, maxval(Hout(1:n,:)) )
        !::Find the relative change in the solution
        tmp = -1
        do j=1,n
            where( H_old(j,:) .ne. 0 )
                tmp(j,:) = abs( (Hout(j,:) - H_old(j,:))/H_old(j,:) )
            endwhere
        enddo
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
        
        
    enddo
    
    
    !::Having obtained the solution we compute the magnetic field at each of the required points
    call getSolution( solPts, l, Mout, dims, pos, n, spacedim, Hout(n+1:n+l,:), rotMat, rotMatInv )
    
endif

end subroutine calcH

!::Returns the solution at the required points
subroutine getSolution( solPts, l, M, dims, pos, n, spacedim, Hout, rotMat, rotMatInv )
real,intent(in),dimension(l,3) :: solPts
real,intent(inout),dimension(l,3) :: Hout
real,intent(in),dimension(n,3) :: M, dims, pos
real,intent(in),dimension(n,3,3) :: rotMat,rotMatInv
integer,intent(in) :: l,n, spacedim

real,dimension(3) :: diffPos, dotProd
real,dimension(3,3) :: Nout
integer :: i,j


!::Loop over each point where the solution is required
do i=1,l
    !::Loop over each prism
    do j=1,n
        diffPos = pos(j,:) - solPts(i,:)
        !::Rotate the position vector to accomodate rotation of the j'th prism
        diffPos = matmul( rotMat(j,:,:), diffPos )
        if ( spacedim .eq. 3 ) then
           call getN_3D(dims(j,1),dims(j,2),dims(j,3), diffPos(1), diffPos(2), diffPos(3), Nout )
        endif
        if ( spacedim .eq. 2 ) then
           call getN_2D(dims(j,1),dims(j,2),dims(j,3), diffPos(1), diffPos(2), diffPos(3), Nout )
        endif
        
        call getDotProd( Nout, M(j,:), dotProd )
        
        !::Rotate the solution back
        dotProd = matmul( rotMatInv(j,:,:), dotProd )
        
        Hout(i,:) = Hout(i,:) - dotProd
    enddo
    
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

    lambda = lambda * 0.9
    lCh = .true.
endif


end subroutine updateLambda

!::Returns the magnitude of the magnetization in each prism
subroutine getM( H, T, Mout, dir, magType, n, stateFcn, stateFcnIndices, m)
real,intent(in),dimension(n,3) :: H,dir
real,intent(in),dimension(n) :: T
real,intent(inout),dimension(n,3) :: Mout
integer,intent(in),dimension(n) :: magType
integer,intent(in) :: n,m
type(stateFunction),intent(in),dimension(m) :: stateFcn
integer,intent(in),dimension(n) :: stateFcnIndices

real,dimension(:),allocatable :: Hnorm

real :: Mtmp
integer :: i

allocate(Hnorm(n))

Hnorm = sqrt( H(:,1)**2 + H(:,2)**2 + H(:,3)**2 )

!call writeDebugString( 'Hnorm',1,minval(Hnorm) )
!call writeDebugString( 'Hnorm',1,maxval(Hnorm) )
!::For each prism
do i=1,n
    call getBilinInterp( stateFcn(stateFcnIndices(i))%M, stateFcn(stateFcnIndices(i))%T, &
                         stateFcn(stateFcnIndices(i))%H, stateFcn(stateFcnIndices(i))%nT, &
                         stateFcn(stateFcnIndices(i))%nH, T(i), Hnorm(i), Mtmp )    
    !call writeDebugString('Hnorm', i, Hnorm(i) )
    !call writeDebugString('Mtmp ', i, Mtmp )
    if ( magType(i) .eq. magTypeSoft ) then
        !::Then the magnetization is along the direction of the local field
        if ( Hnorm(i) .ne. 0 ) then
            Mout(i,:) = Mtmp * H(i,:) / Hnorm(i)      
        else
            Mout(i,1) = Mtmp
            Mout(i,2:3) = 0        
            !call writeDebugString('Hnorm', i, Hnorm(i) )
        endif        
    elseif ( magType(i) .eq. magTypeHard ) then
        !::Then the magnetization is along the direction pre-specified
        Mout(i,:) = Mtmp * dir(i,:)
    endif
    
    
enddo


end subroutine getM

subroutine getHdem(dims,pos,M,n,spacedim,Hdem,rotMat,rotMatInv)
real,intent(in),dimension(n,3) :: dims,pos,M
real,intent(in),dimension(n,3,3) :: rotMat,rotMatInv
integer,intent(in) :: n, spacedim
real,intent(inout),dimension(n,3) :: Hdem

real,dimension(3,3) :: N_out
real,dimension(3) :: diffPos,dotProd
integer :: i,j

Hdem(:,:) = 0
call writeDebugStringArr2D('dims ', dims )
!$OMP PARALLEL DO PRIVATE(i,j,diffPos,N_out,dotProd)
do i=1,n
    !::For each prism find the contribution from all other prisms
    do j=1,n        
        diffPos(1) = pos(i,1)-pos(j,1)
        diffPos(2) = pos(i,2)-pos(j,2)
        diffPos(3) = pos(i,3)-pos(j,3)
        
        !::Rotate the postion vector
        diffPos = matmul( rotMat(j,:,:), diffPos )
        
        call writeDebugStringArr1D('PosAft',diffPos)
        call writeDebugStringArr1D('dims ', dims(j,:) )
        
        !::Find the demag. tensor
        if ( spacedim .eq. 3 ) then
           call getN_3D(dims(j,1),dims(j,2),dims(j,3),diffPos(1),diffPos(2),diffPos(3), N_out )
        endif
        if ( spacedim .eq. 2 ) then
           call getN_2D(dims(j,1),dims(j,2),dims(j,3),diffPos(1),diffPos(2),diffPos(3), N_out )
        endif
        call writeDebugString('ind i', i)
        call writeDebugString('ind j', j)
        !call writeDebugStringArr2D( 'N_out ', N_out )
        !::Compute the dot product between N and M
        call getDotProd( N_out, M(j,:), dotProd )
        
        !::Rotate the solution back
        dotProd = matmul( rotMatInv(j,:,:), dotProd )        
        
        !::Add to the solution
        Hdem(i,:) = Hdem(i,:) - dotProd
        
        
                
    enddo        
enddo
!$OMP END PARALLEL DO
!call writeDebugString('Hdem ',1,minval(Hdem))
!call writeDebugString('Hdem ',1,maxval(Hdem))

end subroutine getHdem

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
    call getRotX( rotAngles(i,1), RotX )
    call getRotY( rotAngles(i,2), RotY )
    call getRotZ( rotAngles(i,3), RotZ )
    
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
!::Returns the rotation matrix for a rotation of angle radians about the y-axis
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
subroutine getN_3D( a, b, c, x, y, z, N_out )
real,intent(in) :: a,b,c,x,y,z
real,intent(out),dimension(3,3) :: N_out
real :: nom,denom,nom_l,nom_h,denom_l,denom_h
!real,parameter :: zero_limit=1.006258
real,parameter :: lim_scl_h=1.0001,lim_scl_l=0.9999
real :: xl,xh,yl,yh,zl,zh,al,ah,bl,bh,cl,ch
real,dimension(6) :: tmp

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
    xh = x * lim_scl_h
    yh = y * lim_scl_h
    zh = z * lim_scl_h
    
    ah = a * lim_scl_h
    bh = b * lim_scl_h
    ch = c * lim_scl_h
    
    xl = x * lim_scl_l
    yl = y * lim_scl_l
    zl = z * lim_scl_l
    
    al = a * lim_scl_l
    bl = b * lim_scl_l
    cl = c * lim_scl_l
    
    nom_l = FF_3D(al,bl,cl,xl,yl,zl)  * FF_3D(-al,-bl,cl,xl,yl,zl) * FF_3D(al,-bl,-cl,xl,yl,zl) * FF_3D(-al,bl,-cl,xl,yl,zl)
    nom_h = FF_3D(ah,bh,ch,xh,yh,zh)  * FF_3D(-ah,-bh,ch,xh,yh,zh) * FF_3D(ah,-bh,-ch,xh,yh,zh) * FF_3D(-ah,bh,-ch,xh,yh,zh)
    
    denom_l = FF_3D(al,-bl,cl,xl,yl,zl) * FF_3D(-al,bl,cl,xl,yl,zl)  * FF_3D(al,bl,-cl,xl,yl,zl)  * FF_3D(-al,-bl,-cl,xl,yl,zl)
    denom_h = FF_3D(ah,-bh,ch,xh,yh,zh) * FF_3D(-ah,bh,ch,xh,yh,zh)  * FF_3D(ah,bh,-ch,xh,yh,zh)  * FF_3D(-ah,-bh,-ch,xh,yh,zh)
    
    nom = 0.5 * ( nom_l + nom_h )
    
    denom = 0.5 * ( denom_l + denom_h )
    
    tmp(1) = nom_l
    tmp(2) = nom_h
    tmp(3) = denom_l
    tmp(4) = denom_h
    tmp(5) = nom
    tmp(1) = denom
    call writeDebugStringArr1D('tmpF ', tmp )
endif

N_out(1,2) = -1./(4.*pi) * log( nom / denom )

N_out(2,1) = N_out(1,2)

nom = GG_3D(a,b,c,x,y,z)  * GG_3D(-a,-b,c,x,y,z) * GG_3D(a,-b,-c,x,y,z) * GG_3D(-a,b,-c,x,y,z)
denom = GG_3D(a,-b,c,x,y,z) * GG_3D(-a,b,c,x,y,z)  * GG_3D(a,b,-c,x,y,z)  * GG_3D(-a,-b,-c,x,y,z)

if ( denom .eq. 0 .or. nom .eq. 0 ) then
!Find the limit
    xh = x * lim_scl_h
    yh = y * lim_scl_h
    zh = z * lim_scl_h
    
    ah = a * lim_scl_h
    bh = b * lim_scl_h
    ch = c * lim_scl_h
    
    xl = x * lim_scl_l
    yl = y * lim_scl_l
    zl = z * lim_scl_l
    
    al = a * lim_scl_l
    bl = b * lim_scl_l
    cl = c * lim_scl_l
    
    nom_l = GG_3D(al,bl,cl,xl,yl,zl)  * GG_3D(-al,-bl,cl,xl,yl,zl) * GG_3D(al,-bl,-cl,xl,yl,zl) * GG_3D(-al,bl,-cl,xl,yl,zl)
    nom_h = GG_3D(ah,bh,ch,xh,yh,zh)  * GG_3D(-ah,-bh,ch,xh,yh,zh) * GG_3D(ah,-bh,-ch,xh,yh,zh) * GG_3D(-ah,bh,-ch,xh,yh,zh)
    
    denom_l = GG_3D(al,-bl,cl,xl,yl,zl) * GG_3D(-al,bl,cl,xl,yl,zl)  * GG_3D(al,bl,-cl,xl,yl,zl)  * GG_3D(-al,-bl,-cl,xl,yl,zl)
    denom_h = GG_3D(ah,-bh,ch,xh,yh,zh) * GG_3D(-ah,bh,ch,xh,yh,zh)  * GG_3D(ah,bh,-ch,xh,yh,zh)  * GG_3D(-ah,-bh,-ch,xh,yh,zh)
    
    nom = 0.5 * ( nom_l + nom_h )
    
    denom = 0.5 * ( denom_l + denom_h )
    
    tmp(1) = nom_l
    tmp(2) = nom_h
    tmp(3) = denom_l
    tmp(4) = denom_h
    tmp(5) = nom
    tmp(1) = denom
    call writeDebugStringArr1D('tmpG ', tmp )
        
endif
N_out(2,3) = -1./(4.*pi) * log( nom/denom )

N_out(3,2) = N_out(2,3)

nom = HH_3D(a,b,c,x,y,z)  * HH_3D(-a,-b,c,x,y,z) * HH_3D(a,-b,-c,x,y,z) * HH_3D(-a,b,-c,x,y,z)
denom = HH_3D(a,-b,c,x,y,z) * HH_3D(-a,b,c,x,y,z)  * HH_3D(a,b,-c,x,y,z)  * HH_3D(-a,-b,-c,x,y,z)

if ( denom .eq. 0 .or. nom .eq. 0 ) then
    !Find the limit
    xh = x * lim_scl_h
    yh = y * lim_scl_h
    zh = z * lim_scl_h
    
    ah = a * lim_scl_h
    bh = b * lim_scl_h
    ch = c * lim_scl_h
    
    xl = x * lim_scl_l
    yl = y * lim_scl_l
    zl = z * lim_scl_l
    
    al = a * lim_scl_l
    bl = b * lim_scl_l
    cl = c * lim_scl_l
    
    nom_l = HH_3D(al,bl,cl,xl,yl,zl)  * HH_3D(-al,-bl,cl,xl,yl,zl) * HH_3D(al,-bl,-cl,xl,yl,zl) * HH_3D(-al,bl,-cl,xl,yl,zl)
    nom_h = HH_3D(ah,bh,ch,xh,yh,zh)  * HH_3D(-ah,-bh,ch,xh,yh,zh) * HH_3D(ah,-bh,-ch,xh,yh,zh) * HH_3D(-ah,bh,-ch,xh,yh,zh)
    
    denom_l = HH_3D(al,-bl,cl,xl,yl,zl) * HH_3D(-al,bl,cl,xl,yl,zl)  * HH_3D(al,bl,-cl,xl,yl,zl)  * HH_3D(-al,-bl,-cl,xl,yl,zl)
    denom_h = HH_3D(ah,-bh,ch,xh,yh,zh) * HH_3D(-ah,bh,ch,xh,yh,zh)  * HH_3D(ah,bh,-ch,xh,yh,zh)  * HH_3D(-ah,-bh,-ch,xh,yh,zh)
    
    nom = 0.5 * ( nom_l + nom_h )
    
    denom = 0.5 * ( denom_l + denom_h )
    
    tmp(1) = nom_l
    tmp(2) = nom_h
    tmp(3) = denom_l
    tmp(4) = denom_h
    tmp(5) = nom
    tmp(1) = denom
    call writeDebugStringArr1D('tmpH ', tmp )
    
endif

N_out(1,3) = -1./(4.*pi) * log( nom / denom )

N_out(3,1) = N_out(1,3)


!call writeDebugStringArr2D('Nout ', N_out )

end subroutine getN_3D

!:: Dims is defined as the dimensions (a/2.,b/2.,c/2.) for each prism
!:: As the code below assumes that dims is center +/- dims, we have to divide dims by two

!::function f from eq. 11
function f_3D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: f_3D

!::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
!::So we need to check for this and return a large number to ensure numerical stability

if ( a/2. - x .eq. 0 ) then
    f_3D = sign(1e10, (b/2. - y) * (c/2. - z))
else    
    f_3D = (b/2. - y) * (c/2. - z) / ( (a/2. - x) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
endif

return
end function f_3D


!::function g from eq. 11
function g_3D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: g_3D

!::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
!::So we need to check for this and return a large number to ensure numerical stability
if ( b/2. - y .eq. 0 ) then
    g_3D = sign(1e10,(a/2. - x) * (c/2. - z))
else    
    g_3D = (a/2. - x) * (c/2. - z) / ( (b/2. - y) * sqrt( (a/2. - x)**2 + (b/2. - y)**2 + (c/2. - z)**2 ) )
endif
return
end function g_3D


!::function h from eq. 11
function h_3D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: h_3D
!::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
!::So we need to check for this and return a large number to ensure numerical stability
if ( c/2. - z .eq. 0 ) then
    h_3D = sign(1e10,(a/2. - x) * (b/2. - y))
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

!::Diagonal elements
N_out(1,1) = 1./(4.*pi) * ( atan(f_2D(a,b,c,x,y,z))   + atan(f_2D(a,b,c,-x,y,z))  + atan(f_2D(a,b,c,x,-y,z)) + &
                            atan(f_2D(a,b,c,x,y,-z))  + atan(f_2D(a,b,c,-x,-y,z)) + atan(f_2D(a,b,c,x,-y,-z)) + &
                            atan(f_2D(a,b,c,-x,y,-z)) + atan(f_2D(a,b,c,-x,-y,-z)) )



N_out(2,2) = 1./(4.*pi) * ( atan(g_2D(a,b,c,x,y,z))   + atan(g_2D(a,b,c,-x,y,z))  + atan(g_2D(a,b,c,x,-y,z)) + &
                            atan(g_2D(a,b,c,x,y,-z))  + atan(g_2D(a,b,c,-x,-y,z)) + atan(g_2D(a,b,c,x,-y,-z)) + &
                            atan(g_2D(a,b,c,-x,y,-z)) + atan(g_2D(a,b,c,-x,-y,-z)) )
                            
                            

N_out(3,3) = 0                           
                            

!::Off-diagonal elements
N_out(1,2) = -1./(4.*pi) * log( (((a/2. + x)**2)/2. + ((b/2. - y)**2)/2.) * (((a/2. - x)**2)/2. + ((b/2. + y)**2)/2.) / &
                               ((((a/2. - x)**2)/2. + ((b/2. - y)**2)/2.) * (((a/2. + x)**2)/2. + ((b/2. + y)**2)/2.)) )
                                
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
real :: f_2D

!::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
!::So we need to check for this and return a large number to ensure numerical stability
if ( a/2. - x .eq. 0 ) then
    f_2D = sign(1e10,(b/2. - y))
else    
    f_2D = (b/2. - y) / (a/2. - x)
endif
return
end function f_2D


!::function g from eq. 11
function g_2D( a, b, c, x, y, z )
real,intent(in) :: a,b,c,x,y,z
real :: g_2D
!::Numerical problem: atan(inf) = pi/2 and when a/2 = x then f_3D returns NaN (inf). 
!::So we need to check for this and return a large number to ensure numerical stability
if ( b/2. - y .eq. 0 ) then
    g_2D = sign(1e10,(a/2. - x))
else    
    g_2D = (a/2. - x) / (b/2. - y)
endif
return
end function g_2D

   
END MODULE MagPrism_CALL
