module TileTriangle
    use SPECIALFUNCTIONS
    implicit none
    
    real,parameter :: pi=3.14159265359
    
    private :: pi
    
    contains
    
    
    
    !---------------------------------------------------------------------------    
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Returns the local-coordinate system demagnetization tensor field at the point r for a triangle
    !> defined through the vertices in v.
    !> @param[in] v is [3,4] where v(:,1) is the first vertex, v(:,2) and v(:,3) the second and third vertices.
    !> v(:,4) is the fourth vertex that defines the orientation of the normal to the triangle face such that v(:,4) is behind
    !> the surface with respect to the normal. 
    !> @param[in] r is (3,1) the point at which the field is to be evaluated
    !> @param[inout] N the 3x3 tensor which only has non-zero components at N(:,3) since it is defined in the local coordinate system
    !> @param[input] P the change-of-basis matrix from the local to the global coordinate system
    !---------------------------------------------------------------------------        
    subroutine getN_Triangle( v_in, r, N, P )
    ! Calculates and returns the demag tensor, N, given the four vertices,
    ! v (size 3,4) and at the positions r (size 3,1)) from a triangular face
    ! defined by vertices 1-3. The last vertex defines the orientation of the normal to the triangular
    ! surface such that the normal points away from this vertex.
    ! the array v is the four vertices organized as column vectors, 
    ! i.e. v(:,1) is the first vertex etc.
    real,dimension(3,4),intent(in) :: v_in         !> vertices
    real,dimension(3),intent(in) :: r           !> evaluation point
    real,dimension(3,3),intent(inout) :: N,P    !> output tensor and change-of-basis matrix
    
    real,parameter :: numErr = 1e-25,threshold = 1e-20 !> Numerical tolerance error and threshold for defining if the point of interest is too close to x,y or z = 0     
    real :: d1,d2    
    real,dimension(3) :: angles, angles_srt,tmp !> Angles of the triangle at each vertex
    integer,dimension(3) :: srt_ind             !> Indices of the sorted array
    real,dimension(3) :: e1,e2,e3,D,r_t         !> Unit vectors and translational point and transformed point of interest to the local system
    real,dimension(3,4) :: v                    !> copy of the vertices
    real,dimension(3,3) :: Pinv,vp              !> Inverse change-of-basis, i.e. from global to local system, local coordinates of the vertices
    integer :: sgn                              !> sign of a given coordinate
    real,dimension(3,3) :: N_loc                !> Local demag tensor
    integer :: i                                !> Loop variable
    
    v = v_in
    
    !test if the three vertices are collinear. If so, return without further calculation
    d1 = norm2( v(:,1)-v(:,2) )
    d2 = norm2( v(:,2)-v(:,3) )
    if ( abs( norm2(v(:,1)-v(:,3)) - (d1+d2)) .lt. numErr ) then
        write(*,*) 'Err. The vertices are collinear to within 1e-15'
        stop
    endif
    
    

    !ensure the vertices are ordered such that the largest angle is at the
    !middle vertex    
    angles(1) = acos( dot_product(v(:,1)-v(:,2),v(:,1)-v(:,3)) / (norm2(v(:,1)-v(:,2))*norm2(v(:,1)-v(:,3))) )
    angles(2) = acos( dot_product(v(:,2)-v(:,1),v(:,2)-v(:,3)) / (norm2(v(:,2)-v(:,1))*norm2(v(:,2)-v(:,3))) )
    angles(3) = acos( dot_product(v(:,3)-v(:,2),v(:,3)-v(:,1)) / (norm2(v(:,3)-v(:,2))*norm2(v(:,3)-v(:,1))) )
    
    !order the vertices after angle where the largest angle is in the middle
    call simple_sort( angles, angles_srt, srt_ind )
    v(:,1:3) = v(:,srt_ind)
    tmp = v(:,3)
    v(:,3) = v(:,2)
    v(:,2) = tmp

    angles = angles_srt
    tmp(1) = angles(3)
    angles(3) = angles(2)
    angles(2) = tmp(1)
    
    !find the transformation matrix from the local to the global system
    !first unit vector
    e1 = v(:,1)-v(:,3)
    e1 = e1/norm2(e1)
    !third unit vector
    e3 = cross(e1, v(:,2)-v(:,3) )
    e3 = e3 / norm2(e3)

    !ensure the fourth vertex is not in the triangular plane
    if ( abs( dot_product( (v(:,4)-v(:,1)), e3 ) ) .le. numErr ) then
        write(*,*) 'Vertices are all in the same plane'
        stop
    endif

    !check which way the surface normal should point
    if ( dot_product( e3, v(:,4)-v(:,1) ) .gt. 0. ) then
        !switch around v1 and v3
        tmp = v(:,1)
        v(:,1) = v(:,3)
        v(:,3) = tmp
        
        tmp(1) = angles(1)
        angles(1) = angles(3)
        angles(3) = tmp(1)
        
        !find the two first unit vectors again
        e1 = v(:,1)-v(:,3)
        e1 = e1/norm2(e1)
        !third unit vector
        e3 = cross( e1, v(:,2)-v(:,3) )
        e3 = e3 / norm2(e3)
    endif

    !second unit vector
    e2 = cross( e3, e1 )

    !p-matrix
    P(:,1) = e1
    P(:,2) = e2
    P(:,3) = e3
    !the inverse is the transpose as P is orthogonal
    Pinv = transpose(P)

    !find the point D in order to translate the triangle to coincide with the
    !Origin
    D = cos(angles(3)) * norm2( v(:,2)-v(:,3) ) * e1 + v(:,3)

    !transformed coordinates (from global to local)    
    do i=1,3
        vp(:,i) = matmul( Pinv, ( v(:,i)-D ) )
    enddo
    
    !transform the point of interest to the local system
    r_t = matmul( Pinv, ( r - D ) )

    do i=1,3
        if ( abs( r_t(i) ) .lt. threshold ) then
            r_t(i) = sign(threshold,r_t(i))            
        endif
    enddo


    !Local demag tensor
    N_loc(:,:) = 0.
    
    !The minus sign reflects the change in integration limits on the integral over the triangle in the second quadrant
    N_loc(1,3) = Nxz( r_t, vp(1,1), vp(2,2) ) - Nxz( r_t, vp(1,3), vp(2,2) )
    
    N_loc(2,3) = Nyz( r_t, vp(1,1), vp(2,2) ) - Nyz( r_t, vp(1,3), vp(2,2) )
    
    N_loc(3,3) = Nzz( r_t, vp(1,1), vp(2,2) ) - Nzz( r_t, vp(1,3), vp(2,2) )


    !Apply change of basis
    N = matmul( matmul( P, N_loc ), Pinv )


    end subroutine getN_Triangle

    !---------------------------------------------------------------------------    
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Nxz component of the local tensor of the right triangle
    !> @param[in] r the point of interest, size 3,1
    !> @param[in] l the length along the x'-direction of the triangle
    !> @param[in] h the length along the y'-direction of the triangle
    !---------------------------------------------------------------------------        
    function Nxz( r, l, h )
    real :: Nxz
    real,dimension(3),intent(in) :: r
    real,intent(in) :: l,h

    
    ! Returns the Nxz tensor component in the local coordinate system
        Nxz = -1./(4.*pi) * ( F_Nxz(r,h,l,h) - F_Nxz(r,0.,l,h) - ( G_Nxz(r,h) - G_Nxz(r,0.) ) )
        
    end function Nxz

    function F_Nxz( r, yp, l, h )
    real :: F_Nxz
    real,dimension(3),intent(in) :: r
    real,intent(in) :: l,h,yp
    
    !RUBI solution (Int):
        F_Nxz = h / sqrt( h**2 + l**2 ) * atanh( (l**2 - l*r(1) + h*r(2) - h*yp*(1+l**2/h**2)) / &
            ( sqrt(h**2+l**2)*sqrt( l**2 - 2*l*r(1) + r(1)**2 + r(2)**2 - 2*yp*(l**2-l*r(1)+h*r(2))/h + yp**2*(1.+l**2/h**2) + r(3)**2) ) )

    
    end function F_Nxz

    function G_Nxz( r, yp )
    real :: G_nxz
    real,dimension(3),intent(in) :: r
    real,intent(in) :: yp
        G_Nxz = atanh( ( r(2)-yp )/ (sqrt(r(1)**2+(r(2)-yp)**2+r(3)**2)) )
    end function G_Nxz

    !---------------------------------------------------------------------------    
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Nyz component of the local tensor of the right triangle
    !> @param[in] r the point of interest, size 3,1
    !> @param[in] l the length along the x'-direction of the triangle
    !> @param[in] h the length along the y'-direction of the triangle
    !---------------------------------------------------------------------------    
    function Nyz( r, l, h )
    real :: Nyz
    real,dimension(3),intent(in) :: r
    real,intent(in) :: l,h
       ! returns the Nyz component of the tensor in the local coordinate system
       Nyz = -1./(4.*pi) * ( K_Nyz(r,l,l,h) - K_Nyz(r,0.,l,h) - ( L_Nyz(r,l) - L_Nyz(r,0.) ) )
      
    end function Nyz
   
    function K_Nyz( r, xp, l, h )
    real :: K_Nyz
    real,dimension(3),intent(in) :: r
    real,intent(in) :: l,h,xp
        K_Nyz = l/sqrt(h**2+l**2) * atanh( (h**2-h*r(2)+l*r(1)-l*xp*(1+h**2/l**2) ) / &
            (sqrt(h**2+l**2)*sqrt(h**2-2*h*r(2)+r(1)**2 + r(2)**2 - 2*xp*(h**2+l*r(1)-h*r(2))/l +xp**2*(1+h**2/l**2) + r(3)**2)) )

    end function K_Nyz


    function L_Nyz( r, xp )
    real :: L_Nyz
    real,dimension(3),intent(in) :: r
    real,intent(in) :: xp
        L_Nyz = atanh( (r(1) - xp) / (sqrt((r(1)-xp)**2+r(2)**2+r(3)**2)) )
    end function L_Nyz

    !---------------------------------------------------------------------------    
    !> @author Kaspar K. Nielsen, kasparkn@gmail.com, DTU, 2019
    !> @brief
    !> Nzz component of the local tensor of the right triangle
    !> @param[in] r the point of interest, size 3,1
    !> @param[in] l the length along the x'-direction of the triangle
    !> @param[in] h the length along the y'-direction of the triangle
    !---------------------------------------------------------------------------  
    function Nzz( r, l, h )
    real :: Nzz
    real,dimension(3),intent(in) :: r
    real,intent(in) :: l,h
        Nzz = -1./(4.*pi) * ( P_Nzz( r, l, l, h ) - P_Nzz( r, 0., l, h ) - ( Q_Nzz(r,l) - Q_Nzz(r,0.) ) )        
    end function Nzz

    function P_Nzz( r, xp, l, h )
    real :: P_Nzz
    real,dimension(3),intent(in) :: r
    real,intent(in) :: l,h,xp
    
        P_Nzz = atan( ( r(1)*(h-r(2)) - xp*(h*(1-r(1)/l)-r(2)) - h*(r(1)**2+r(3)**2)/l ) / &
            (r(3)*sqrt(h**2+r(1)**2 + xp**2*(1+h**2/l**2) - 2*h*r(2) + r(2)**2 - 2*xp*(h**2+l*r(1)-h*r(2))/l + r(3)**2 )) )
    
    end function P_Nzz

    function Q_Nzz( r, xp )
    real :: Q_Nzz
    real,intent(in) :: xp
    real,dimension(3),intent(in) :: r
        Q_Nzz = -atan( (r(1)-xp)*r(2) / ( r(3) * sqrt( (r(1)-xp)**2 + r(2)**2 + r(3)**2 ) ) )
        
    end function Q_Nzz

    
    
end module TileTriangle
    