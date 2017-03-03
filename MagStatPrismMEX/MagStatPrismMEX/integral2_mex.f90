#include "fintrf.h"  
      module model_glob
      USE MagPrism_Call
      USE Parameters_CALL
      USE NumInt
      implicit none
      
      type modelStruct
      integer,dimension(3) :: inds
      real*8 :: z0,r0
      integer :: n,spacedim,intInd
      real*8,dimension(:,:),allocatable :: dims,pos,M,rotAngles
      real*8,dimension(:,:,:),allocatable :: rotMat,rotMatInv
      integer*8 :: cylPlane
      integer*4,dimension(3) :: retForces
      endtype
      integer*8,parameter :: phi_z_plane=1,r_phi_plane=2
      type(modelStruct) :: model
      
      real*8,parameter :: mu0=4*pi*1e-7
            
       procedure (func), pointer :: f_tens_ptr => null()
       
       procedure (func), pointer :: f_tens_ptr_cyl1 => null()
       procedure (func), pointer :: f_tens_ptr_cyl2 => null()
       procedure (func), pointer :: f_tens_ptr_cyl3 => null()
      
      contains
      !::
      !::Integrates Maxwell's stress tensor over a cylindrical volume
      !::defines from rMax and zPos (zPos(1) = zMin and zPos(2) = zMax)
      !::eps_abs defines the minimum absolute error
      !::eps_rel defines the minimum relative error
      !::F is the resulting force
      !::ier is an integer describing the error state 
      !::of the integration (0 = no error)
      !:: abserr_tot is the total absolute error
      !::neval is the number of function evaluations used by the integrator
      subroutine cylIntegral( rMax, zPos, eps_abs, eps_rel, 
     +                       F, ier,abserr_tot, neval)
      real*8,intent(in) :: rMax,eps_abs,eps_rel
      real*8,dimension(2),intent(in) :: zPos
      real*8,dimension(3),intent(inout) :: F
      integer*8,intent(inout) :: ier,neval
      real*8,intent(inout) :: abserr_tot
      integer*8 :: j
      real*8 :: res,abserr
      real*8,parameter :: zero=0.0
       !::Reset input array
      F(:) = 0.
      
      !set to a large number
      abserr_tot = 1e9
      neval = 0
      !::Get the rotation matrices
      call findRotationMatrices( model%rotAngles, model%n, 
     +                            model%rotMat, model%rotMatInv)
      
     
       do j=3,3
          
      if ( j .eq. 1 ) then        
        !::The r-component of the force
        f_tens_ptr_cyl1 => F_int_diag_cyl
        f_tens_ptr_cyl2 => F_int_off_diag_cyl
        f_tens_ptr_cyl3 => F_int_off_diag_cyl
        model%intInd = 3
        model%inds(3) = 1
      elseif ( j .eq. 2 ) then
        !:: the phi-component of the force
        f_tens_ptr_cyl1 => F_int_off_diag_cyl
        f_tens_ptr_cyl2 => F_int_off_diag_cyl
        f_tens_ptr_cyl3 => F_int_off_diag_cyl        
        model%inds(3) = 2
        model%intInd = 1
      elseif ( j .eq. 3 ) then
        !:: the z-component of the force
        f_tens_ptr_cyl1 => F_int_off_diag_cyl
        f_tens_ptr_cyl2 => F_int_diag_cyl
        f_tens_ptr_cyl3 => F_int_diag_cyl
        model%intInd = 1
        model%inds(3) = 3
      endif    
          
      !::The phi-z plane
      model%cylPlane = phi_z_plane
      model%r0 = rMax
      !call integral2( f_tens_ptr_cyl1, zero, 2*pi, 
      !+                                zPos(1), zPos(2), 
      !+                                eps_abs, eps_rel,
      !+                                res, abserr, neval, ier )
      !    if ( abserr .lt. abserr_tot ) then
      !        abserr_tot = abserr
      !    endif      
      !    if ( ier .gt. 0 ) then
      !        exit
      !    endif
      !F(j) = res
      !::Lower r-phi plane
      model%cylPlane = r_phi_plane
      model%z0 = zPos(1)
      
      if ( j .eq. 2 ) then
          model%intInd = 3
      endif      
      call integral2( f_tens_ptr_cyl2, zero, rMax, 
     +                                zero, 2*pi, 
     +                                eps_abs, eps_rel,
     +                                res, abserr, neval, ier )
          if ( abserr .lt. abserr_tot ) then
              abserr_tot = abserr
          endif      
          if ( ier .gt. 0 ) then
              exit
          endif
      F(j) = F(j) - res
      
      !::Upper plane
      model%cylPlane = r_phi_plane
      model%z0 = zPos(2)
      call integral2( f_tens_ptr_cyl3, zero, rMax, 
     +                                zero, 2*pi, 
     +                                eps_abs, eps_rel,
     +                                res, abserr, neval, ier )
          if ( abserr .lt. abserr_tot ) then
              abserr_tot = abserr
          endif      
          if ( ier .gt. 0 ) then
              exit
          endif
      F(j) = F(j) + res
      enddo
      end subroutine cylIntegral
      
      !::Returns the integrand for the off-diagonal part of the Maxwell stress
      !::tensor dotted with the surface area
      function F_int_off_diag_cyl( v1, v2 )
      real*8 :: F_int_off_diag_cyl
      real*8,intent(in) :: v1,v2
      real*8 :: x,y,z,r,phi          
      real*8,dimension(3) :: solPts,B,Hsol
      
      Hsol(:) = 0
      
      call car2cyl( v1, v2, r, phi, x, y, z )
      
      !::setup the points where we wish to get the solution
      solPts(1) = x
      solPts(2) = y
      solPts(3) = z

      !::call the subroutine to find the solution at the specified points
      call getSolution( solPts, 1, model%M, model%dims, 
     +                  model%pos, model%n, model%spacedim,
     +                  Hsol, model%rotMat, model%rotMatInv )
      
      !::Convert to cylindrical coordiantes                                  
      call getCylVec( phi, Hsol, B )
    
       F_int_off_diag_cyl = 
     + 1. /mu0 * r * B(model%intInd) * B(model%inds(3))
    
       return

      end function F_int_off_diag_cyl
      
      !::Returns the integrand of the Maxwell stress tensor along the diagonal over
      !::the given surface area. 
      !::            
      function F_int_diag_cyl( v1, v2 )
      real*8 :: F_int_diag_cyl
      real*8,intent(in) :: v1,v2
      real*8 :: x,y,z,r,phi      
      real*8,dimension(3) :: solPts,B,Hsol
      real*8 :: Bnorm
            
      Hsol(:) = 0
       
      call car2cyl( v1, v2, r, phi, x, y, z )
      
      solPts(model%inds(1)) = x
      solPts(model%inds(2)) = y
      solPts(model%inds(3)) = z

    
      call getSolution( solPts, 1, model%M, model%dims, 
     +                  model%pos, model%n, model%spacedim,
     +                  Hsol, model%rotMat, model%rotMatInv )
	
      call getCylVec( phi, Hsol, B )
      Bnorm = sqrt( B(1)**2 + B(2)**2 + B(3)**2 )
                                 
      F_int_diag_cyl = 1./mu0 * 
     +                   ( B(model%inds(3))**2 - 0.5 * Bnorm**2 )
    
        
      end function F_int_diag_cyl
      !::
      !::Converts v1 and v2 to cylindrical and carthesian coordinates
      !::based on which plane is currently selected in the model
      !::
      subroutine car2cyl( v1, v2, r, phi, x, y, z )
      real*8,intent(in) :: v1,v2
      real*8,intent(inout) :: r, phi, x, y, z
      !::Make coordinate transformation
      if ( model%cylPlane .eq. phi_z_plane ) then
      !::Integrate over the phi-z plane, i.e. r = const
      !::v1 is then assumed to be phi and v2 to be z
          r = model%r0
          phi = v1      
          z = v2
      else
      !::Integrate over the r-phi plane, i.e. z = const
      !::v1 is assumed to be r and v2 to be phi
          r = v1
          phi = v2      
          z = model%z0
      endif
      x = r * cos( phi )
      y = r * sin( phi )
      
      end subroutine car2cyl
      
      subroutine getCylVec( phi, Hsol, B )
      real*8,intent(in) :: phi
      real*8,dimension(3),intent(in) :: Hsol
      real*8,dimension(3),intent(inout) :: B
       !::Convert to cylindrical coordinates
       B(1) = mu0 * (  cos(phi) * Hsol(1) + sin(phi) * Hsol(2) )
       B(2) = mu0 * ( -sin(phi) * Hsol(1) + cos(phi) * Hsol(2) )
       B(3) = mu0 * Hsol(3)
      
      end subroutine getCylVec
      
      
      
      
      !::
      !::Performs the box integral around the region specified 
      !::and over the Maxwell stress tensor using integral2 via quadpack
      !::xPos(2) defines min and max of the x range
      !::yPos(2) defines min and max of the y range
      !::zPos(2) defines min and max of the z range
      !::eps_abs defines the minimum absolute error
      !::eps_rel defines the minimum relative error      
      !::F is the resulting force
      !::ier is an integer describing the error state 
      !::of the integration (0 = no error)
      !:: abserr_tot is the total absolute error
      !::neval is the number of function evaluations used by the integrator
      subroutine boxIntegral(xPos, yPos, zPos, eps_abs, eps_rel, 
     +                       F, ier,abserr_tot, neval)
      
      real*8,dimension(2),intent(in) :: xPos,yPos,zPos
      real*8,dimension(3),intent(inout) :: F
      real*8,intent(inout) :: abserr_tot
      real*8,intent(in) :: eps_abs,eps_rel
      integer*8,intent(inout) :: ier,neval
      real*8 :: res,abserr
      integer*8 :: i,j
      integer*8,dimension(3,3,3) :: indices
      real*8,dimension(3,2,3) :: xp,yp,zp
     
      !::Reset input array
      F(:) = 0.
      
      !set to a large number
      abserr_tot = 1e9
      neval = 0
      
      !::set the indices for the x-component of the force
      indices(:,1,1) = (/ 2, 3, 1 /) !yz-plane
      indices(:,2,1) = (/ 1, 3, 2 /) !xz-plane
      indices(:,3,1) = (/ 1, 2, 3 /) !xy-plane
      
      xp(1,:,1) = yPos !yz-plane
      xp(2,:,1) = xPos !xz-plane
      xp(3,:,1) = xPos !xy-plane
      
      yp(1,:,1) = zPos !yz-plane
      yp(2,:,1) = zPos !xz-plane
      yp(3,:,1) = yPos !xy-plane
      
      zp(1,:,1) = xPos !yz-plane
      zp(2,:,1) = yPos !xz-plane
      zp(3,:,1) = zPos !xy-plane
      
      
      !::set the indices for the y-component of the force
      indices(:,1,2) = (/ 1, 3, 2 /) !xz-plane
      indices(:,2,2) = (/ 2, 3, 1 /) !yz-plane
      indices(:,3,2) = (/ 1, 2, 3 /) !xy-plane
      
      xp(1,:,2) = xPos !xz-plane
      xp(2,:,2) = yPos !yz-plane
      xp(3,:,2) = xPos !xy-plane
      
      yp(1,:,2) = zPos !xz-plane
      yp(2,:,2) = zPos !xz-plane
      yp(3,:,2) = yPos !xy-plane
     
      zp(1,:,2) = yPos !xz-plane
      zp(2,:,2) = xPos !yz-plane
      zp(3,:,2) = zPos !xy-plane
      
      !::set the indices for the z-component of the force
      indices(:,1,3) = (/ 1, 2, 3 /) !xy-plane
      indices(:,2,3) = (/ 1, 3, 2 /) !xz-plane
      indices(:,3,3) = (/ 2, 3, 1 /) !yz-plane
      
      xp(1,:,3) = xPos !xy-plane      
      xp(2,:,3) = xPos !xz-plane      
      xp(3,:,3) = yPos !yz-plane
            
      yp(1,:,3) = yPos !xy-plane
      yp(2,:,3) = zPos !xz-plane
      yp(3,:,3) = zPos !yz-plane
      
      zp(1,:,3) = zPos !xy-plane
      zp(2,:,3) = yPos !xz-plane
      zp(3,:,3) = xPos !yz-plane
      
      !::find the rotation matrices once and for all
      call findRotationMatrices( model%rotAngles, model%n, 
     +                            model%rotMat, model%rotMatInv)
      !::j loops over the three components
      do j=1,3
        if ( model%retForces(j) .eq. 0 ) then
          !:: i loops over the three surface pairs
          do i=1,3
          
              if ( i .eq. 1 ) then
                  !::The diagonal components are always the first surface pair
                  f_tens_ptr => F_int_diag
              elseif ( i .gt. 1 ) then
                  !::And then comes the two surface pairs with off-diagonal components
                  f_tens_ptr => F_int_off_diag
              endif    
              !::Set the plane indices
              model%inds(1:3) = indices(:,i,j)
              !::set the index of the first integration function in the integrand
              model%intInd = j
              !::lower plane 
               model%z0 = zp(i,1,j)
          
              !::direction of the other xy plane
!              call integral2( f_tens_ptr, xp(i,1,j), xp(i,2,j), 
!     +                                yp(i,1,j), yp(i,2,j), 
!     +                                eps_abs, eps_rel,
!     +                                res, abserr, neval, ier )
!              if ( abserr .lt. abserr_tot ) then
!                  abserr_tot = abserr
!              endif      
!              if ( ier .gt. 0 ) then
!                  exit
!              endif
          
              !::Note the sign accounts for the normal vector pointing the opposite
!              F(j) = F(j) - res
              !::upper plane
              model%z0 = zp(i,2,j)
           
              call integral2( f_tens_ptr, xp(i,1,j), xp(i,2,j), 
     +                                yp(i,1,j), yp(i,2,j), 
     +                                eps_abs, eps_rel,
     +                                res, abserr, neval, ier )
              if ( abserr .lt. abserr_tot ) then
                  abserr_tot = abserr
              endif      
              if ( ier .gt. 0 ) then
                  exit
              endif
     
              F(j) = F(j) + res
          
          
          enddo
          endif
      enddo
      
      end subroutine boxIntegral
      
      
      !::Returns the integrand for the off-diagonal part of the Maxwell stress
      !::tensor dotted with the surface area
      function F_int_off_diag( x, y )
      real*8 :: F_int_off_diag
      real*8,intent(in) :: x,y
      
      real*8,dimension(3) :: solPts,B,Hsol
      Hsol(:) = 0
      !::setup the points where we wish to get the solution
      solPts(model%inds(1)) = x
      solPts(model%inds(2)) = y
      solPts(model%inds(3)) = model%z0;

      !::call the subroutine to find the solution at the specified points
      call getSolution( solPts, 1, model%M, model%dims, 
     +                  model%pos, model%n, model%spacedim,
     +                  Hsol, model%rotMat, model%rotMatInv )
      
                                 
       B = mu0 * Hsol
    
       F_int_off_diag = 
     + 1. /mu0 * B(model%intInd) * B(model%inds(3))
    
       return

      end function F_int_off_diag
      
        !::Returns the integrand of the Maxwell stress tensor along the diagonal over
        !::the given surface area. 
        !::
        !::integration variables while .zInd gives the index of the plane
        !::1 = x, 2 = y and 3 = z
        function F_int_diag( x, y )
        real*8 :: F_int_diag
        real*8,intent(in) :: x,y
        
        real*8,dimension(3) :: solPts,B,Hsol
        real*8 :: Bnorm
            
        Hsol(:) = 0
        
        solPts(model%inds(1)) = x
        solPts(model%inds(2)) = y
        solPts(model%inds(3)) = model%z0;

    
        call getSolution( solPts, 1, model%M, model%dims, 
     +                  model%pos, model%n, model%spacedim,
     +                  Hsol, model%rotMat, model%rotMatInv )
	
        B = mu0 * Hsol
        Bnorm = sqrt( B(1)**2 + B(2)**2 + B(3)**2 )
                                 
        F_int_diag = 1./mu0 * 
     +                   ( B(model%inds(3))**2 - 0.5 * Bnorm**2 )
    
        
        end function F_int_diag
      
      end module model_glob


      
      
!     Gateway routine to do a box integral using data points from the MagStatPrism library call
!     The prhs is a pointer array with the following content
!     prhs(1) rotAngles (n,3), rotation angles for each prism
!     prhs(2) M (n,3) magnetization for each prism
!     prhs(3) dims (n,3) dimensions of each prism
!     prhs(4) pos(n,3) positions of each prism
!     prhs(5) n, integer
!     prhs(6) spacedim, integer
!     prhs(7) xlim (2), x bounds
!     prhs(8) ylim (2), y bounds
!     prhs(9) zlim (2), z bounds
!     prhs(10) eps_abs, real*8
!     prhs(11) eps_rel, real*8
!     prhs(12) retForces (optional), integer [3,1] array with 0 indicates false and otherwise true

!     plhs(1) Fout(3) output force vector [Fx,Fy,Fz]
!     plhs(2) ier, integer returned
!     plhs(3) max_err, real
!     plhs(4) neval, integer, number of function evaluations

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use model_glob

      implicit none
      
!     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

!     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxCreateNumericArray
      integer mxIsNumeric
      integer mxIsDouble
      integer mxIsInt32
      integer mxIsStruct
      integer mxIsCell
      integer mxClassIDFromClassName
      integer writeToConsole
      mwPointer mxGetField
      mwPointer mxGetFieldByNumber
      
      mwPointer mxGetM, mxGetN
      mwSize mxGetNumberOfDimensions,dimIteOut(1)
      mwPointer mxGetDimensions
      
      integer*4 ComplexFlag,classid

!     Pointers to input/output mxArrays:
      mwPointer :: stPr      
      real*8,dimension(2) :: xlim,ylim,zlim     
      real*8,dimension(3) :: Fout
      integer*8 :: ier,neval
      real*8 :: maxerr_tot,eps_abs, eps_rel
      
      
	  maxerr_tot = 0
	  !-----------------------------------------------------------------------	  	       
	     open(11,file='version_integral2_mex.txt',
     +   status='unknown',access='sequential',form='formatted',
     +                       position='rewind',action='write')     
	  write(11,*) "integral2_mex compiled on:"
      write(11,*) __DATE__
      write(11,*) __TIME__
      close(11)
      !-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------------
!     Check for proper number of arguments. 
      if(nrhs .lt. 11) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nInput',
     +                           'at least 11 inputs are required.')
      elseif(nlhs .gt. 4) then
         call mexErrMsgIdAndTxt ('MATLAB:magStat_mex:nOutput',
     +                           'Too many output arguments.')
      endif

!Check the type of the inputs      
      if ( .NOT. mxIsDouble(prhs(1)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input one should be double')
      elseif ( .NOT. mxIsDouble(prhs(2)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input two should be double')      
      elseif ( .NOT. mxIsDouble(prhs(3)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input three should be double')
      elseif ( .NOT. mxIsDouble(prhs(4)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input four should be Double')
      elseif ( .NOT. mxIsInt32(prhs(5)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input five should be an integer')
      elseif ( .NOT. mxIsInt32(prhs(6)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input six should be an integer')
      elseif ( .NOT. mxIsDouble(prhs(7)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input seven should be a double')
      elseif ( .NOT. mxIsDouble(prhs(8)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input eight should be a double')      
      elseif ( .NOT. mxIsDouble(prhs(9)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input nine should be a double')      
      elseif ( .NOT. mxIsDouble(prhs(10)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input ten should be a double')      
      elseif ( .NOT. mxIsDouble(prhs(11)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input 11 should be a double')      
      elseif (nrhs .ge. 12 .and. .NOT. mxIsInt32(prhs(12)) ) then
          call mexErrMsgIdAndTxt ('MATLAB:Matlab_single_mex:DataType',
     +                           'Input 12 should be an integer')
      endif
      
      !::Allocate input data
      call mxCopyPtrToInteger4(mxGetPr(prhs(5)), model%n, 1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(6)), model%spacedim, 1)
      
      allocate( model%M(model%n,3), model%dims(model%n,3), 
     +                model%pos(model%n,3), model%rotAngles(model%n,3) )
      !::Copy input data
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), model%rotAngles, 
     +                      model%n * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), model%M, model%n * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), model%dims, model%n * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(4)), model%pos, model%n * 3 )
      
      call mxCopyPtrToReal8(mxGetPr(prhs(7)), xlim, 2 )
      call mxCopyPtrToReal8(mxGetPr(prhs(8)), ylim, 2 )
      call mxCopyPtrToReal8(mxGetPr(prhs(9)), zlim, 2 )
      
      call mxCopyPtrToReal8(mxGetPr(prhs(10)), eps_abs, 2 )
      call mxCopyPtrToReal8(mxGetPr(prhs(11)), eps_rel, 2 )
      
      if ( nrhs .ge. 12 ) then         
          call mxCopyPtrToInteger4(mxGetPr(prhs(12)), model%retForces,3) 
      else
          model%retForces(:) = 1
      endif
     
      !::Load the result back to Matlab
      call boxIntegral( xlim, ylim, zlim, eps_abs, eps_rel, Fout, 
     +                  ier, maxerr_tot, neval )
     
      !call cylIntegral( xlim(1), zlim, eps_abs, eps_rel, Fout, 
      !+                  ier, maxerr_tot, neval )
      ComplexFlag = 0;
      
      plhs(1) = mxCreateDoubleMatrix( 3, 1, ComplexFlag )
           
      call mxCopyReal8ToPtr( Fout, mxGetPr( plhs(1) ), 3 )
      
      if (nlhs .gt. 1) then
          classid = mxClassIDFromClassName('int64')
          plhs(2) = mxCreateNumericArray( 1, 1, classid, ComplexFlag )
      call mxCopyInteger4ToPtr( int( ier, 4 ), mxGetPr(plhs(2)), 1 )
      endif
      if ( nlhs .gt. 2) then          
          plhs(3) = mxCreateDoubleMatrix( 1, 1, ComplexFlag )
          call mxCopyInteger4ToPtr( maxerr_tot, mxGetPr(plhs(3)), 1 )
      endif
      if (nlhs .gt. 3) then
          classid = mxClassIDFromClassName('int64')
          plhs(4) = mxCreateNumericArray( 1, 1, classid, ComplexFlag )
      call mxCopyInteger4ToPtr( int( neval, 4 ), mxGetPr(plhs(4)), 1 )
      endif
      
      
      deallocate( model%rotAngles, model%M, model%dims, model%pos )
      deallocate( model%rotMat, model%rotMatInv )
      
      
      end subroutine
      
      
      
            
      
      
      