#include "fintrf.h"  
      module model_glob
      USE MagPrism_Call
      USE Parameters_CALL
      USE NumInt
      implicit none
      
      !::Data model for the integration of the Maxwell stress tensor
      !type, extends( dataCollectionBase ) :: dataCollectionExt
      !    class( magStatModel ), pointer :: model
      !end type dataCollectionExt
      
      type, extends( dataCollectionModelBase) :: magStatModel
          integer :: n,spacedim
          real*8,dimension(:,:),allocatable :: dims,pos,M,rotAngles
          real*8,dimension(:,:,:),allocatable :: rotMat,rotMatInv         
      end type magStatModel
      
      
      real*8,parameter :: mu0=4*pi*1e-7
            
       
      contains
      
      
      
      subroutine handleError( dat, abserr_tot )      
      class( dataCollectionBase ), intent(in), target :: dat
      real*8,intent(inout),dimension(2) :: abserr_tot
      
      if ( dat%abserr_x .lt. abserr_tot(1) ) then
        abserr_tot(1) = dat%abserr_x
      endif      
      if ( dat%abserr_y .lt. abserr_tot(2) ) then
        abserr_tot(2) = dat%abserr_y
      endif      
      if ( dat%ier_x .gt. 0 .or. dat%ier_y .gt. 0 ) then
        !stop
      endif
      end subroutine handleError
    
      !::The nine tensor components of the Maxwell stress tensor  
      function F_int_12( dat )
      real*8 :: F_int_12
      class( dataCollectionBase ), intent(in), target :: dat
      
      real*8,dimension(3) :: B
          
      !::Get the field from the model solution
      call getB_field( dat, B )
                                 
      !::Calculate the off diagonal component. 
       F_int_12 = 1. /mu0 * B(1) * B(2)
    
       return

      end function F_int_12
      
      function F_int_13( dat )
      real*8 :: F_int_13
      class( dataCollectionBase ), intent(in), target :: dat
      
      real*8,dimension(3) :: B
          
      !::Get the field from the model solution
      call getB_field( dat, B )
                                 
      !::Calculate the off diagonal component
       F_int_13 = 1. /mu0 * B(1) * B(3)
    
       return
      end function F_int_13
      
      function F_int_23( dat )
      real*8 :: F_int_23
      class( dataCollectionBase ), intent(in), target :: dat
      
      real*8,dimension(3) :: B
          
      !::Get the field from the model solution
      call getB_field( dat, B )
                                 
      !::Calculate the off diagonal component
       F_int_23 = 1. /mu0 * B(2) * B(3)
    
       return

      end function F_int_23
      
      function F_int_21( dat )
      real*8 :: F_int_21
      class( dataCollectionBase ), intent(in), target :: dat
      
      F_int_21 = F_int_12( dat )
      
      end function F_int_21
      
      function F_int_31( dat )
      real*8 :: F_int_31
      class( dataCollectionBase ), intent(in), target :: dat
      
      F_int_31 = F_int_13( dat )
      
      end function F_int_31
      
      function F_int_32( dat )
      real*8 :: F_int_32
      class( dataCollectionBase ), intent(in), target :: dat
      
      F_int_32 = F_int_23( dat )
      
      end function F_int_32
                  
      function F_int_11( dat )
      real*8 :: F_int_11
      class( dataCollectionBase ), intent(in), target :: dat
      real*8,dimension(3) :: B
      real*8 :: Bnorm
      
      call getB_field( dat, B )
      
      Bnorm = sqrt( B(1)**2 + B(2)**2 + B(3)**2 )
                                 
      F_int_11 = 1./mu0 * 
     +                   ( B(1)**2 - 0.5 * Bnorm**2 )
      
      end function F_int_11
      
      function F_int_22( dat )
      real*8 :: F_int_22
      class( dataCollectionBase ), intent(in), target :: dat
      real*8,dimension(3) :: B
      real*8 :: Bnorm
      
      call getB_field( dat, B )
      
      Bnorm = sqrt( B(1)**2 + B(2)**2 + B(3)**2 )
                                 
      F_int_22 = 1./mu0 * 
     +                   ( B(2)**2 - 0.5 * Bnorm**2 )
      
      end function F_int_22
      
      function F_int_33( dat )
      real*8 :: F_int_33
      class( dataCollectionBase ), intent(in), target :: dat
      real*8,dimension(3) :: B
      real*8 :: Bnorm
      
      call getB_field( dat, B )
      
      Bnorm = sqrt( B(1)**2 + B(2)**2 + B(3)**2 )
                                 
      F_int_33 = 1./mu0 * 
     +                   ( B(3)**2 - 0.5 * Bnorm**2 )
      
      end function F_int_33
      
        !::Get the field from the model solution
        subroutine getB_field( dat, B )
        class( dataCollectionBase ), intent(in), target :: dat
        class( magStatModel ), pointer :: model
        real*8,dimension(3),intent(inout) :: B
        
        real*8,dimension(3) :: Hsol, solPts
        integer*8,dimension(3) :: inds
        
        !::Make the down cast to get the actual data model set
        call getDataModel( dat%model, model )
        
        !::Make sure the coordinate indices are mapped according to the normal vector
        if ( abs(dat%n_vec(1)) .eq. 1 ) then
            !yz-plane
            inds(1) = 2
            inds(2) = 3
            inds(3) = 1
        else if ( abs(dat%n_vec(2)) .eq. 1 ) then
            !xz-plane
            inds(1) = 1
            inds(2) = 3
            inds(3) = 2
        else if ( abs(dat%n_vec(3)) .eq. 1 ) then
            !xy-plane
            inds(1) = 1
            inds(2) = 2
            inds(3) = 3
        endif
        
        
        !::Reset the solution array
        Hsol(:) = 0
        
        solPts(inds(1)) = dat%x
        solPts(inds(2)) = dat%y
        solPts(inds(3)) = dat%z0;
    
        call getSolution( solPts, 1, model%M, model%dims, model%pos, 
     +                    model%n, model%spacedim, Hsol, 
     +                    model%rotMat, model%rotMatInv )
	
        B = mu0 * Hsol
        end subroutine getB_field
        
        !::Make the down cast
        subroutine getDataModel( modelIn, modelOut )
        class( dataCollectionModelBase ), intent(in), target :: modelIn
        class( magStatModel ), intent(out), pointer :: modelOut
        
            select type( modelIn )
            class is (magStatModel)
                modelOut => modelIn
                class default
                    write(*,*) 'Problem with the type cast!?'
            end select
        end subroutine getDataModel
        
      
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
!     plhs(2) ier(2), integer returned
!     plhs(3) max_err(2), real
!     plhs(4) neval(2), integer, number of function evaluations

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
      real*8,dimension(3) :: Fout
      integer*8,dimension(2) :: ier,neval
      real*8 :: eps_abs, eps_rel
      real*8,dimension(2) :: maxerr_tot
      class( magStatModel ), pointer :: datModel
      type( dat_ptr ), dimension(18) :: dat_arr
      type( surf_carth ) :: surf
      integer*4,dimension(3) :: rV
      integer*8,dimension(3) :: retV
      integer :: i,j,ind
      
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
            
      allocate( datModel )
      
      !::Allocate input data

      call mxCopyPtrToInteger4(mxGetPr(prhs(5)), datModel%n, 1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(6)), datModel%spacedim, 1)
      
      allocate(datModel%M(datModel%n,3),datModel%dims(datModel%n,3))
      allocate(datModel%pos(datModel%n,3) )
      allocate(datModel%rotAngles(datModel%n,3) )
      !::Copy input data
      call mxCopyPtrToReal8(mxGetPr(prhs(1)), datModel%rotAngles, 
     +                      datModel%n * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), datModel%M, 
     +                      datModel%n * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(3)), datModel%dims, 
     +                      datModel%n * 3 )
      call mxCopyPtrToReal8(mxGetPr(prhs(4)), datModel%pos, 
     +                      datModel%n * 3 )
      
      call mxCopyPtrToReal8(mxGetPr(prhs(7)), surf%x, 2 )
      call mxCopyPtrToReal8(mxGetPr(prhs(8)), surf%y, 2 )
      call mxCopyPtrToReal8(mxGetPr(prhs(9)), surf%z, 2 )
      
      call mxCopyPtrToReal8(mxGetPr(prhs(10)), eps_abs, 2 )
      call mxCopyPtrToReal8(mxGetPr(prhs(11)), eps_rel, 2 )
      
      if ( nrhs .ge. 12 ) then         
          call mxCopyPtrToInteger4(mxGetPr(prhs(12)), rV,3)           
          retV(:) = rV(:)
      else
          retV(:) = 0
      endif
      !::Initialize each member of the dat_arr
      do i=1,3
          do j=1,6             
              ind = (i-1) * 6 + j
              
              allocate( dat_arr(ind)%dat )
              
              dat_arr(ind)%dat%model => datModel
              
              dat_arr(ind)%dat%epsabs = eps_abs
              dat_arr(ind)%dat%epsrel = eps_rel
              
          enddo
      enddo
      
      dat_arr(1)%dat%f_ptr => F_int_11
      dat_arr(2)%dat%f_ptr => F_int_11
      dat_arr(3)%dat%f_ptr => F_int_12
      dat_arr(4)%dat%f_ptr => F_int_12
      dat_arr(5)%dat%f_ptr => F_int_13
      dat_arr(6)%dat%f_ptr => F_int_13
      
      dat_arr(7)%dat%f_ptr => F_int_21
      dat_arr(8)%dat%f_ptr => F_int_21
      dat_arr(9)%dat%f_ptr => F_int_22
      dat_arr(10)%dat%f_ptr => F_int_22
      dat_arr(11)%dat%f_ptr => F_int_23
      dat_arr(12)%dat%f_ptr => F_int_23
      
      dat_arr(13)%dat%f_ptr => F_int_31
      dat_arr(14)%dat%f_ptr => F_int_31
      dat_arr(15)%dat%f_ptr => F_int_32
      dat_arr(16)%dat%f_ptr => F_int_32
      dat_arr(17)%dat%f_ptr => F_int_33
      dat_arr(18)%dat%f_ptr => F_int_33

      
      
      !::find the rotation matrices once and for all
      call findRotationMatrices( datModel%rotAngles, datModel%n, 
     +                            datModel%rotMat, datModel%rotMatInv)
     
      call surface_integral_carth( surf, dat_arr, retV, 
     +                             handleError, Fout, ier, neval )
     
      !::Load the result back to Matlab
     
      
      ComplexFlag = 0
      
      plhs(1) = mxCreateDoubleMatrix( 3, 1, ComplexFlag )
           
      call mxCopyReal8ToPtr( Fout, mxGetPr( plhs(1) ), 3 )
      
      if (nlhs .gt. 1) then
          classid = mxClassIDFromClassName('int64')
          plhs(2) = mxCreateNumericArray( 2, 1, classid, ComplexFlag )
      call mxCopyInteger4ToPtr( int( ier, 4 ), mxGetPr(plhs(2)), 2 )
      endif
      if ( nlhs .gt. 2) then          
          plhs(3) = mxCreateDoubleMatrix( 2, 1, ComplexFlag )
          call mxCopyInteger4ToPtr( maxerr_tot, mxGetPr(plhs(3)), 2 )
      endif
      if (nlhs .gt. 3) then
          classid = mxClassIDFromClassName('int64')
          plhs(4) = mxCreateNumericArray( 2, 1, classid, ComplexFlag )
      call mxCopyInteger4ToPtr( int( neval, 4 ), mxGetPr(plhs(4)), 2 )
      endif
      
      
            
      end subroutine
      
      
      
            
      
      
      