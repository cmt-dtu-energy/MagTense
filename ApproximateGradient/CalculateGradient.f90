
    module GetGradient
    
    use TileNComponents
	use FaceComponents

    implicit none
      
    contains
    
    !--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !>
    !!Calculates gradient given a number of (different) tiles
    !!Input arguments:
    !!@param tiles: Array of type MagTile, size n_tiles
    !!@param gradient the average gradient over each cell, size [n_ele,3]
    !!@param n_tiles, the number of tiles
	!!@param faces contains the three important aspects of each boundary between cells:
	!!    %area: The area of the cell face
	!!    %normal: the orientation of the cell face
	!!    %offset: the location of the cell face center
	!!@param n_faces, the number of faces
    !!@param B is a matrix with n_faces rows and n_tiles columns with B_if = +1 if face f faces outwards from cell i, -1 if inwards and 0 if face f does not belong to cell i
	!!@todo: make B sparse.
	!!
    subroutine getGradientFromFaceAverage( tiles, gradient, faces, B, n_tiles, n_faces )
		integer,intent(in) :: n_tiles,n_faces
	    type(MagTile),intent(in),dimension(n_tiles) :: tiles
		type(TileFace),intent(in),dimension(n_faces) :: faces
		integer,dimension(n_tiles,n_faces),intent(in) :: B
		real,dimension(n_tiles,3),intent(inout) :: gradient		
		integer :: dummy
		real :: tileVolume, face_avg
		real,dimension(3) :: grad_tmp
		
		do i=1,n_tiles
		    grad_tmp = 0
		    do j=1,n_faces
			    if (abs(B(i,j)) .eq. 1) then
				    !! Calculates the face value average
				    dummy = 0
					face_avg = 0
				    do k=1,n_tiles
					    if (abs(B(k,j)) .eq. 1) then						
						    face_avg = face_avg + tiles(k)
							dummy = dummy + 1
				        endif						
					enddo
					face_avg = face_avg / dummy
					!! Add face contribution to gradient
					grad_tmp = grad_tmp + B(i,j) * faces(j)%normal * faces(j)%area * face_avg
				endif				
			enddo
		    tileVolume = tiles(i)%a * tiles(i)%b * tiles(i)%c			
			gradient(i) = grad_tmp / tileVolume
		enddo
	end subroutine getGradientFromFaceAverage
	
	subroutine getInverse44(A, AINV)
		!! Performs a direct calculation of the inverse of a 4×4 matrix.
		!! Shamelessly stolen from http://fortranwiki.org/fortran/show/Matrix+inversion
		DOUBLE PRECISION, DIMENSION(4,4), INTENT(IN)  :: A
		DOUBLE PRECISION, DIMENSION(4,4), INTENT(OUT) :: AINV
		DOUBLE PRECISION :: detinv

		! Calculate the inverse determinant of the matrix
		detinv = &
		1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
		- A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
		+ A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
		- A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))
	
		! Calculate the inverse of the matrix
		AINV(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
		AINV(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
		AINV(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
		AINV(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
		AINV(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
		AINV(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
		AINV(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
		AINV(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
		AINV(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
		AINV(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
		AINV(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
		AINV(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
		AINV(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
		AINV(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
		AINV(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
		AINV(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))	
	end subroutine
	
	subroutine PrepareLeastSquaresInterpolation( tiles, LSQRMatrix, faces, n_tiles, n_faces)
		!! Sets up the matrix to be applied for the least squares method
		integer,intent(in) :: n_tiles,n_faces
	    type(MagTile),intent(in),dimension(n_tiles) :: tiles
		type(TileFace),intent(in),dimension(n_faces) :: faces
		real,dimension(n_faces, n_tiles),intent(inout) :: LSQRMatrix
		integer,dimension(n_tiles) :: diffx, diffy, diffz, LSQRVector, weights		
		real,dimension(4,4) :: G, Ginv
		real,dimension(n_tiles,4) :: H !! Note that H is transposed from that seen in documentation
		integer :: dummy
		real :: tileVolume, face_avg
		real,dimension(3) :: grad_tmp
		
		do i=1,n_faces			
			diffx = tiles(:)%offset(1) - faces(i)%offset(1)
			diffy = tiles(:)%offset(2) - faces(i)%offset(2)
			diffz = tiles(:)%offset(3) - faces(i)%offset(3)
			weights = 1.0 / (abs(diffx) + abs(diffy) + abs(diffz))
			G = reshape( &
			    sum(weights),      sum(weights*diffx),      sum(weights*diffy),      sum(weights*diffz), &
				sum(weights*diffx),sum(weights*diffx*diffx),sum(weights*diffx*diffy),sum(weights*diffx*diffz), &
				sum(weights*diffy),sum(weights*diffx*diffy),sum(weights*diffy*diffy),sum(weights*diffy*diffz), &
				sum(weights*diffz),sum(weights*diffx*diffz),sum(weights*diffy*diffz),sum(weights*diffz*diffz), &
			    (/4,4/))
			H(:,1) = weights
			H(:,2) = weights*diffx
			H(:,3) = weights*diffy
			H(:,4) = weights*diffz
			call getInverse44(G,Ginv)
			LSQRMatrix(i,:) = matmul(H,Ginv(1,:))
		enddo		
	end subroutine PrepareLeastSquaresInterpolation
	
	subroutine getGradientFromLeastSquaresInterpolation( tile, gradient, faces, B, LSQRMatrix, n_tiles, n_faces )
		integer,intent(in) :: n_tiles,n_faces
		type(MagTile),intent(in),dimension(n_tiles) :: tiles
		real,intent(in),dimension(n_faces,n_tiles) :: LSQRMatrix
		type(TileFace),intent(in),dimension(n_faces) :: faces
		integer,dimension(n_tiles,n_faces),intent(in) :: B
		real,dimension(n_tiles,3),intent(inout) :: gradient
		real :: tileVolume
		real,dimension(n_faces) :: face_est
		real,dimension(3) :: grad_tmp
		
		!! Calculates the face value estimate
		face_est = matmul(LSQRMatrix,tiles)
	    
		do i=1,n_tiles
		    grad_tmp = 0
		    do j=1,n_faces
			    if (abs(B(i,j)) .eq. 1) then
				    !! Add contribution to gradient
					grad_tmp = grad_tmp + B(i,j) * faces(j)%normal * faces(j)%area * face_est(j)
				endif				
			enddo
		    tileVolume = tiles(i)%a * tiles(i)%b * tiles(i)%c			
			gradient(i) = grad_tmp / tileVolume
		enddo
	end subroutine getGradientFromLeastSquaresInterpolation

end module GetGradient