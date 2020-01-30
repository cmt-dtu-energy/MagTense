

#include <stdlib.h>
#include <stdio.h>
#include <cusparse.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>


 
 /**
 Demag matrices to be stored in GPU memory
 */
 float* d_Kxx = NULL;
 float* d_Kxy = NULL;
 float* d_Kxz = NULL;
 float* d_Kyy = NULL;
 float* d_Kyz = NULL;
 float* d_Kzz = NULL;
 
 
 //Magnetization vectors on GPU
 float* d_Mx = NULL;
 float* d_My = NULL;
 float* d_Mz = NULL;
 
 //Magnetic field vectors on GPU
 float* d_Hx = NULL;
 float* d_Hy = NULL;
 float* d_Hz = NULL;
 
 //size of the matrices (n_K x n_K)
 int n_K;
 
 /**
 Sparse demag matrices to be stored in the GPU memory
 */
 
 struct CUSparse
 {

	//data values
	float* values = NULL;
	//column indices (size nnz)
	int* cols;
	//row indices (size n + 1  with the last element equal to nnz)
	int* rows;
	//no of matrix elements (matrices are n x n )
	int n;
	//no. of non zero element
	int nnz;
 };
	
 
 CUSparse spKxx;
 CUSparse spKxy;
 CUSparse spKxz;
 CUSparse spKyy;
 CUSparse spKyz;
 CUSparse spKzz;
 
 //general handle for cuBlas (initalized once)
 cublasHandle_t handle = NULL;
 
 
 //handle to the sparse matrix description
cusparseMatDescr_t sparse_descr = NULL;

//handle to the sparse matrix blas
cusparseHandle_t sparse_handle = NULL; 
 
 void loadSparseToDevice( float* values, int* colInds, int* rowInds, CUSparse* mat, int nnz_ );
 void freeSparseMatrix( CUSparse* mat );
 /**
 Kaspar K. Nielsen, kasparkn@gmail.com.dk, 2019
 Initializes the sparse matrices for later use in the matrix-vector multiplications
 n is the no. of rows and cols in the matrices
 nnz is the no. of non-zero elements
 mat_no 1...6 identifies which matrix to load (Kxx = 1, Kxy = 2, Kxz = 3, Kyy = 4, Kyz = 5, Kzz = 6 )
 values float array with the non-zero values of the matrix (length nnz)
 colInds int array with the column indices of the CSR stored matrix (length nnz)
 rowInds the indices of the first non-zero element in each row (size nrows+1 and the last element has to be equal to nnz + rowInds[0] )
 */
 void cu_initDemagMatrices_sparse( const int* n, const int* nnz, const int* mat_no, float* values, int* colInds, int* rowInds )
 {
	
	
	 //should only be done once
	 if ( sparse_handle == NULL )
	 {
		 //handle to the sparse matrix multiplier in CUDA
		cusparseStatus_t status = cusparseCreate(&sparse_handle);
		
		//no. of rows and columns
		spKxx.n = *n;
		spKxy.n = *n;
		spKxz.n = *n;
		spKyy.n = *n;
		spKyz.n = *n;
		spKzz.n = *n;
		cudaError_t err;
		//allocate memory for the magnetization vectors
		err = cudaMalloc((void**) &d_Mx, *n * sizeof(float));
		err = cudaMalloc((void**) &d_My, *n * sizeof(float));
		err = cudaMalloc((void**) &d_Mz, *n * sizeof(float));
		
		//allocate memory for the magnetic field vectors
		err = cudaMalloc((void**) &d_Hx, *n * sizeof(float));
		err = cudaMalloc((void**) &d_Hy, *n * sizeof(float));
		err = cudaMalloc((void**) &d_Hz, *n * sizeof(float));
			 
		status = cusparseCreateMatDescr(&sparse_descr);
		cusparseSetMatType(sparse_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(sparse_descr,CUSPARSE_INDEX_BASE_ONE);
		
		/*Kaspar: This code will not work on Windows before the next major release thanks to nvidia
		//The handle on the device for the magnetization vectors
		cusparseCreateDnVec( &d_Mx_descr, n_sparse, d_Mx, CUDA_R_32F);
		cusparseCreateDnVec( &d_My_descr, n_sparse, d_My, CUDA_R_32F);
		cusparseCreateDnVec( &d_Mz_descr, n_sparse, d_Mz, CUDA_R_32F);
		
		//handles for the field on the device
		cusparseCreateDnVec( &d_Hx_descr, n_sparse, d_Hx, CUDA_R_32F);
		cusparseCreateDnVec( &d_Hy_descr, n_sparse, d_Hy, CUDA_R_32F);
		cusparseCreateDnVec( &d_Hz_descr, n_sparse, d_Hz, CUDA_R_32F);
			*/	
		
		
	 }
	
	switch( *mat_no )
	{
		case 1:
			
			loadSparseToDevice( values, colInds, rowInds, &spKxx, *nnz );
								
			break;
		case 2:
			loadSparseToDevice( values, colInds, rowInds, &spKxy, *nnz );
			
			break;
		case 3:
			loadSparseToDevice( values, colInds, rowInds, &spKxz, *nnz );
			
			break;
		case 4:
			loadSparseToDevice( values, colInds, rowInds, &spKyy, *nnz );
			
			break;
		case 5:
			loadSparseToDevice( values, colInds, rowInds, &spKyz, *nnz );
			
			break;
		case 6:
			loadSparseToDevice( values, colInds, rowInds, &spKzz, *nnz );
			
			break;
			
	}
	
	
	
 }
 
 
 void loadSparseToDevice( float* values, int* colInds, int* rowInds, CUSparse* mat, int nnz_ )
 {
	cudaError_t err;
	//set the no. of non-zero entries
 	mat->nnz = nnz_;
	 //allocate the row inds (+1 in size as the last element contains nnz + rowInds(0)
	err = cudaMalloc((void**) &(mat->rows), (mat->n + 1) * sizeof(int));
	//allocate the column inds
	err = cudaMalloc((void**) &(mat->cols), mat->nnz * sizeof(int));
	//allocate the values array
	err = cudaMalloc((void**) &(mat->values), mat->nnz * sizeof(float));
	
	//copy to device
	err = cudaMemcpy( mat->rows, rowInds, (mat->n+1) * sizeof(int),cudaMemcpyHostToDevice);
	err = cudaMemcpy( mat->cols, colInds, mat->nnz * sizeof(int),cudaMemcpyHostToDevice);
	err = cudaMemcpy( mat->values, values, mat->nnz * sizeof(float),cudaMemcpyHostToDevice);
	
	
	
	//init the sparse matrix handles
	/*Kaspar: This code will not work before the next major release thanks to nvidia
	cusparseCreateCsr(&mat, n_sparse, n_sparse, nnz_, d_rowInds, d_colInds, d_values,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ONE, CUDA_R_32F);
									  */
 }
 
 /**
 Kaspar K. Nielsen, kasparkn@gmail.com.dk, 2019
 Does the matrix-vector product of a sparse matrix and a dense vector
 and returns the demag field on the form:
 Hx = Kxx * Mx + Kxy * My + Kxz * Mz
 Hy = Kyz * Mx + Kyy * My + Kzz * Mz
 Hz = Kzx * Mx + Kzy * My + Kzz * Mz
 
 */
 void cu_MVMult_GetH_sparse(const float* Mx, const float* My, const float* Mz, float* Hx, float* Hy, float* Hz, const float* pref)
 {	 
	
	cusparseStatus_t status;
	
	
	
	int n = spKxx.n;
	//Copy Mx, My and Mz to device memory	
	cudaMemcpy( d_Mx, Mx, n * sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( d_My, My, n * sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( d_Mz, Mz, n * sizeof(float), cudaMemcpyHostToDevice );
	
	float alpha = *pref;
	float beta = 0.0;
	//by setting beta = 0 in the first call we ensure that the previous values of d_Hx are irrelevant
	//as this operation does the following:
	//d_Hx = alpha * d_Kxx * d_Mx + beta * d_Hx
	status = cusparseScsrmv( sparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, spKxx.nnz,
                           &alpha, sparse_descr, spKxx.values, spKxx.rows, spKxx.cols,
                           d_Mx, &beta, d_Hx);
	
	//set beta = 1 so as to keep the previous set values of d_Hx
	beta = 1.0;
	//Hx = Kxy * My
	status = cusparseScsrmv( sparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, spKxy.nnz,
                           &alpha, sparse_descr, spKxy.values, spKxy.rows, spKxy.cols,
                           d_My, &beta, d_Hx);
	//Hx = Kxz * Mz
	status = cusparseScsrmv( sparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, spKxz.nnz,
                           &alpha, sparse_descr, spKxz.values, spKxz.rows, spKxz.cols,
                           d_Mz, &beta, d_Hx);						   
						   
						   
	//reset beta = 0 such that any previous values of d_Hy are erased
	beta = 0.;
	status = cusparseScsrmv( sparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, spKxy.nnz,
                           &alpha, sparse_descr, spKxy.values, spKxy.rows, spKxy.cols,
                           d_Mx, &beta, d_Hy);
	
	//set beta = 1 so as to keep the previous set values of d_Hy
	beta = 1.0;
	//Hy = Kyy * My
	status = cusparseScsrmv( sparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, spKyy.nnz,
                           &alpha, sparse_descr, spKyy.values, spKyy.rows, spKyy.cols,
                           d_My, &beta, d_Hy);
	//Hy = Kyz * Mz
	status = cusparseScsrmv( sparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, spKyz.nnz,
                           &alpha, sparse_descr, spKyz.values, spKyz.rows, spKyz.cols,
                           d_Mz, &beta, d_Hy);						   
						   
						   
	//reset beta = 0 such that any previous values of d_Hz are erased
	beta = 0.;
	status = cusparseScsrmv( sparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, spKxz.nnz,
                           &alpha, sparse_descr, spKxz.values, spKxz.rows, spKxz.cols,
                           d_Mx, &beta, d_Hz);
	
	//set beta = 1 so as to keep the previous set values of d_Hy
	beta = 1.0;
	//Hz = Kzy * My
	status = cusparseScsrmv( sparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, spKyz.nnz,
                           &alpha, sparse_descr, spKyz.values, spKyz.rows, spKyz.cols,
                           d_My, &beta, d_Hz);
	//Hz = Kzz * Mz
	status = cusparseScsrmv( sparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, spKzz.nnz,
                           &alpha, sparse_descr, spKzz.values, spKzz.rows, spKzz.cols,
                           d_Mz, &beta, d_Hz);									   
						   
	//kaspar: All the calls below will not work before the next major release by nvidia
	
	//Hx = Kxx * Mx + Kxy * My + Kxz * Mz	
	//Kxx * Mx
	/*
	cusparseSpMV( sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, d_Kxx_descr, d_Mx_descr, &beta, d_Hx_descr, CUDA_R_32F,
                                 CUSPARSE_MV_ALG_DEFAULT, dBuffer);
			
	beta = 1.0;
	//Kxy * My
	cusparseSpMV( sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, d_Kxy_descr, d_My_descr, &beta, d_Hx_descr, CUDA_R_32F,
                                 CUSPARSE_MV_ALG_DEFAULT, dBuffer);
								 
	cusparseSpMV( sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, d_Kxz_descr, d_Mz_descr, &beta, d_Hx_descr, CUDA_R_32F,
                                 CUSPARSE_MV_ALG_DEFAULT, dBuffer);

			   
	//Hy = Kxy * Mx + Kyy * My + Kyz * Mz	
	beta = 0.;
	cusparseSpMV( sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, d_Kxy_descr, d_Mx_descr, &beta, d_Hy_descr, CUDA_R_32F,
                                 CUSPARSE_MV_ALG_DEFAULT, dBuffer);
			
	beta = 1.0;	
	cusparseSpMV( sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, d_Kyy_descr, d_My_descr, &beta, d_Hy_descr, CUDA_R_32F,
                                 CUSPARSE_MV_ALG_DEFAULT, dBuffer);
								 
	cusparseSpMV( sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, d_Kyz_descr, d_Mz_descr, &beta, d_Hy_descr, CUDA_R_32F,
                                 CUSPARSE_MV_ALG_DEFAULT, dBuffer);
								 
	//Hz = Kzx * Mx + Kzy * My + Kzz * Mz		
	beta = 0.;
	cusparseSpMV( sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, d_Kxz_descr, d_Mx_descr, &beta, d_Hz_descr, CUDA_R_32F,
                                 CUSPARSE_MV_ALG_DEFAULT, dBuffer);
			
	beta = 1.0;
	//Kxy * My
	cusparseSpMV( sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, d_Kyz_descr, d_My_descr, &beta, d_Hz_descr, CUDA_R_32F,
                                 CUSPARSE_MV_ALG_DEFAULT, dBuffer);
								 
	cusparseSpMV( sparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, d_Kzz_descr, d_Mz_descr, &beta, d_Hz_descr, CUDA_R_32F,
                                 CUSPARSE_MV_ALG_DEFAULT, dBuffer);
			   
			   */
	//copy the solution back
	cudaMemcpy( Hx, d_Hx, n * sizeof(float), cudaMemcpyDeviceToHost );
	cudaMemcpy( Hy, d_Hy, n * sizeof(float), cudaMemcpyDeviceToHost );
	cudaMemcpy( Hz, d_Hz, n * sizeof(float), cudaMemcpyDeviceToHost );
 }
 
 
 /**
 Kaspar K. Nielsen, kasparkn@gmail.com.dk, DTU 2019
 Method for transfering the demagnetization matrices from Fortran
 @param Kxx, Kxy, Kxz, Kyy, Kyz and Kzz are 2D arrays of size n x n implemented as an array of n pointers each pointing to an array of n values
 @param n is the size of the square matrix
 */
 void cu_initDemagMatrices( const float* Kxx, const float* Kxy, const float* Kxz, const float* Kyy, const float* Kyz, const float* Kzz, int* n)
 {	 
	/* FILE *fp;

	fp = fopen("test.txt", "w+");
	fprintf(fp,"%f %d\n", Kxx[0], *n );
	//fprintf(fp, "This is testing for fprintf...\n");
	//fputs("This is testing for fputs...\n", fp);
	fclose(fp);
	 */
	 
	 cublasStatus_t stat;
	 n_K = *n;
	 
	 //Allocate the device (GPU) arrays
	 size_t bytes = n_K * n_K * sizeof(float);
	 
	 cudaMalloc( &d_Kxx, bytes );
	 cudaMalloc( &d_Kxy, bytes );
	 cudaMalloc( &d_Kxz, bytes );
	 cudaMalloc( &d_Kyy, bytes );
	 cudaMalloc( &d_Kyz, bytes );
	 cudaMalloc( &d_Kzz, bytes );
	 
	 
	 
	 //copy the demag tensors to the device	 
	 stat = cublasSetMatrix (n_K,n_K, sizeof (float), Kxx, n_K, d_Kxx, n_K );
	 stat = cublasSetMatrix (n_K,n_K, sizeof (float), Kxy, n_K, d_Kxy, n_K );
	 stat = cublasSetMatrix (n_K,n_K, sizeof (float), Kxz, n_K, d_Kxz, n_K );
	 stat = cublasSetMatrix (n_K,n_K, sizeof (float), Kyy, n_K, d_Kyy, n_K );
	 stat = cublasSetMatrix (n_K,n_K, sizeof (float), Kyz, n_K, d_Kyz, n_K );
	 stat = cublasSetMatrix (n_K,n_K, sizeof (float), Kzz, n_K, d_Kzz, n_K );
	 
	//allocate the internal M and H vectors
	bytes = n_K * sizeof(float);
	cudaMalloc( &d_Mx, bytes );
	cudaMalloc( &d_My, bytes );
	cudaMalloc( &d_Mz, bytes );

	cudaMalloc( &d_Hx, bytes );
	cudaMalloc( &d_Hy, bytes );
	cudaMalloc( &d_Hz, bytes );
	
	//initialize the cuBlas handle
	cublasCreate(&handle);
	 
 }
 
 /**
 Get the resulting demag field through
 Hx = Kxx * Mx + Kxy * My + Kxz * Mz
 Hy = Kxy * Mx + Kyy * My + Kyz * Mz
 Hz = Kxz * Mx + Kyz * My + Kzz * Mz
 
 Noting that K is symmetric, i.e. Kxy = Kyx etc.
 
 @param Mx, My and Mz are n,1 arrays assuming that n == n_K (no checking, so use it correctly!)
 @param Hx, Hy and Hz are the output field vectors (also size n,1)
 @param pref is the factor that goes in front of the multiplication, i.e.
 Hx = pref * ( Kxx * Mx + Kxy * My + Kxz * Mz ) etc.
 */
 void cu_MVMult_GetH( const float* Mx, const float* My, const float* Mz, float* Hx, float* Hy, float* Hz, int* n, float* pref )
 {
	 cublasStatus_t stat;
	 //copy the M vectors to the GPU card
	 stat = cublasSetVector (n_K, sizeof (float), Mx, 1, d_Mx, 1);
	 stat = cublasSetVector (n_K, sizeof (float), My, 1, d_My, 1);
	 stat = cublasSetVector (n_K, sizeof (float), Mz, 1, d_Mz, 1);
	 
	 float beta = 0.;
	 
	 //Kxx * Mx
	 cublasSgemv(handle, CUBLAS_OP_N, n_K, n_K, pref, d_Kxx, n_K, d_Mx, 1, &beta, d_Hx, 1);
	 //Kxy * My
	 beta = 1.; //change beta = 1 so that Hx is updated and not overwritten
	 cublasSgemv(handle, CUBLAS_OP_N, n_K, n_K, pref, d_Kxy, n_K, d_My, 1, &beta, d_Hx, 1);
	 //Kxz * Mz
	 cublasSgemv(handle, CUBLAS_OP_N, n_K, n_K, pref, d_Kxz, n_K, d_Mz, 1, &beta, d_Hx, 1);
	 
	 beta = 0.;
	 //Kyx * Mx
	 cublasSgemv(handle, CUBLAS_OP_N, n_K, n_K, pref, d_Kxy, n_K, d_Mx, 1, &beta, d_Hy, 1);
	 //Kyy * My
	 beta = 1.; //change beta = 1 so that Hy is updated and not overwritten
	 cublasSgemv(handle, CUBLAS_OP_N, n_K, n_K, pref, d_Kyy, n_K, d_My, 1, &beta, d_Hy, 1);
	 //Kyz * Mz
	 cublasSgemv(handle, CUBLAS_OP_N, n_K, n_K, pref, d_Kyz, n_K, d_Mz, 1, &beta, d_Hy, 1);
	 
	 beta = 0.;
	 //Kzx * Mx
	 cublasSgemv(handle, CUBLAS_OP_N, n_K, n_K, pref, d_Kxz, n_K, d_Mx, 1, &beta, d_Hz, 1);
	 //Kzy * My
	 beta = 1.; //change beta = 1 so that Hx is updated and not overwritten
	 cublasSgemv(handle, CUBLAS_OP_N, n_K, n_K, pref, d_Kyz, n_K, d_My, 1, &beta, d_Hz, 1);
	 //Kzz * Mz
	 cublasSgemv(handle, CUBLAS_OP_N, n_K, n_K, pref, d_Kzz, n_K, d_Mz, 1, &beta, d_Hz, 1);
	 
	 //copy the resulting field vector back	 
	 stat = cublasGetVector (n_K, sizeof (float), d_Hx, 1, Hx, 1);
	 stat = cublasGetVector (n_K, sizeof (float), d_Hy, 1, Hy, 1);
	 stat = cublasGetVector (n_K, sizeof (float), d_Hz, 1, Hz, 1);
	 
 }
 
 void cu_test( const float* K, int* n)
 {
	 // FILE *fp;

	//fp = fopen("test.txt", "w+");
	for ( int i=0; i < *n; i++ )
		printf("%f %d\n", K[i], 12 );
	
	//fclose(fp);
 }
 
 void cu_destroy()
 {
	 if ( d_Kxx != NULL )
	 {
		 cudaFree( d_Kxx );
		 d_Kxx = NULL;
		 cudaFree( d_Kxy );
		 d_Kxy = NULL;
		 cudaFree( d_Kxz );
		 d_Kxz = NULL;
		 cudaFree( d_Kyy );
		 d_Kyy = NULL;
		 cudaFree( d_Kyz );
		 d_Kyz = NULL;
		 cudaFree( d_Kzz );
		 d_Kzz = NULL;
	 }
	 
	 if ( d_Hx != NULL )
	 {	 
		 cudaFree( d_Hx );
		 d_Hx = NULL;
		 cudaFree( d_Hy );
		 d_Hy = NULL;
		 cudaFree( d_Hz );
		 d_Hz = NULL;
	 }
	 
	 if ( d_Mx != NULL )
	 {
		 cudaFree( d_Mx );
		 d_Mx = NULL;
		 cudaFree( d_My );
		 d_My = NULL;
		 cudaFree( d_Mz );
		 d_Mz = NULL;
	 }
	 
	 if ( handle != NULL )
	{
		cublasDestroy(handle);
		handle = NULL;
	}
	
		
	freeSparseMatrix( &spKxx );
	freeSparseMatrix( &spKxy );
	freeSparseMatrix( &spKxz );
	freeSparseMatrix( &spKyy );
	freeSparseMatrix( &spKyz );
	freeSparseMatrix( &spKzz );

	if ( sparse_descr != NULL )
	{
		cusparseDestroyMatDescr(sparse_descr);
		sparse_descr = NULL;
	}

	if (sparse_handle != NULL )
	{
		cusparseDestroy(sparse_handle);
		sparse_handle = NULL;
	}

 }



void freeSparseMatrix( CUSparse* mat )
{
	if ( mat->values != NULL )
	{
		cudaFree( mat->values );
		mat->values = NULL;
		cudaFree( mat->cols );
		mat->cols = NULL;
		cudaFree( mat->rows );
		mat->rows = NULL;

	}
}
