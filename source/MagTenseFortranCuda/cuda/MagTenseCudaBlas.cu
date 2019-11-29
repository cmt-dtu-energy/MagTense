
#include <stdlib.h>
#include <stdio.h>
#include "cublas_v2.h"



 
 /**
 Demag matrices to be stored in GPU memory
 */
 float* d_Kxx;
 float* d_Kxy;
 float* d_Kxz;
 float* d_Kyy;
 float* d_Kyz;
 float* d_Kzz;
 
 //Magnetization vectors on GPU
 float* d_Mx;
 float* d_My;
 float* d_Mz;
 
 //Magnetic field vectors on GPU
 float* d_Hx;
 float* d_Hy;
 float* d_Hz;
 
 //size of the matrices (n_K x n_K)
 int n_K;
 
 //general handle for cuBlas (initalized once)
 cublasHandle_t handle;
 
 /**
 Kaspar K. Nielsen, kaki@dtu.dk, DTU 2019
 Method for transfering the demagnetization matrices from Fortran
 @param Kxx, Kxy, Kxz, Kyy, Kyz and Kzz are 2D arrays of size n x n implemented as an array of n pointers each pointing to an array of n values
 @param n is the size of the square matrix
 */
 void cu_initDemagMatrices( const float* Kxx, const float* Kxy, const float* Kxz, const float* Kyy, const float* Kyz, const float* Kzz, int* n)
 {
	 FILE *fp;

	fp = fopen("test.txt", "w+");
	fprintf(fp,"%f %d\n", Kxx[0], *n );
	//fprintf(fp, "This is testing for fprintf...\n");
	//fputs("This is testing for fputs...\n", fp);
	fclose(fp);
	 
	 
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
	 cudaFree( d_Kxx );
	 cudaFree( d_Kxy );
	 cudaFree( d_Kxz );
	 cudaFree( d_Kyy );
	 cudaFree( d_Kyz );
	 cudaFree( d_Kzz );
	 
	 cudaFree( d_Hx );
	 cudaFree( d_Hy );
	 cudaFree( d_Hz );
	 
	 cudaFree( d_Mx );
	 cudaFree( d_My );
	 cudaFree( d_Mz );
	 
	 cublasDestroy(handle);
 }