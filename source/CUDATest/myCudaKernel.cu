#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
// CUDA kernel. Each thread takes care of one element of c
__global__ void vecAdd(float *a, float *b, float *c, int n)
{
    // Get our global thread ID
    int id = blockIdx.x*blockDim.x+threadIdx.x;
 
    // Make sure we do not go out of bounds
    if (id < n)
        c[id] = a[id] + b[id];
}
 
float* internalBuffer;
int nBuf;

void cu_initBuf( float* buf, int* n)
{
	nBuf = *n;
	
	internalBuffer = (float*)malloc(nBuf*sizeof(float));
	
	for ( int i = 0; i < nBuf; i++ )
	{
		internalBuffer[i] = buf[i];
	}
}
 
 void cu_destroy()
 {
	free(internalBuffer);
 }
 
void cu_vecAdd_wrapper(float* val, int* nn)
{
	for ( int uu=0; uu < *nn; uu++ )
	{
		val[uu] = internalBuffer[uu] * 3;
	}
		
    // Size of vectors
    int n = 100000;
 
    // Host input vectors
    float *h_a;
    float *h_b;
    //Host output vector
    float *h_c;
 
    // Device input vectors
    float *d_a;
    float *d_b;
    //Device output vector
    float *d_c;
 
    // Size, in bytes, of each vector
    size_t bytes = n*sizeof(float);
 
    // Allocate memory for each vector on host
    h_a = (float*)malloc(bytes);
    h_b = (float*)malloc(bytes);
    h_c = (float*)malloc(bytes);
 
    // Allocate memory for each vector on GPU
    cudaMalloc(&d_a, bytes);
    cudaMalloc(&d_b, bytes);
    cudaMalloc(&d_c, bytes);
 
    int i;
    // Initialize vectors on host
    for( i = 0; i < n; i++ ) {
        h_a[i] = sinf(i)*sinf(i);
        h_b[i] = cosf(i)*cosf(i);
    }
 
    // Copy host vectors to device
    cudaMemcpy( d_a, h_a, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy( d_b, h_b, bytes, cudaMemcpyHostToDevice);
 
    int blockSize, gridSize;
 
    // Number of threads in each thread block
    blockSize = 1024;
 
    // Number of thread blocks in grid
    gridSize = (int)ceil((float)n/blockSize);
 
    // Execute the kernel
    vecAdd<<<gridSize, blockSize>>>(d_a, d_b, d_c, n);
 
    // Copy array back to host
    cudaMemcpy( h_c, d_c, bytes, cudaMemcpyDeviceToHost );
 
    // Sum up vector c and print result divided by n, this should equal 1 within error
    float sum = 0;
    for(i=0; i<n; i++)
        sum += h_c[i];
    printf("final result: %f\n", sum/n);
 	
 
    // Release device memory
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);
 
    // Release host memory
    free(h_a);
    free(h_b);
    free(h_c);
	

}