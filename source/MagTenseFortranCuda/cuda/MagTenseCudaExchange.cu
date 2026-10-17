// Exchange field on the GPU.
//
// The exchange field is a sparse matrix-vector product H_j = diag(s) * A_exch * m, applied to the
// three components of m. On the CPU this is three MKL products, which inside the OpenMP task that
// computeHeff runs it in are sequential, and which read the matrix three times. Here the row
// scaling s = 2*Jfact is folded into the values once at initialization, and one kernel reads each
// row once and produces all three components. The arithmetic is FP64, as on the CPU: the exchange
// field is the Laplacian of a nearly uniform m, which loses too many digits in FP32.
//
// The kernel runs on its own non-blocking stream, because computeHeff calls it from an OpenMP task
// while the main thread may be running the dip-fmm demag field on the same GPU.
//
// Set MAGTENSE_EXCHANGE_GPU=0 to keep the exchange field on the CPU in a CUDA run.

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {

struct ExchangeGpu {
    int n = 0;
    int nnz = 0;
    int *rowptr = nullptr;
    int *cols = nullptr;
    double *vals = nullptr;
    double *M = nullptr;    // [mx; my; mz], 3n
    double *H = nullptr;    // [hx; hy; hz], 3n
    cudaStream_t stream = nullptr;
    bool ready = false;
};

ExchangeGpu g_exch;

bool check( cudaError_t err, const char* what )
{
    if ( err == cudaSuccess ) return true;
    fprintf( stderr, "MagTense CUDA exchange: %s failed: %s\n", what, cudaGetErrorString(err) );
    return false;
}

// One thread per row; a row of the exchange operator holds the cell and its neighbours, i.e. 7
// entries on a uniform grid and a few tens on an unstructured mesh, so a row is short enough that
// one thread is enough to keep the memory system busy.
__global__ void exchangeCsr3( int n, const int* __restrict__ rowptr, const int* __restrict__ cols,
                              const double* __restrict__ vals, const double* __restrict__ M,
                              double* __restrict__ H )
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if ( i >= n ) return;
    double hx = 0.0, hy = 0.0, hz = 0.0;
    const int k_end = rowptr[i+1];
    for ( int k = rowptr[i]; k < k_end; ++k ) {
        const int c = cols[k];
        const double v = vals[k];
        hx += v * M[c];
        hy += v * M[n + c];
        hz += v * M[2*n + c];
    }
    H[i] = hx;
    H[n + i] = hy;
    H[2*n + i] = hz;
}

} // namespace

extern "C" {

void cu_exchange_destroy()
{
    if ( g_exch.stream ) cudaStreamSynchronize( g_exch.stream );
    cudaFree( g_exch.rowptr );
    cudaFree( g_exch.cols );
    cudaFree( g_exch.vals );
    cudaFree( g_exch.M );
    cudaFree( g_exch.H );
    if ( g_exch.stream ) cudaStreamDestroy( g_exch.stream );
    g_exch = ExchangeGpu();
}

// Uploads the exchange operator in zero-based three-array CSR form, with the row scaling already
// applied to the values. Returns 0 on success, 1 when disabled by MAGTENSE_EXCHANGE_GPU=0 and
// -1 on a CUDA error, in which case nothing stays allocated and the caller uses the CPU.
int cu_exchange_init( const int* n, const int* nnz, const int* rowptr, const int* cols, const double* vals )
{
    const char* env = getenv( "MAGTENSE_EXCHANGE_GPU" );
    if ( env && strcmp( env, "0" ) == 0 ) return 1;

    cu_exchange_destroy();
    g_exch.n = *n;
    g_exch.nnz = *nnz;
    const size_t n3 = 3 * (size_t)g_exch.n;
    bool ok = check( cudaStreamCreateWithFlags( &g_exch.stream, cudaStreamNonBlocking ), "cudaStreamCreate" )
        && check( cudaMalloc( &g_exch.rowptr, (g_exch.n + 1) * sizeof(int) ), "cudaMalloc rowptr" )
        && check( cudaMalloc( &g_exch.cols, g_exch.nnz * sizeof(int) ), "cudaMalloc cols" )
        && check( cudaMalloc( &g_exch.vals, g_exch.nnz * sizeof(double) ), "cudaMalloc vals" )
        && check( cudaMalloc( &g_exch.M, n3 * sizeof(double) ), "cudaMalloc M" )
        && check( cudaMalloc( &g_exch.H, n3 * sizeof(double) ), "cudaMalloc H" )
        && check( cudaMemcpy( g_exch.rowptr, rowptr, (g_exch.n + 1) * sizeof(int), cudaMemcpyHostToDevice ), "copy rowptr" )
        && check( cudaMemcpy( g_exch.cols, cols, g_exch.nnz * sizeof(int), cudaMemcpyHostToDevice ), "copy cols" )
        && check( cudaMemcpy( g_exch.vals, vals, g_exch.nnz * sizeof(double), cudaMemcpyHostToDevice ), "copy vals" );
    if ( !ok ) {
        cu_exchange_destroy();
        return -1;
    }
    g_exch.ready = true;
    return 0;
}

// H = diag(s) * A_exch * M for the three components. Returns 0 on success and -1 on a CUDA error.
int cu_exchange_apply( const double* Mx, const double* My, const double* Mz, double* Hx, double* Hy, double* Hz )
{
    if ( !g_exch.ready ) return -1;
    const int n = g_exch.n;
    const size_t bytes = n * sizeof(double);
    cudaStream_t s = g_exch.stream;
    cudaMemcpyAsync( g_exch.M,         Mx, bytes, cudaMemcpyHostToDevice, s );
    cudaMemcpyAsync( g_exch.M + n,     My, bytes, cudaMemcpyHostToDevice, s );
    cudaMemcpyAsync( g_exch.M + 2 * n, Mz, bytes, cudaMemcpyHostToDevice, s );
    const int block = 256;
    exchangeCsr3<<< (n + block - 1) / block, block, 0, s >>>( n, g_exch.rowptr, g_exch.cols, g_exch.vals, g_exch.M, g_exch.H );
    // A launch failure (e.g. no kernel image for the device) is reported here and not by the
    // synchronization below
    if ( !check( cudaGetLastError(), "exchange kernel launch" ) ) return -1;
    cudaMemcpyAsync( Hx, g_exch.H,         bytes, cudaMemcpyDeviceToHost, s );
    cudaMemcpyAsync( Hy, g_exch.H + n,     bytes, cudaMemcpyDeviceToHost, s );
    cudaMemcpyAsync( Hz, g_exch.H + 2 * n, bytes, cudaMemcpyDeviceToHost, s );
    return check( cudaStreamSynchronize( s ), "exchange product" ) ? 0 : -1;
}

}
