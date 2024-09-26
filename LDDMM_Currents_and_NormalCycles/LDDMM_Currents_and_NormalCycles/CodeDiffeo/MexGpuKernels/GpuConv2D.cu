#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <cuda.h>
#include "Pack.h"


using namespace std;

template <typename TYPE, int DIMVECT>
__global__ void reduce0(TYPE* in, TYPE* out, int sizeY,int nx)
{
	TYPE res = 0;
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if(tid < nx*DIMVECT)
    {
		for (int i = 0; i < sizeY; i++) 
            res += in[tid + i*nx*DIMVECT];
		/*res = in[tid+ nx* DIMVECT];*/
		out[tid] = res;
	}
}








// thread kernel: computation of x1i = sum_j k(x2i,x3i,...,y1j,y2j,...) for index i given by thread id.
template < typename TYPE, class KER, class FUN >
__global__ void GpuConv2DOnDevice(KER Ker, FUN fun, int nx, int ny, TYPE** px, TYPE** py)
{

    typedef typename FUN::DIMSX DIMSX;
    typedef typename FUN::DIMSY DIMSY;
    const int DIMX = DIMSX::SUM;
    const int DIMY = DIMSY::SUM;        
    const int DIMX1 = DIMSX::FIRST;

    extern __shared__ char yj_char[];
    TYPE* const yj = reinterpret_cast<TYPE*>(yj_char);

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    TYPE xi[DIMX];
    if(i<nx)  // we will compute x1i only if i is in the range
    {
        for(int k=0; k<DIMX1; k++)
            xi[k] = 0.0f;
        // load xi from device global memory
	DIMSX::NEXT::load(i,xi+DIMX1,px+1);
    }
    
    int j = blockIdx.y * blockDim.x + threadIdx.x;
    if(j<ny) // we load yj from device global memory only if j<ny
	DIMSY::load(j,yj+threadIdx.x*DIMY,py);    	
    __syncthreads();
        
    if(i<nx) // we compute x1i only if needed
    {
    	TYPE* yjrel = yj;
        for(int jrel = 0; (jrel<blockDim.x) && ((blockDim.x*blockIdx.y+jrel)< ny); jrel++, yjrel+=DIMY)
		DIMSX::template call<DIMSY>(fun,xi,yjrel,Ker);	
        __syncthreads();
    }

    // Save the result in global memory.
    if(i<nx)
        for(int k=0; k<DIMX1; k++)
            (*px)[blockIdx.y*DIMX1*nx+i*DIMX1+k] = xi[k];
}
///////////////////////////////////////////////////


template < typename TYPE, class KER, class FUN >
int GpuConv2D(KER Ker, FUN fun, int nx, int ny, TYPE** px_h, TYPE** py_h)
{
    typedef typename FUN::DIMSX DIMSX;
    typedef typename FUN::DIMSY DIMSY;
    const int DIMX = DIMSX::SUM;
    const int DIMY = DIMSY::SUM;
    const int DIMX1 = DIMSX::FIRST;
    const int SIZEX = DIMSX::SIZE;
    const int SIZEY = DIMSY::SIZE;

    // Data on the device.
    TYPE *x1B, *x_d, *y_d, **px_d, **py_d;

    cudaHostAlloc((void**)&px_d, SIZEX*sizeof(TYPE*), cudaHostAllocMapped);
    cudaHostAlloc((void**)&py_d, SIZEY*sizeof(TYPE*), cudaHostAllocMapped);

    // Allocate arrays on device.
    cudaMalloc((void**)&x_d, sizeof(TYPE)*(nx*DIMX));
    cudaMalloc((void**)&y_d, sizeof(TYPE)*(ny*DIMY));

    // Send data from host to device.

    int nvals;
    px_d[0] = x_d;
    nvals = nx*DIMSX::VAL(0);
    for(int k=1; k<SIZEX; k++)
    {
        px_d[k] = px_d[k-1] + nvals;
        nvals = nx*DIMSX::VAL(k);
        cudaMemcpy(px_d[k], px_h[k], sizeof(TYPE)*nvals, cudaMemcpyHostToDevice);
    }
    py_d[0] = y_d;
    nvals = ny*DIMSY::VAL(0);
    cudaMemcpy(py_d[0], py_h[0], sizeof(TYPE)*nvals, cudaMemcpyHostToDevice);
    for(int k=1; k<SIZEY; k++)
    {
        py_d[k] = py_d[k-1] + nvals;
        nvals = ny*DIMSY::VAL(k);
        cudaMemcpy(py_d[k], py_h[k], sizeof(TYPE)*nvals, cudaMemcpyHostToDevice);
    }

    // Compute on device.
    dim3 blockSize;
    blockSize.x = 192; // number of threads in each block
    int blockSizey = blockSize.x;
    dim3 gridSize;
    gridSize.x =  nx / blockSize.x + (nx%blockSize.x==0 ? 0 : 1);
	gridSize.y =  ny / blockSizey + (ny%blockSizey==0 ? 0 : 1);

    // Reduce  : grid and block are 1d
    dim3 blockSize2;
    blockSize2.x = 192; // number of threads in each block
    dim3 gridSize2;
    gridSize2.x =  (nx*DIMX1) / blockSize2.x + ((nx*DIMX1)%blockSize2.x==0 ? 0 : 1);

    cudaMalloc((void**)&x1B, sizeof(TYPE)*(nx*DIMX1*gridSize.y));
    px_d[0] = x1B;

    GpuConv2DOnDevice<TYPE><<<gridSize,blockSize,blockSize.x*(DIMY)*sizeof(TYPE)>>>(Ker,fun,nx,ny,px_d,py_d);
	

    reduce0<TYPE,DIMX1><<<gridSize2, blockSize2>>>(x1B, x_d, gridSize.y,nx);
    
    // block until the device has completed
    cudaThreadSynchronize();

    // Send data from device to host.
    cudaMemcpy(*px_h, x_d, sizeof(TYPE)*(nx*DIMX1),cudaMemcpyDeviceToHost);

    // Free memory.
    cudaFree(x_d);
    cudaFree(y_d);
    cudaFree(x1B);
    cudaFree(px_d);
    cudaFree(py_d);

    return 0;
}


template < typename TYPE, class KER, class FUN, typename... Args >
int GpuConv2D(KER Ker, FUN fun, int nx, int ny, TYPE* x1_h, Args... args)
{

    typedef typename FUN::DIMSX DIMSX;
    typedef typename FUN::DIMSY DIMSY;
    const int SIZEX = DIMSX::SIZE;
    const int SIZEY = DIMSY::SIZE;

    TYPE *px_h[SIZEX];
    TYPE *py_h[SIZEY];
    DIMSX::getlist(px_h,x1_h,args...);
    DIMSX::template getlist<DIMSY>(py_h,x1_h,args...);
 
	return GpuConv2D(Ker,fun,nx,ny,x1_h,px_h,py_h);

}



