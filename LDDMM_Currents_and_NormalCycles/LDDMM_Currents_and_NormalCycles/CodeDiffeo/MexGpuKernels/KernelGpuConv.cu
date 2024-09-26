#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#include <mex.h>
#include "GpuConv2D.cu"
#include "CudaKernels.h"
#include "CudaNCSurfKernels.h"
#include "CudaVarSurfKernels.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    typedef RADIALFUN<__TYPE__ > RadialFun;

#define SCALARRADIAL 1
#if KERNEL==SCALARRADIAL
typedef ScalarRadialKernel<__TYPE__,__DIMPOINT__,__DIMVECT__,RadialFun> KER;
#endif

#define NCSURF 2
#if KERNEL==NCSURF
typedef NCSurfKernel<__TYPE__,RadialFun> KER;
#endif

#define VARSURF 3
#if KERNEL==VARSURF
typedef VarSurfKernel<__TYPE__,RadialFun> KER;
#endif

struct KER::EVAL funeval;

    const int DIMX1 = KER::EVAL::DIMSX::FIRST;
    typedef typename KER::EVAL::DIMSX DIMSX;   
    typedef typename KER::EVAL::DIMSY DIMSY;
    const int SIZEX = DIMSX::SIZE;
    const int SIZEY = DIMSY::SIZE;        

    if(nrhs<SIZEX+SIZEY-1)
        mexErrMsgIdAndTxt("KernelGpuConv:nrhs","At least %d inputs required.",SIZEX+SIZEY-1);

    if(nrhs>SIZEX+SIZEY)
        mexErrMsgIdAndTxt("KernelGpuConv:nrhs","Too many inputs.");

    double *x[SIZEX];
    int nx[SIZEX];

    double *y[SIZEY];
    int ny[SIZEY];

    for(int k=1; k<SIZEX; k++)
    {
        x[k] = mxGetPr(prhs[k-1]);
        if(mxGetM(prhs[k-1])!=DIMSX::VAL(k))
            mexErrMsgIdAndTxt("KernelGpuConv:inputdims","Dimensions of inputs are not compatible.");
        nx[k] = mxGetN(prhs[k-1]);
        if(nx[k]!=nx[1])
            mexErrMsgIdAndTxt("KernelGpuConv:inputdims","Dimensions of inputs are not compatible.");
    }

    for(int k=0; k<SIZEY; k++)
    {
        y[k] = mxGetPr(prhs[SIZEX+k-1]);
        if(mxGetM(prhs[SIZEX+k-1])!=DIMSY::VAL(k))
            mexErrMsgIdAndTxt("KernelGpuConv:inputdims","Dimensions of inputs are not compatible.");
        ny[k] = mxGetN(prhs[SIZEX+k-1]);
        if(ny[k]!=ny[0])
            mexErrMsgIdAndTxt("KernelGpuConv:inputdims","Dimensions of inputs are not compatible.");
    }

    nx[0] = nx[1];
    plhs[0] = mxCreateDoubleMatrix((mwSize)DIMX1,(mwSize)nx[0],mxREAL);
    x[0] = mxGetPr(plhs[0]);

    if(nrhs==SIZEX+SIZEY)
    {
        int deviceID = *mxGetPr(prhs[SIZEX+SIZEY-1]);
        cudaSetDevice(deviceID);
    }

	GpuConv2D(KER(RadialFun()),funeval,nx[0],ny[0],x,y);

}

