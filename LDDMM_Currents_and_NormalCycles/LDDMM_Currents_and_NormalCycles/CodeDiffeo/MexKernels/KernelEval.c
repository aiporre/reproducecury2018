#include "mex.h"
#include "math.h"
#include "Fun.h"

void KernelEval(double *gamma, double *y, double *x, double *alpha, mwSize nx, mwSize ny, mwSize dx, mwSize da)
{
    mwSize i, j, k, indi;
    double u, ku, dxy;
    for (i=0; i<ny; i++)
    {
        for (k=0; k<da; k++)
            gamma[k*ny+i] = 0.0;
        for (j=0; j<nx; j++)
        {
            u = 0.0;
            for (k=0; k<dx; k++)
            {
                dxy = y[k*ny+i]-x[k*nx+j];
                u += dxy*dxy;
            }
            ku = FUNEVAL(u);
            for (k=0; k<da; k++)
                gamma[k*ny+i] += ku * alpha[k*nx+j];
        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs!=3)
        mexErrMsgIdAndTxt("KernelEval:nrhs","Three inputs required.");

    double *y = mxGetPr(prhs[0]);
    size_t ny = mxGetM(prhs[0]);
    size_t dy = mxGetN(prhs[0]);

    double *x = mxGetPr(prhs[1]); 
    size_t nx = mxGetM(prhs[1]);
    size_t dx = mxGetN(prhs[1]);
 
    double *alpha = mxGetPr(prhs[2]);  
    size_t na = mxGetM(prhs[2]);
    size_t da = mxGetN(prhs[2]);

     if(dy!=dx || nx!=na)
        mexErrMsgIdAndTxt("KernelEval:inputdims","Dimensions of inputs are not compatible.");

    if(nlhs!=1)
        mexErrMsgIdAndTxt("KernelEval:nlhs","One output required.");

    plhs[0] = mxCreateDoubleMatrix((mwSize)ny,(mwSize)da,mxREAL);
    double *gamma = mxGetPr(plhs[0]);

    KernelEval(gamma,y,x,alpha,(mwSize)nx,(mwSize)ny,(mwSize)dx,(mwSize)da);
}
