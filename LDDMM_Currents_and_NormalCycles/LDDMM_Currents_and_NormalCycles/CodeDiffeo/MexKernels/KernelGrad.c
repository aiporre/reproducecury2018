
#include "mex.h"
#include "math.h"
#include "Fun.h"

void KernelGrad(double *G, double *x, double *alpha, double *beta, mwSize n, mwSize dx, mwSize da)
{
    mwSize i, j, k, indi, indj;
    double u, dku, dxij, ab, T;
    for (i=0; i<n; i++)
    {
        for (k=0; k<dx; k++)
            G[k*n+i] = 0.0;
        for (j=0; j<n; j++) 
        {
            u = 0.0;
            for (k=0; k<dx; k++)
            {
                indi = k*n+i;
                indj = k*n+j;
                dxij = x[indi]-x[indj];
                u += dxij*dxij;
            }
            ab = 0;
            for (k=0; k<da; k++)
            {
                indi = k*n+i;
                indj = k*n+j;
                ab += alpha[indi]*beta[indj]+alpha[indj]*beta[indi];
            }
            dku = FUNDIFF(u);
            T = dku * ab;
            for (k=0; k<dx; k++)
                G[k*n+i] += 2 * T * (x[k*n+i]-x[k*n+j]);
        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs!=3)
        mexErrMsgIdAndTxt("KernelGrad:nrhs","Three inputs required.");
    if(nlhs!=1)
        mexErrMsgIdAndTxt("KernelGrad:nlhs","One output required.");
   
    double *x = mxGetPr(prhs[0]);
    size_t nx = mxGetM(prhs[0]);
    size_t dx = mxGetN(prhs[0]);

    double *alpha = mxGetPr(prhs[1]);
    size_t na = mxGetM(prhs[1]);
    size_t da = mxGetN(prhs[1]);

    double *beta = mxGetPr(prhs[2]);
    size_t nb = mxGetM(prhs[2]);
    size_t db = mxGetN(prhs[2]);

    if( nx!=na || nx!=nb || da!=db )
        mexErrMsgIdAndTxt("KernelGrad:inputdims","Dimensions of inputs are not compatible.");
            
    plhs[0] = mxCreateDoubleMatrix((mwSize)nx,(mwSize)dx,mxREAL);
    double *G = mxGetPr(plhs[0]);

    KernelGrad(G,x,alpha,beta,(mwSize)nx,(mwSize)dx,(mwSize)da);
}
