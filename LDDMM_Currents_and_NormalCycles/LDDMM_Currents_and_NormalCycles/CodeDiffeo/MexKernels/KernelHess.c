
#include "mex.h"
#include "math.h"
#include "Fun.h"

void KernelHess(double *H, double *x, double *alpha, double *eta, mwSize nx, mwSize dx, mwSize da)
{
    mwSize i, j, k, indi, indj;
    double a, u, dku, d2ku, diffx, xe, T;
    for (i=0; i<nx; i++) 
    {
        for (k=0; k<dx; k++)
            H[k*nx+i] = 0.0;
        for (j=0; j<nx; j++) 
        {
            u = 0.0;
            xe = 0;
            a = 0;
            for (k=0; k<da; k++)
                a += alpha[k*nx+i]*alpha[k*nx+j];
            for (k=0; k<dx; k++)
            {
                indi = k*nx+i;
                indj = k*nx+j;
                diffx = x[indi]-x[indj];
                u += diffx*diffx;
                xe += diffx*(eta[indi]-eta[indj]);
            }
            dku = FUNDIFF(u);
            d2ku = FUNDIFF2(u);
            T = 2 * d2ku * xe;
            for (k=0; k<dx; k++)
	    {
		indi = k*nx+i;
                indj = k*nx+j;
                H[indi] += a * (T * (x[indi]-x[indj]) + dku * (eta[indi]-eta[indj]));
	    }
        }
        for (k=0; k<dx; k++)
            H[k*nx+i] *= 4;
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs!=3)
        mexErrMsgIdAndTxt("KernelHess:nrhs","Three inputs required.");
    if(nlhs!=1)
        mexErrMsgIdAndTxt("KernelHess:nlhs","One output required.");
   
    double *x = mxGetPr(prhs[0]);
    size_t nx = mxGetM(prhs[0]);
    size_t dx = mxGetN(prhs[0]);

    double *alpha = mxGetPr(prhs[1]);
    size_t na = mxGetM(prhs[1]);
    size_t da = mxGetN(prhs[1]);

    double *eta = mxGetPr(prhs[2]);
    size_t ne = mxGetM(prhs[2]);
    size_t de = mxGetN(prhs[2]);

    if( nx!=na || nx!=ne || dx!=de )
        mexErrMsgIdAndTxt("KernelHess:inputdims","Dimensions of inputs are not compatible.");
            
    plhs[0] = mxCreateDoubleMatrix((mwSize)nx,(mwSize)dx,mxREAL);
    double *H = mxGetPr(plhs[0]);

    KernelHess(H,x,alpha,eta,(mwSize)nx,(mwSize)dx,(mwSize)da);
}

