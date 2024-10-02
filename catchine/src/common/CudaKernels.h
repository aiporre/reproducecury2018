#include "CudaFunctions.h"

template < typename TYPE, int DIMPOINT, int DIMVECT, class RADIAL_FUN >
class ScalarRadialKernel
{
    
    RADIAL_FUN Rfun;
    
public:
    
    ScalarRadialKernel(RADIAL_FUN rfun)
    {
        Rfun = rfun;
    }
    
    __device__ __forceinline__ void Eval(TYPE* gammai, TYPE* xi, TYPE* yj, TYPE* betaj)
    {
        TYPE r2 = 0.0f;
        TYPE temp;
        for(int k=0; k<DIMPOINT; k++)
        {
            temp =  yj[k]-xi[k];
            r2 += temp*temp;
        }
        TYPE s = Rfun.Eval(r2);
        for(int k=0; k<DIMVECT; k++)
            gammai[k] += s * betaj[k];
    }
    
    __device__ __forceinline__ void Grad1(TYPE* gammai, TYPE* alphai, TYPE* xi, TYPE* yj, TYPE* betaj)
    {
        TYPE r2 = 0.0f, sga = 0.0f;
        TYPE xmy[DIMPOINT];
        for(int k=0; k<DIMPOINT; k++)
        {
            xmy[k] =  xi[k]-yj[k];
            r2 += xmy[k]*xmy[k];
        }
        for(int k=0; k<DIMVECT; k++)
            sga += betaj[k]*alphai[k];
        TYPE s = 2.0 * sga * Rfun.Diff(r2);
        for(int k=0; k<DIMPOINT; k++)
            gammai[k] += s * xmy[k];
    }
    
    __device__ __forceinline__ void Grad(TYPE* gammai, TYPE* xi, TYPE* xj, TYPE* alphai, TYPE* alphaj, TYPE* betai, TYPE* betaj)
    {
        TYPE r2 = 0.0f, sga = 0.0f;
        TYPE ximxj[DIMPOINT];
        for(int k=0; k<DIMPOINT; k++)
        {
            ximxj[k] =  xi[k]-xj[k];
            r2 += ximxj[k]*ximxj[k];
        }
        for(int k=0; k<DIMVECT; k++)
            sga += betaj[k]*alphai[k] + betai[k]*alphaj[k];
        TYPE s = 2.0 * sga * Rfun.Diff(r2);
        for(int k=0; k<DIMPOINT; k++)
            gammai[k] += s * ximxj[k];
    }
    
    __device__ __forceinline__ void GradDiff(TYPE* gammai, TYPE* xi, TYPE* xj, TYPE* betai, TYPE* betaj, TYPE* etai, TYPE* etaj)
    {
        TYPE r2 = 0.0f, bidbj = 0.0f, dotex = 0.0f;
        TYPE ximxj[DIMPOINT], eimej[DIMPOINT];
        for(int k=0; k<DIMPOINT; k++)
        {
            ximxj[k] =  xi[k]-xj[k];
            r2 += ximxj[k]*ximxj[k];
            eimej[k] =  etai[k]-etaj[k];
            dotex += ximxj[k]*eimej[k];
        }
        for(int k=0; k<DIMVECT; k++)
            bidbj += betai[k]*betaj[k];
        TYPE d1, d2;
        Rfun.DiffDiff2(r2,&d1,&d2);
        d1 *= 4 * bidbj;
        d2 *= 8 * bidbj * dotex;
        for(int k=0; k<DIMPOINT; k++)
            gammai[k] += d1 * eimej[k] + d2 * ximxj[k];
    }
    
    __device__ __forceinline__ void Diff(TYPE* gammai, TYPE* xi, TYPE* xj, TYPE* betaj, TYPE* etai, TYPE* etaj)
    {
        TYPE r2 = 0.0f, dotex = 0.0f;
        TYPE ximxj[DIMPOINT], eimej[DIMPOINT];
        for(int k=0; k<DIMPOINT; k++)
        {
            ximxj[k] =  xi[k]-xj[k];
            r2 += ximxj[k]*ximxj[k];
            eimej[k] =  etai[k]-etaj[k];
            dotex += ximxj[k]*eimej[k];
        }
        TYPE s = Rfun.Diff(r2) * 2.0 * dotex;
        for(int k=0; k<DIMVECT; k++)
            gammai[k] += s * betaj[k];
    }
};
