#ifndef GAUSSGPUKERNEL
#define GAUSSGPUKERNEL

// Class for the Gauss kernel with convolutions performed on the GPU
// Gauss kernel is K(x,y)=1/(1+|x-y|^2/sigma^2)

#include "SqDistScalarKernel.h"
#include "CommonFunctions.h"
#include "Utils.h"

////////////////////////////////////
// declarations of Cuda functions //
////////////////////////////////////

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuEvalConv1D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGrad1Conv1D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGradConv1D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGradDiffConv1D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuDiffConv1D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuEvalConv2D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGrad1Conv2D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGradConv2D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGradDiffConv2D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuDiffConv2D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


template < typename TYPE, int DIMPOINT, int DIMVECT >
class GaussGpuKernel : public SqDistScalarKernel< TYPE, DIMPOINT, TinyVector<TYPE,DIMVECT> >
{
    
    using ScalarKernel< TYPE, DIMPOINT, TinyVector<TYPE,DIMVECT> >::Fct;
    using ScalarKernel< TYPE, DIMPOINT, TinyVector<TYPE,DIMVECT> >::testdelete;
    
    typedef TinyVector<TYPE,DIMVECT> Vect;
    typedef Array<Vect,1> ArrVect;
    typedef TinyVector<TYPE,DIMPOINT> Point;
    typedef Array<Point,1> ArrPoint;
    typedef TinyVector<TYPE,DIMPOINT> VectPoint;
    typedef Array<VectPoint,1> ArrVectPoint;
    
    TYPE Sigma, ooSigma2;
    
    public:
    
    void Write(ofstream &f)
    {
        f << "GaussGpu,sigma=" << endl << Sigma << endl;
    }
    
    GaussGpuKernel(TYPE sigma)
    {
        Sigma = sigma;
        Fct = new GaussFunction<TYPE>(Sigma);
        testdelete = true;
        ooSigma2 = 1.0/pow2(Sigma);
    }
    
    ArrVect EvalConv(ArrVectPoint &x, ArrVectPoint &y, ArrVect &beta)
    {
        int nx = x.extent(firstDim);
        int ny = y.extent(firstDim);
        ArrVect gamma(Range(1,nx));
        gamma = 0.0;
        TYPE *xp, *yp, *betap, *gammap;
        xp = (TYPE*)x.data();
        yp = (TYPE*)y.data();
        betap = (TYPE*)beta.data();
        gammap = (TYPE*)gamma.data();
        GaussGpuEvalConv2D<TYPE,DIMPOINT,DIMVECT>(Sigma,xp,yp,betap,gammap,nx,ny);
        return gamma;
    }



    
    ArrVectPoint Grad1Conv(ArrVect &alpha, ArrPoint &x, ArrPoint &y, ArrVect &beta)
    {
        int nx = x.extent(firstDim);
        int ny = y.extent(firstDim);
        ArrVectPoint gamma(Range(1,nx));
        TYPE *xp, *yp, *alphap, *gammap, *betap;
        xp = (TYPE*)x.data();
        yp = (TYPE*)y.data();
        alphap = (TYPE*)alpha.data();
        gammap = (TYPE*)gamma.data();
        betap = (TYPE*)beta.data();
        GaussGpuGrad1Conv2D<TYPE,DIMPOINT,DIMVECT>(Sigma,alphap,xp,yp,betap,gammap, nx, ny);
        return gamma;
    }
    
    ArrVectPoint GradConv(ArrVect &alpha, ArrPoint &x, ArrVect &beta)
    {
        int nx = x.extent(firstDim);
        ArrVectPoint gamma(Range(1,nx));
        TYPE *xp, *alphap, *gammap, *betap;
        xp = (TYPE*)x.data();
        alphap = (TYPE*)alpha.data();
        gammap = (TYPE*)gamma.data();
        betap = (TYPE*)beta.data();
        GaussGpuGradConv2D<TYPE,DIMPOINT,DIMVECT>(Sigma,alphap,xp,betap,gammap, nx);
        return gamma;
    }
	
	ArrVect DiffConv(ArrPoint &x, ArrVect &beta, ArrVectPoint &eta)
	{
		int nx = x.extent(firstDim);
		ArrVect gamma(Range(1,nx));
		TYPE *xp, *etap, *gammap, *betap;
		xp = (TYPE*)x.data();
		etap = (TYPE*)eta.data();
		gammap = (TYPE*)gamma.data();
		betap = (TYPE*)beta.data();
		GaussGpuDiffConv2D<TYPE,DIMPOINT,DIMVECT>(Sigma,xp,betap,etap,gammap, nx);
		return gamma;
	}
	
    
    ArrVectPoint GradDiffConv(ArrPoint &x, ArrVect &beta, ArrVectPoint &eta)
    {
        int nx = x.extent(firstDim);
        ArrVectPoint gamma(Range(1,nx));
        TYPE *xp, *etap, *gammap, *betap;
        xp = (TYPE*)x.data();
        etap = (TYPE*)eta.data();
        gammap = (TYPE*)gamma.data();
        betap = (TYPE*)beta.data();
        GaussGpuGradDiffConv2D<TYPE,DIMPOINT,DIMVECT>(Sigma,xp,betap,etap,gammap, nx);
        return gamma;
    }
    
};

#endif // GAUSSGPUKERNEL
