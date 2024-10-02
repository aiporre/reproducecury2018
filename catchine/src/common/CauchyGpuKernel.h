#ifndef CAUCHYGPUKERNEL
#define CAUCHYGPUKERNEL

// Class for the Cauchy kernel with convolutions performed on the GPU
// Cauchy kernel is K(x,y)=1/(1+|x-y|^2/sigma^2)

#include "SqDistScalarKernel.h"
#include "CommonFunctions.h"
#include "Utils.h"

////////////////////////////////////
// declarations of Cuda functions //
////////////////////////////////////

template < typename TYPE, int DIMPOINT, int DIMVECT >
int CauchyGpuEvalConv1D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int CauchyGpuGrad1Conv1D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int CauchyGpuGradConv1D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int CauchyGpuGradDiffConv1D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int CauchyGpuDiffConv1D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);




template < typename TYPE, int DIMPOINT, int DIMVECT >
int CauchyGpuEvalConv2D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int CauchyGpuGrad1Conv2D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int CauchyGpuGradConv2D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int CauchyGpuGradDiffConv2D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int CauchyGpuDiffConv2D(TYPE, TYPE*, TYPE*, TYPE*, TYPE*, int);



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


template < typename TYPE, int DIMPOINT, int DIMVECT >
class CauchyGpuKernel : public SqDistScalarKernel< TYPE, DIMPOINT, TinyVector<TYPE,DIMVECT> >
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
            f << "CauchyGpu,sigma=" << endl << Sigma << endl;
        }

        CauchyGpuKernel(TYPE sigma)
        {
            Sigma = sigma;
            Fct = new CauchyFunction<TYPE>(Sigma);
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
            CauchyGpuEvalConv2D<TYPE,DIMPOINT,DIMVECT>(Sigma,xp,yp,betap,gammap,nx,ny);
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
            CauchyGpuGrad1Conv2D<TYPE,DIMPOINT,DIMVECT>(Sigma,alphap,xp,yp,betap,gammap, nx, ny);
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
            CauchyGpuGradConv2D<TYPE,DIMPOINT,DIMVECT>(Sigma,alphap,xp,betap,gammap, nx);
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
		CauchyGpuDiffConv2D<TYPE,DIMPOINT,DIMVECT>(Sigma,xp,betap,etap,gammap, nx);
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
        CauchyGpuGradDiffConv2D<TYPE,DIMPOINT,DIMVECT>(Sigma,xp,betap,etap,gammap, nx);
        return gamma;
    }
			
};

#endif // CAUCHYGPUKERNEL
