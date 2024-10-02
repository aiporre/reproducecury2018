#ifndef SUMCAUCHYGPUKERNEL
#define SUMCAUCHYGPUKERNEL

// Class for the SumCauchy kernel with convolutions performed on the GPU
// SumCauchy kernel is K(x,y)=1/(1+|x-y|^2/sigma^2)

#include "SqDistScalarKernel.h"
#include "CommonFunctions.h"
#include "Utils.h"

////////////////////////////////////
// declarations of Cuda functions //
////////////////////////////////////

template < typename TYPE, int DIMPOINT, int DIMVECT >
int SumCauchyGpuEvalConv1D(int, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int SumCauchyGpuGrad1Conv1D(int, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int SumCauchyGpuGradConv1D(int, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int SumCauchyGpuGradDiffConv1D(int, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int SumCauchyGpuDiffConv1D(int, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int SumCauchyGpuEvalConv2D(int, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int SumCauchyGpuGrad1Conv2D(int, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int SumCauchyGpuGradConv2D(int, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int SumCauchyGpuGradDiffConv2D(int, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int SumCauchyGpuDiffConv2D(int, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, int);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


template < typename TYPE, int DIMPOINT, int DIMVECT >
class SumCauchyGpuKernel : public SqDistScalarKernel< TYPE, DIMPOINT, TinyVector<TYPE,DIMVECT> >
{
    
    using ScalarKernel< TYPE, DIMPOINT, TinyVector<TYPE,DIMVECT> >::Fct;
    using ScalarKernel< TYPE, DIMPOINT, TinyVector<TYPE,DIMVECT> >::testdelete;
    
    typedef TinyVector<TYPE,DIMVECT> Vect;
    typedef Array<Vect,1> ArrVect;
    typedef TinyVector<TYPE,DIMPOINT> Point;
    typedef Array<Point,1> ArrPoint;
    typedef TinyVector<TYPE,DIMPOINT> VectPoint;
    typedef Array<VectPoint,1> ArrVectPoint;
    typedef Array<TYPE,1> ArrType;

    int Nfuns;
    TYPE *Weightsp, *Sigmasp;
    ArrType Weights, Sigmas;
    
    public:
    
    void Write(ofstream &f)
    {
        f << "SumCauchyGpu,Weights=" << endl ;
	WriteArr(Weights,f);
	f << "Sigmas=" << endl;
	WriteArr(Sigmas,f);
    }
    
    SumCauchyGpuKernel(ArrType& weights, ArrType& sigmas)
    {
        Nfuns = weights.rows();
        Weights.resize(Range(1,Nfuns));
        Weights = weights.copy();
        Sigmas.resize(Range(1,Nfuns));
        Sigmas = sigmas.copy();
        Weightsp = Weights.data();
        Sigmasp = Sigmas.data();
    }
    
    ArrVect EvalConv(ArrVectPoint &x, ArrVectPoint &y, ArrVect &beta)
    {
        int nx = x.extent(firstDim);
        int ny = y.extent(firstDim);
        ArrVect gamma(Range(1,nx));
        gamma = 0.0;
        TYPE *xp, *yp, *betap, *gammap, *weightsp, *sigmasp;
        xp = (TYPE*)x.data();
        yp = (TYPE*)y.data();
        betap = (TYPE*)beta.data();
        gammap = (TYPE*)gamma.data();
	SumCauchyGpuEvalConv2D<TYPE,DIMPOINT,DIMVECT>(Nfuns,Weightsp,Sigmasp,xp,yp,betap,gammap,nx,ny);
//SHOW("E")
//SHOW(gamma)
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
        SumCauchyGpuGrad1Conv2D<TYPE,DIMPOINT,DIMVECT>(Nfuns,Weightsp,Sigmasp,alphap,xp,yp,betap,gammap, nx, ny);
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
        SumCauchyGpuGradConv2D<TYPE,DIMPOINT,DIMVECT>(Nfuns,Weightsp,Sigmasp,alphap,xp,betap,gammap, nx);
//SHOW("G")
//SHOW(gamma)
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
		SumCauchyGpuDiffConv2D<TYPE,DIMPOINT,DIMVECT>(Nfuns,Weightsp,Sigmasp,xp,betap,etap,gammap, nx);
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
        SumCauchyGpuGradDiffConv2D<TYPE,DIMPOINT,DIMVECT>(Nfuns,Weightsp,Sigmasp,xp,betap,etap,gammap, nx);
        return gamma;
    }
    
};

#endif // SUMCAUCHYGPUKERNEL
