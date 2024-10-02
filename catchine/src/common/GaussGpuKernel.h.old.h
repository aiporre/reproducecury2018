#ifndef GAUSSGPUKERNEL
#define GAUSSGPUKERNEL

#include "SqDistScalarKernel.h"
#include "CommonFunctions.h"
#include "Utils.h"

////////////////////////////
// wrappers for Cuda code //
////////////////////////////

// for EvalConv

extern "C" int GaussGpuEvalConv_float(float ooSigma2,
                                   float* x_h, float* y_h, float* beta_h, float* gamma_h,
                                   int dimPoint, int dimVect, int nx, int ny);

extern "C" int GaussGpuEvalConv_double(double ooSigma2,
                                    double* x_h, double* y_h, double* beta_h, double* gamma_h,
                                    int dimPoint, int dimVect, int nx, int ny);

int GaussGpuEvalConv(float ooSigma2,
                  float* x_h, float* y_h, float* beta_h, float* gamma_h,
                  int dimPoint, int dimVect, int nx, int ny)
{
    return GaussGpuEvalConv_float(ooSigma2, x_h, y_h, beta_h, gamma_h, dimPoint, dimVect, nx, ny);
}

int GaussGpuEvalConv(double ooSigma2,
                  double* x_h, double* y_h, double* beta_h, double* gamma_h,
                  int dimPoint, int dimVect, int nx, int ny)
{
    return GaussGpuEvalConv_double(ooSigma2, x_h, y_h, beta_h, gamma_h, dimPoint, dimVect, nx, ny);
}

// for Grad1Conv

extern "C" int GaussGpuGrad1Conv_float(float ooSigma2,
                                       float* alpha_h, float* x_h, float* y_h, float* beta_h, float* gamma_h,
                                       int dimPoint, int dimVect, int nx, int ny);

extern "C" int GaussGpuGrad1Conv_double(double ooSigma2,
                                        double* alpha_h, double* x_h, double* y_h, double* beta_h, double* gamma_h,
                                        int dimPoint, int dimVect, int nx, int ny);

int GaussGpuGrad1Conv(float ooSigma2,
                      float* alpha_h, float* x_h, float* y_h, float* beta_h, float* gamma_h,
                      int dimPoint, int dimVect, int nx, int ny)
{
    return GaussGpuGrad1Conv_float(ooSigma2, alpha_h, x_h, y_h, beta_h, gamma_h, dimPoint, dimVect, nx, ny);
}

int GaussGpuGrad1Conv(double ooSigma2,
                      double* alpha_h, double* x_h, double* y_h, double* beta_h, double* gamma_h, 
                      int dimPoint, int dimVect, int nx, int ny)
{
    return GaussGpuGrad1Conv_double(ooSigma2, alpha_h, x_h, y_h, beta_h, gamma_h, dimPoint, dimVect, nx, ny);
}

// for GradConv

extern "C" int GaussGpuGradConv_float(float ooSigma2,
                                       float* alpha_h, float* x_h, float* beta_h, float* gamma_h,
                                       int dimPoint, int dimVect, int nx);

extern "C" int GaussGpuGradConv_double(double ooSigma2,
                                        double* alpha_h, double* x_h, double* beta_h, double* gamma_h,
                                        int dimPoint, int dimVect, int nx);

int GaussGpuGradConv(float ooSigma2,
                      float* alpha_h, float* x_h, float* beta_h, float* gamma_h,
                      int dimPoint, int dimVect, int nx)
{
    return GaussGpuGradConv_float(ooSigma2, alpha_h, x_h, beta_h, gamma_h, dimPoint, dimVect, nx);
}

int GaussGpuGradConv(double ooSigma2,
                      double* alpha_h, double* x_h, double* beta_h, double* gamma_h, 
                      int dimPoint, int dimVect, int nx)
{
    return GaussGpuGradConv_double(ooSigma2, alpha_h, x_h, beta_h, gamma_h, dimPoint, dimVect, nx);
}


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
            TYPE *xp, *yp, *betap, *gammap;
            xp = (TYPE*)x.data();
            yp = (TYPE*)y.data();
            betap = (TYPE*)beta.data();
            gammap = (TYPE*)gamma.data();
            GaussGpuEvalConv(ooSigma2,xp,yp,betap,gammap,DIMPOINT, DIMVECT, nx, ny);
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
            GaussGpuGrad1Conv(ooSigma2,alphap,xp,yp,betap,gammap,DIMPOINT, DIMVECT, nx, ny);
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
            GaussGpuGradConv(ooSigma2,alphap,xp,betap,gammap,DIMPOINT, DIMVECT, nx);
            return gamma;
        }
};

#endif // GAUSSGPUKERNEL
