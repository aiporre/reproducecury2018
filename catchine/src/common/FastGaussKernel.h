#ifndef FASTGAUSSKERNEL
#define FASTGAUSSKERNEL

#include "Utils.h"
#include "CommonFunctions.h"
#include "figtree.h"
#include "SqDistScalarKernel.h"

template < typename TYPE, int DIMPOINT, int DIMVECT >
class FastGaussKernel : public SqDistScalarKernel< TYPE, DIMPOINT, TinyVector<TYPE,DIMVECT> >
{

		//using SqDistScalarKernel< TYPE, DIMPOINT, TinyVector<TYPE,DIMVECT> >::Fct;
        //using SqDistScalarKernel< TYPE, DIMPOINT, TinyVector<TYPE,DIMVECT> >::testdelete;
		
        typedef TinyVector<TYPE,DIMPOINT> Point;
        typedef Array<Point,1> ArrPoint;
        typedef TinyVector<TYPE,DIMVECT> Vect;
        typedef Array<Vect,1> ArrVect;
        typedef TinyVector<TYPE,DIMPOINT> VectPoint;
        typedef Array<VectPoint,1> ArrVectPoint;

        TYPE Sigma, ooSigma2, Epsilon;
    
 public:

        void Write(ofstream &f)
        {
            f << "FastGauss,sigma,epsilon=" << endl << Sigma << " " << Epsilon << endl;
        }

        FastGaussKernel(TYPE sigma, TYPE epsilon)
        {
            if(sizeof(TYPE)==sizeof(float))
                cout << "Warning: float type will be converted to double type in Fast Gauss computations" << endl;
            Sigma = sigma;
            this->Fct = new GaussFunction<TYPE>(Sigma);
			this->testdelete = true;
            ooSigma2 = 1.0/pow2(Sigma);
            Epsilon = epsilon;

        }

        ArrVect EvalConv(ArrPoint &y, ArrPoint &x, ArrVect &alpha)
        {
            int nx = x.extent(firstDim);
            int ny = y.extent(firstDim);
            TYPE *ptRadius;
            TYPE dummy = 0;
            ptRadius = &dummy;
            Array<double,1> Coeffs(nx*DIMVECT);
            for(int d=0; d<DIMVECT; d++)
                for(int i=0; i<nx; i++)
                    Coeffs(d*nx+i) = alpha(i+1)(d);
            Array<double,1> gammacont(ny*DIMVECT);
            figtree(DIMPOINT,nx,ny,DIMVECT,contig(x,double()),(double)Sigma,Coeffs.data(),contig(y,double()),(double)Epsilon,gammacont.data(),FIGTREE_EVAL_AUTO,FIGTREE_PARAM_NON_UNIFORM,FIGTREE_TRUNC_CLUSTER,0);
            ArrVect gamma(Range(1,ny));
            for(int j=1; j<ny+1; j++)
                for(int d=0; d<DIMVECT; d++)
                    gamma(j)(d) = gammacont(d*ny+j-1);
            return gamma;
        }

        ArrVectPoint Grad1Conv(ArrVect &alpha, ArrVectPoint &x, ArrVectPoint &y, ArrVect &gamma)
        {
            int nx = x.extent(firstDim);
            int ny = y.extent(firstDim);
            ArrVectPoint beta(Range(1,nx));

            TYPE *ptRadius;
            TYPE dummy = 0;
            ptRadius = &dummy;

            int dimc = (DIMVECT*(1+DIMPOINT));
            Array<double,1> Coeffs(ny*dimc);

            for(int d=0; d<DIMVECT; d++)
            {
                for(int i=0; i<ny; i++)
                {
                    Coeffs(d*ny+i) = gamma(i+1)(d);
                    for(int dp=0; dp<DIMPOINT; dp++)
                        Coeffs((DIMVECT+DIMPOINT*d+dp)*ny+i) = gamma(i+1)(d)*y(i+1)(dp);
                }
            }

            Array<double,1> betacont(nx*dimc);
            figtree(DIMPOINT,ny,nx,dimc,contig(y,double()),(double)Sigma,Coeffs.data(),contig(x,double()),(double)Epsilon,betacont.data(),FIGTREE_EVAL_AUTO,FIGTREE_PARAM_NON_UNIFORM,FIGTREE_TRUNC_CLUSTER,0);

            int base1, base2;
            for(int i=1; i<nx+1; i++)
            {
                base1 = (i-1)*dimc;
                base2 = base1+DIMVECT;
                TYPE scalprod = 0;
                for(int dv=0; dv<DIMVECT; dv++)
                    scalprod += alpha(i)(dv)*betacont(dv*nx+i-1);
                for(int dp=0; dp<DIMPOINT; dp++)
                {
                    beta(i)(dp) = - scalprod * x(i)(dp);
                    for(int dv=0; dv<DIMVECT; dv++)
                        beta(i)(dp) += betacont((DIMVECT+dv*DIMPOINT+dp)*nx+i-1)*alpha(i)(dv);
                }
                beta(i) *= 2*ooSigma2;
            }
            return beta;
        }

        ArrVectPoint GradConv(ArrVect &alpha, ArrVectPoint &x, ArrVect &eta)
        {
            int nx = x.extent(firstDim);
            ArrVectPoint beta(Range(1,nx));

            TYPE *ptRadius;
            TYPE dummy = 0;
            ptRadius = &dummy;

            int dimc = 2*(DIMVECT*(1+DIMPOINT));
            Array<double,1> Coeffs(nx*dimc);


            for(int d=0; d<DIMVECT; d++)
            {
                for(int i=0; i<nx; i++)
                {
                    Coeffs(d*nx+i) = eta(i+1)(d);
                    for(int dp=0; dp<DIMPOINT; dp++)
                        Coeffs((DIMVECT+DIMPOINT*d+dp)*nx+i) = eta(i+1)(d)*x(i+1)(dp);
                    Coeffs((DIMVECT+DIMPOINT*DIMVECT+d)*nx+i) = alpha(i+1)(d);
                    for(int dp=0; dp<DIMPOINT; dp++)
                        Coeffs((2*DIMVECT+DIMPOINT*DIMVECT+DIMPOINT*d+dp)*nx+i) = alpha(i+1)(d)*x(i+1)(dp);
                }
            }

            Array<double,1> betacont(nx*dimc);
            figtree(DIMPOINT,nx,nx,dimc,contig(x,double()),(double)Sigma,Coeffs.data(),contig(x,double()),(double)Epsilon,betacont.data(),FIGTREE_EVAL_AUTO,FIGTREE_PARAM_NON_UNIFORM,FIGTREE_TRUNC_CLUSTER,0);

            for(int i=1; i<nx+1; i++)
            {
                int base = (i-1)*dimc;
                TYPE scalprod = 0;
                for(int dv=0; dv<DIMVECT; dv++)
                    scalprod += alpha(i)(dv)*betacont(dv*nx+i-1);
                base += DIMVECT;
                for(int dp=0; dp<DIMPOINT; dp++)
                {
                    beta(i)(dp) = 0;
                    for(int dv=0; dv<DIMVECT; dv++)
                        beta(i)(dp) += betacont((DIMVECT+dv*DIMPOINT+dp)*nx+i-1)*alpha(i)(dv);
                }
                base = (i-1)*dimc + (DIMPOINT+1)*DIMVECT;
                for(int dv=0; dv<DIMVECT; dv++)
                    scalprod += eta(i)(dv)*betacont((DIMVECT+DIMVECT*DIMPOINT+dv)*nx+i-1);
                base += DIMVECT;
                for(int dp=0; dp<DIMPOINT; dp++)
                {
                    beta(i)(dp) -= scalprod * x(i)(dp);
                    for(int dv=0; dv<DIMVECT; dv++)
                        beta(i)(dp) += betacont((2*DIMVECT+DIMVECT*DIMPOINT+dv*DIMPOINT+dp)*nx+i-1)*eta(i)(dv);
                }

                beta(i) *= 2*ooSigma2;
            }
            return beta;
        }
};

#endif // FASTGAUSSKERNEL
