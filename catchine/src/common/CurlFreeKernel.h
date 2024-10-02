#ifndef CURLFREEKERNEL
#define CURLFREEKERNEL

// Implementation of Curl Free kernels
// This file is old and should be checked and/or

#include "Kernel.h"
#include "Function.h"

template < typename TYPE, int DIMPOINT, int DIMVECT >
class CurlFreeKernel : public Kernel< TinyVector<TYPE,DIMPOINT> , TinyVector<TYPE,DIMVECT> , TinyVector<TYPE,DIMPOINT> >
{

        typedef TinyVector<TYPE,DIMPOINT> Point;
        typedef TinyVector<TYPE,DIMVECT> Vect;
        typedef TinyVector<TYPE,DIMPOINT> VectPoint;

    protected:

        Function<TYPE> *Fct;

    public:

        CurlFreeKernel(Function<TYPE> *fct)
        {
            Fct = fct;
        }

        void Write(ofstream &f)
        {
            f << "CurlFree,function=" << endl;
            Fct->Write(f);
        }

        Vect Eval(Point &x, Point &y, Vect &alpha)
        {
            Vect gamma;
            VectPoint ymx = y-x;
            TYPE r2 = sum(ymx*ymx);
            gamma = Fct->Eval(r2) * alpha + 2 * Fct->Diff(r2) * sum(ymx*alpha) * ymx;
            return gamma;
        }

        Vect Diff(Point &x, VectPoint &beta, Point &y, Vect &alpha)
        {
            Vect gamma;
            VectPoint ymx = y-x;
            TYPE r2 = sum(ymx*ymx);
            TYPE ymxbeta = sum(ymx*beta);
            TYPE ymxalpha = sum(ymx*alpha);
            TYPE alphabeta = sum(beta*alpha);
            gamma = - 2 * Fct->Diff(r2) *
                    ( ymxbeta * alpha + alphabeta * ymx + ymxalpha * beta )
                    - 4 * Fct->Diff2(r2) * ymxalpha * ymxbeta * ymx;
            return gamma;
        }

        VectPoint Grad(Vect &gamma, Point &x, Point &y, Vect &alpha)
        {
            VectPoint beta;
            VectPoint ymx = y-x;
            TYPE r2 = sum(ymx*ymx);
            TYPE ymxgamma = sum(ymx*gamma);
            TYPE ymxalpha = sum(ymx*alpha);
            TYPE alphagamma = sum(gamma*alpha);
            beta = - 2 * Fct->Diff(r2) *
                   ( alphagamma * ymx + ymxgamma * alpha + ymxalpha * gamma )
                   - 4 * Fct->Diff2(r2) * ymxalpha * ymxgamma * ymx;
            return beta;
        }
};

#endif // CURLFREEKERNEL

