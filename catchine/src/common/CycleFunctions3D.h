#ifndef CYCLEFUNCTIONS3D
#define CYCLEFUNCTIONS3D

#include "Function.h"
#include "Utils.h"

template < typename TYPE >
Array<TYPE,1> NormAssocLegendre(Array<TYPE,1> &x, int m, int l)
{
    // computes fully normalized associated Legendre function N_m^l
    // evaluated at points in array x

    Range r = Range(x.base(0),x.rows()-x.base(0)-1);
    Array<TYPE,1> P1(r), P2(r), P(r);
    P1 = 0.0;
    P2 = pow((1.0-pow2(x)),(double)l/2);
    Array<int,1> p(2), a(2), q(1), b(1);
    p = 2*l-1;
    a = 2;
    q = 2*l;
    b = 1;
    P2 *= sqrt((double)(2*l+1)/2*quotfact(p,a,q,b));
    for(int k=l+1; k<m+1; k++)
    {
        P = (sqrt((double)(2*k+1)*(2*k-1)*(k-l)/(double)(k+l))*x*P2 -
             sqrt((double)(2*k+1)*(k+l-1)*(k-l)*(k-l-1)/(double)(2*k-3)/(k+l))*P1)/(double)(k-l);
        P1 = P2;
        P2 = P;
    }
    return P;
}

template < typename TYPE >
class ConeFunction3D : public Function< TYPE >
{
        Array<TYPE,1> Coeff, CoeffDiff;
        unsigned long int Nmax, N, Nd;
        TYPE Acc, Order, PI, PI2;
        std::complex<TYPE> I;

        TYPE Cml2(int m, int l)
        {
            int N = 100*m;
            Array<TYPE,1> x(N), P(N);
            firstIndex i;
            x = (i+.5)/N;
            P = NormAssocLegendre(x,m,l);
            return pow2(4.0*sum(P)/N/l)/PI;
        }

        TYPE funcoeff(int n)
        {
            int mmax = 50;
            int l = 2*n+1;
            TYPE c = 0.0;
            for(int m=l; m<mmax; m+=2)
                c += Cml2(m,l)/pow((TYPE)1.0+l*l*(l+1)*(l+1),Order);
            return c;
        };
        TYPE funcoeffdiff(int n)
        {
            return (2*n+1)*funcoeff(n);
        };
        using Function<TYPE>::minArg;
        using Function<TYPE>::maxArg;

    public:

        ConeFunction3D(TYPE order=2.0, double acc=1e-6)
        {
            PI = 3.14159265358979323846;
            PI2 = PI*PI;
            minArg = 0.0;
            maxArg = PI;
            I = std::complex<TYPE>(0.0,1.0);
            Order = order;
            Acc = acc;
            Nmax = 1000000;
            TYPE s=0, sd=0, sp=0, sdp=0;
            unsigned long int n;
            // get order of series truncation for required precision
            for(n=0; n<Nmax && funcoeff(n)>0; n++)
            {
                s += funcoeff(n);
                sd += funcoeffdiff(n);
            }
            for(n=0; s-sp>.1*Acc; n++)
                sp += funcoeff(n);
            N = n;
            N = 20;
            for(n=0; sd-sdp>.1*Acc; n++)
                sdp += funcoeffdiff(n);
            Nd = n;
            Nd = 20;
            // store expansion coeffs
            Coeff.resize(Range(0,N-1));
            for(n=0; n<N; n++)
                Coeff(n) = funcoeff(n);
            CoeffDiff.resize(Range(0,Nd-1));
            for(n=0; n<Nd; n++)
                CoeffDiff(n) = funcoeffdiff(n);
        }

        void Write(ofstream &f)
        {
            f << "Cone3D,order=" << endl << Order << endl << "acc=" << Acc << endl;
        }

        TYPE Eval(TYPE theta)
        {
            unsigned long int n;
            std::complex<TYPE> z, s, z2np1, z2;
            z=exp(I*theta);
            s=0.0;
            z2np1=z;
            z2 = z*z;
            for(n=0; n<N; n++)
            {
                s += z2np1*Coeff(n);
                z2np1 *= z2;
            }
            return real(s);
        }

        TYPE Diff(TYPE theta)
        {
            unsigned long int n;
            std::complex<TYPE> z, s, z2np1, z2;
            z = exp(I*theta);
            s = 0.0;
            z2np1 = z;
            z2 = z*z;
            for(n=0; n<Nd; n++)
            {
                s += z2np1*CoeffDiff(n);
                z2np1 *= z2;
            }
            return -imag(s);
        }

        TYPE Diff2(TYPE z)
        {
            assert(0);
            return z-z;
        }
};





template < typename TYPE >
class EdgeFunction3D : public Function< TYPE >
{
        Array<TYPE,1> Coeff, CoeffDiff;
        unsigned long int Nmax, N, Nd;
        TYPE Acc, Order, PI, PI2;
        std::complex<TYPE> I;

        TYPE Dml2(int m, int l)
        {
            int N = 100*m;
            Array<TYPE,1> x(N), P(N);
            firstIndex i;
            x = cos(.5*PI*(i+.5)/N);
            P = NormAssocLegendre(x,m,l);
//SHOW(m)
//SHOW(l)
//SHOW(pow2(sum(P)*2.0/N)*PI)
            return pow2(sum(P)*2.0/N)*PI;
        }

        TYPE funcoeff(int n)
        {
            int mmax = 50;
            int l = 2*n;
            TYPE c = 0.0;
            if(l)
            for(int m=l; m<mmax; m+=2)
                c += Dml2(m,l)/pow((TYPE)1.0+l*l*(l+1)*(l+1),Order);
//SHOW(n)
//SHOW(c)
            return c;
        };
        TYPE funcoeffdiff(int n)
        {
            return 2*n*funcoeff(n);
        };
        using Function<TYPE>::minArg;
        using Function<TYPE>::maxArg;

    public:

        EdgeFunction3D(TYPE order=2.0, double acc=1e-6)
        {
            PI = 3.14159265358979323846;
            PI2 = PI*PI;
            minArg = 0.0;
            maxArg = PI;
            I = std::complex<TYPE>(0.0,1.0);
            Order = order;
            Acc = acc;
            Nmax = 1000000;
            TYPE s=0, sd=0, sp=0, sdp=0;
            unsigned long int n;
            // get order of series truncation for required precision
            for(n=0; n<Nmax && funcoeff(n)>0; n++)
            {
                s += funcoeff(n);
                sd += funcoeffdiff(n);
            }
            for(n=0; s-sp>.1*Acc; n++)
                sp += funcoeff(n);
            N = n;
            N = 20;
            for(n=0; sd-sdp>.1*Acc; n++)
                sdp += funcoeffdiff(n);
            Nd = n;
            Nd = 20;
            // store expansion coeffs
            Coeff.resize(Range(0,N-1));
            for(n=0; n<N; n++)
                Coeff(n) = funcoeff(n);
            CoeffDiff.resize(Range(0,Nd-1));
            for(n=0; n<Nd; n++)
                CoeffDiff(n) = funcoeffdiff(n);
        }

        void Write(ofstream &f)
        {
            f << "Edge3D,order=" << endl << Order << endl << "acc=" << Acc << endl;
        }

        TYPE Eval(TYPE theta)
        {
            unsigned long int n;
            std::complex<TYPE> z, s, z2n, z2;
            z=exp(I*theta);
            s=0.0;
            z2n=1.0;
            z2 = z*z;
            for(n=0; n<N; n++)
            {
                s += z2n*Coeff(n);
                z2n *= z2;
            }
            return real(s);
        }

        TYPE Diff(TYPE theta)
        {
            unsigned long int n;
            std::complex<TYPE> z, s, z2n, z2;
            z = exp(I*theta);
            s = 0.0;
            z2n = z;
            z2 = z*z;
            for(n=0; n<Nd; n++)
            {
                s += z2n*CoeffDiff(n);
                z2n *= z2;
            }
            return -imag(s);
        }

        TYPE Diff2(TYPE z)
        {
            assert(0);
            return z-z;
        }
};






#endif  // CYCLEFUNCTIONS3D


