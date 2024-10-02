#ifndef CYCLEFUNCTIONS
#define CYCLEFUNCTIONS

#include "Function.h"

template < typename TYPE >
class ConeFunction : public Function< TYPE >
{
        Array<TYPE,1> Coeff, CoeffDiff;
        unsigned long int Nmax, N, Nd;
        TYPE Acc, Order, PI, PI2;
        std::complex<TYPE> I;
        TYPE funcoeff(TYPE n)
        {
            return 1.0/(pow2(2*n+1)*pow(pow2(2*n+1)+1,Order));
        };
        TYPE funcoeffdiff(TYPE n)
        {
            return (2*n+1)*funcoeff(n);
        };
        using Function<TYPE>::minArg;
        using Function<TYPE>::maxArg;

    public:

        ConeFunction(TYPE order=2.0, double acc=1e-6)
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
            f << "Cone,order=" << endl << Order << endl << "acc=" << Acc << endl;
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
            return 8*real(s)+PI2;
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
            return -8*imag(s);
        }

        TYPE Diff2(TYPE z)
        {
            assert(0);
            return z-z;
        }
};





template < typename TYPE >
class EdgeFunction : public Function< TYPE >
{
        Array<TYPE,1> Coeff, CoeffDiff;
        unsigned long int Nmax, N, Nd;
        TYPE Acc, Order, PI, Twop2mO;
        std::complex<TYPE> I;
        TYPE funcoeff(TYPE n)
        {
            return 1.0/(pow(pow2(2*n+1)+1,Order))+1.0/(pow(pow2(2*n-1)+1,Order));
        };
        TYPE funcoeffdiff(TYPE n)
        {
            return 2*n*funcoeff(n);
        };
        using Function<TYPE>::minArg;
        using Function<TYPE>::maxArg;

    public:

        EdgeFunction(TYPE order=2.0, double acc=1e-6)
        {
            PI = 3.14159265358979323846;
            minArg = 0.0;
            maxArg = PI;
            I = std::complex<TYPE>(0.0,1.0);
            Order = order;
            Twop2mO = pow(2.0-Order,2);
            Acc = acc;
            Nmax = 1000000;
            TYPE s=0, sd=0, sp=0, sdp=0;
            unsigned long int n;
            // get order of series truncation for required precision
            for(n=1; n<Nmax && funcoeff(n)>0; n++)
            {
                s += funcoeff(n);
                sd += funcoeffdiff(n);
            }
            for(n=1; s-sp>.1*Acc; n++)
                sp += funcoeff(n);
            N = n;
            N = 20;
            for(n=1; sd-sdp>.1*Acc; n++)
                sdp += funcoeffdiff(n);
            Nd = n;
            Nd = 20;
            // store expansion coeffs
            Coeff.resize(Range(1,N));
            for(n=1; n<N+1; n++)
                Coeff(n) = funcoeff(n);
            CoeffDiff.resize(Range(1,Nd));
            for(n=1; n<Nd+1; n++)
                CoeffDiff(n) = funcoeffdiff(n);
            cout << "done" << endl;
        }

        void Write(ofstream &f)
        {
            f << "Edge,order=" << endl << Order << endl << "acc=" << Acc << endl;
        }

        TYPE Eval(TYPE theta)
        {
            unsigned long int n;
            std::complex<TYPE> z, s, z2n, z2;
            z=exp(I*theta);
            s=0.0;
            z2 = z*z;
            z2n = z2;
            for(n=1; n<N+1; n++)
            {
                s += z2n*Coeff(n);
                z2n *= z2;
            }
            return 4*real(s)+Twop2mO;
        }

        TYPE Diff(TYPE theta)
        {
            unsigned long int n;
            std::complex<TYPE> z, s, z2n, z2;
            z=exp(I*theta);
            s=0.0;
            z2 = z*z;
            z2n = z2;
            for(n=1; n<Nd+1; n++)
            {
                s += z2n*CoeffDiff(n);
                z2n *= z2;
            }
            return -4*imag(s);
        }

        TYPE Diff2(TYPE z)
        {
            assert(0);
            return z-z;
        }
};






#endif  // CYCLEFUNCTIONS


