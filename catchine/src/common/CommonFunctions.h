#ifndef COMMONFUNCTIONS
#define COMMONFUNCTIONS

#include "Function.h"

// Functions used in some common kernels.

// Gauss function is h(u) = exp(-u/sigma^2) (because input u is assumed to be a squared distance u=r^2)

template < typename TYPE >
class GaussFunction : public Function<TYPE>
{
        TYPE Sigma, ooSigma2, ooSigma4;
    public:
        GaussFunction(TYPE sigma)
        {
            Sigma = sigma;
            ooSigma2 = 1.0/(Sigma*Sigma);
            ooSigma4 = ooSigma2*ooSigma2;
        }
        TYPE Eval(TYPE r2)
        {
            return exp(-r2*ooSigma2);
        }
        TYPE Diff(TYPE r2)
        {
            return -ooSigma2*exp(-r2*ooSigma2);
        }
        TYPE Diff2(TYPE r2)
        {
            return ooSigma4*exp(-r2*ooSigma2);
        }
        void Write(ofstream &f)
        {
            f << "Gaussian,sigma=" << endl << Sigma << endl;
        }
};


// CauchyFunction is h(u) = 1/(1+u/sigma^2)

template < typename TYPE >
class CauchyFunction : public Function<TYPE>
{
        TYPE Sigma, ooSigma2, ooSigma4;
    public:
        CauchyFunction(TYPE sigma)
        {
            Sigma = sigma;
            ooSigma2 = 1.0/(Sigma*Sigma);
            ooSigma4 = ooSigma2*ooSigma2;
        }
        TYPE Eval(TYPE r2)
        {
            return 1.0/(1.0+r2*ooSigma2);
        }
        TYPE Diff(TYPE r2)
        {
		TYPE u = 1.0+r2*ooSigma2;
            return -ooSigma2*1.0/(u*u);
        }
        TYPE Diff2(TYPE r2)
        {
		TYPE u = 1.0+r2*ooSigma2;
            return ooSigma4*2.0/(u*u*u);
        }       
	void Write(ofstream &f)
        {
            f << "Cauchy,sigma=" << endl << Sigma << endl;
        }
};



// WeightedGaussFunction : h(u) = Coef*exp(-u/sigma^2)

template < typename TYPE >
class WeightedGaussFunction : public Function<TYPE>
{
        TYPE Coef, Sigma, ooSigma2, ooSigma4;
    public:
        WeightedGaussFunction(TYPE coef, TYPE sigma)
        {
	    Coef = coef;
            Sigma = sigma;
            ooSigma2 = 1.0/pow2(Sigma);
            ooSigma4 = pow2(ooSigma2);
        }
        TYPE Eval(TYPE r2)
        {
            return Coef*exp(-r2*ooSigma2);
        }
        TYPE Diff(TYPE r2)
        {
            return -Coef*ooSigma2*exp(-r2*ooSigma2);
        }
        TYPE Diff2(TYPE r2)
        {
            return Coef*ooSigma4*exp(-r2*ooSigma2);
        }
        void Write(ofstream &f)
        {
            f << "WeightedGaussian,coef=" << endl << Coef << endl << "sigma=" << endl << Sigma << endl;
        }
};


// SpecGaussFunction : h(u) = ( Coef0 + Coef1 * u ) * exp(-u/sigma^2)
// This is used for non-scalar TRI kernel types (See file TriKernel.h)

template < typename TYPE >
class SpecGaussFunction : public Function<TYPE>
{
        TYPE Coef0, Coef1, Sigma, ooSigma2, ooSigma4;
    public:
        SpecGaussFunction(TYPE coef0, TYPE coef1, TYPE sigma)
        {
			Coef0 = coef0;
			Coef1 = coef1;
            Sigma = sigma;
            ooSigma2 = 1.0/pow2(Sigma);
            ooSigma4 = pow2(ooSigma2);
        }
        TYPE Eval(TYPE r2)
        {
            return (Coef0+Coef1*r2)*exp(-r2*ooSigma2);
        }
        TYPE Diff(TYPE r2)
        {
            return (Coef1-(Coef0+Coef1*r2)*ooSigma2)*exp(-r2*ooSigma2);
        }
        TYPE Diff2(TYPE r2)
        {
            return (-2*Coef1*ooSigma2+ooSigma4*(Coef0+Coef1*r2))*exp(-r2*ooSigma2);
        }
        void Write(ofstream &f)
        {
            f << "SpecGaussian,coef0=" << endl << Coef0 << endl << "coef1=" << endl << Coef1 << endl << "sigma=" << endl << Sigma << endl;
        }
};


// Cauchy2Function is the quare of the Cauchy function : h(u) = 1/(1+u/sigma^2)^2

template < typename TYPE >
class Cauchy2Function : public Function<TYPE>
{
        TYPE Sigma, ooSigma2, ooSigma4;
    public:
        Cauchy2Function(TYPE sigma)
        {
            Sigma = sigma;
            ooSigma2 = 1.0/pow2(Sigma);
            ooSigma4 = pow2(ooSigma2);
        }
        TYPE Eval(TYPE r2)
        {
            return 1.0/pow2(1.0+r2*ooSigma2);
        }
        TYPE Diff(TYPE r2)
        {
            return -ooSigma2*2.0/pow3(1.0+r2*ooSigma2);
        }
        TYPE Diff2(TYPE r2)
        {
            return ooSigma4*6.0/pow4(1.0+r2*ooSigma2);
        }
        void Write(ofstream &f)
        {
            f << "Cauchy2,sigma=" << endl << Sigma << endl;
        }
};


// CubicFunction is h(u) = -u^1.5 = -r^3

template < typename TYPE >
class CubicFunction : public Function<TYPE>
{
        TYPE minushalf, oneandhalf;
    public:
        CubicFunction()
        {
            minushalf = -0.5f;
            oneandhalf = 1.5f;
        }
        TYPE Eval(TYPE r2)
        {
            return -pow(r2,oneandhalf);
        }
        TYPE Diff(TYPE r2)
        {
            return -1.5*sqrt(r2);
        }
        TYPE Diff2(TYPE r2)
        {
            return -0.75*pow(r2,minushalf);
        }
        void Write(ofstream &f)
        {
            f << "Cubic" << endl;
        }
};

#endif // COMMONFUNCTIONS



