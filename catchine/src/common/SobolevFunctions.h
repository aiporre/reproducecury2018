#ifndef SOBOLEVFUNCTIONS
#define SOBOLEVFUNCTIONS

#include "Function.h"

#include <boost/math/special_functions.hpp>
#include <boost/math/constants/constants.hpp>

#define PI boost::math::constants::pi<double>()

using namespace boost::math;

template < typename TYPE >
class H1Function : public Function<TYPE>
{
        TYPE Sigma, ooSigma2, oo2PISigma2, oo2PISigma4, oo2PISigma6;

    public:

        H1Function(TYPE sigma)
        {
            Sigma = sigma;
            ooSigma2 = 1.0/pow2(Sigma);
            oo2PISigma2 = ooSigma2/(2*PI);
            oo2PISigma4 = pow2(ooSigma2)/(2*PI);
            oo2PISigma6 = pow3(ooSigma2)/(2*PI);
        }
        TYPE Eval(TYPE r2)
        {
            if(r2) return oo2PISigma2*cyl_bessel_k<int,TYPE>(0,sqrt(r2*ooSigma2));
            else return 0.0;
        }
        TYPE Diff(TYPE r2)
        {
            if(r2) return -1.0/(8*PI*pow3(Sigma)*sqrt(r2))*cyl_bessel_k<int,TYPE>(1,sqrt(r2*ooSigma2));
            else return 0.0;
        }
        TYPE Diff2(TYPE r2)
        {
            return 1.0/(16*PI*pow3(Sigma)*r2)*(cyl_bessel_k<int,TYPE>(1,sqrt(r2*ooSigma2))/r2+1.0/(2*Sigma)*(cyl_bessel_k<int,TYPE>(2,sqrt(r2*ooSigma2))+cyl_bessel_k<int,TYPE>(0,sqrt(r2*ooSigma2))));
        }
        void Write(ofstream &f)
        {
            f << "H1,sigma=" << endl << Sigma << endl;
        }
};

template < typename TYPE >
class H2Function : public Function<TYPE>
{
        TYPE Sigma, ooSigma2, oo2PISigma2, oo2PISigma4, oo2PISigma6;

    public:

        H2Function(TYPE sigma)
        {
            Sigma = sigma;
            ooSigma2 = 1.0/pow2(Sigma);
            oo2PISigma2 = ooSigma2/(2*PI);
            oo2PISigma4 = pow2(ooSigma2)/(2*PI);
            oo2PISigma6 = pow3(ooSigma2)/(2*PI);
        }
        TYPE Eval(TYPE r2)
        {
            if(r2) return 1.0/(4*PI*pow3(Sigma))*sqrt(r2)*cyl_bessel_k<int,TYPE>(1,sqrt(r2*ooSigma2));
            else return 1.0/(4*PI*pow2(Sigma));
        }
        TYPE Diff(TYPE r2)
        {
            if(r2) return 1.0/(8*PI*pow3(Sigma))*(cyl_bessel_k<int,TYPE>(1,sqrt(r2*ooSigma2))/sqrt(r2)-1.0/(2*Sigma)*(cyl_bessel_k<int,TYPE>(0,sqrt(r2*ooSigma2))+cyl_bessel_k<int,TYPE>(2,sqrt(r2*ooSigma2))));
            else return 0.0;
        }
        TYPE Diff2(TYPE r2)
        {
            return 1.0/(8*PI*pow3(Sigma))*((-1.0/(4*Sigma*r2))*cyl_bessel_k<int,TYPE>(0,sqrt(r2*ooSigma2))+(-1.0/(4*Sigma*pow2(r2))+1.0/(4*pow2(Sigma)*sqrt(r2)))*cyl_bessel_k<int,TYPE>(1,sqrt(r2*ooSigma2))+(-1.0/(4*Sigma*r2))*cyl_bessel_k<int,TYPE>(2,sqrt(r2*ooSigma2))+(1.0/(8*pow2(Sigma)*sqrt(r2)))*cyl_bessel_k<int,TYPE>(3,sqrt(r2*ooSigma2)));
        }
        void Write(ofstream &f)
        {
            f << "H2,sigma=" << endl << Sigma << endl;
        }
};

#endif // SOBOLEVFUNCTIONS
