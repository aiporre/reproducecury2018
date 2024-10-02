#ifndef READFUNCTION
#define READFUNCTION

#include "Function.h"

template < typename TYPE >
Function<TYPE>* ReadFunction(ifstream &f)
{
    Function<TYPE> *fct = NULL;
    string Tag;
    f >> Tag;
    double sigma;
    if(Tag.compare("Cubic") && Tag.compare("Sum"))
        f >> sigma;
#ifdef COMMONFUNCTIONS
    if(!Tag.compare("Gaussian,sigma="))
        fct = new GaussFunction<TYPE>(sigma);
    else if(!Tag.compare("WeightedGaussian,coef="))
	{
		double coef = sigma;
		f >> Tag;
		f >> sigma;
        fct = new WeightedGaussFunction<TYPE>(coef,sigma);
	}
    else if(!Tag.compare("SpecGaussian,coef0="))
	{
		double coef1, coef0 = sigma;
		f >> Tag;
		f >> coef1;
		f >> Tag;
		f >> sigma;
        fct = new SpecGaussFunction<TYPE>(coef0,coef1,sigma);
	}
    else if(!Tag.compare("Cauchy,sigma="))
        fct = new CauchyFunction<TYPE>(sigma);
    else if(!Tag.compare("Cauchy2,sigma="))
        fct = new Cauchy2Function<TYPE>(sigma);
    else if(!Tag.compare("Cubic"))
        fct = new CubicFunction<TYPE>();
    else
#endif // COMMONFUNCTIONS
#ifdef SOBOLEVFUNCTIONS
        if(!Tag.compare("H1,sigma="))
            fct = new H1Function<TYPE>(sigma);
        else if(!Tag.compare("H2,sigma="))
            fct = new H2Function<TYPE>(sigma);
        else
#endif // SOBOLEVFUNCTIONS
        if(!Tag.compare("Sum"))
            fct = new SumFunction<TYPE>(f);
        else
	{
            cout << "Error constructing function from file: unknown function type '" << Tag << "'." << endl;
	    throw -1;
	}
    return fct;
}

#endif // READFUNCTION

