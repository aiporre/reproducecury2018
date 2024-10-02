#ifndef LOOKUPFUNCTION1D
#define LOOKUPFUNCTION1D

#include "Function.h"

template < typename TYPE >
class LookUpFunction1D : public Function< TYPE >
{
        Array<TYPE,1> TabEval, TabDiff;
        TYPE Acc, maxArg, minArg, CoefEval, CoefDiff, Period;

        Function<TYPE>* Fct;

        void ComputeTable(TYPE(Function<TYPE>::*meth)(TYPE), Array<TYPE,1>& tabeval, TYPE maxder2, TYPE& coefx2t)
        {
            unsigned long int T = (unsigned long int) ceil(1.0/sqrt(.9*Acc*8.0/maxder2));
            tabeval.resize(Range(0,T));
            TYPE coeft2x = (maxArg-minArg)/T;
            for(unsigned long int t=0; t<T+1; t++)
                tabeval(t) = (Fct->*meth)(minArg+coeft2x*t);
            coefx2t = T/(maxArg-minArg);
        }

        inline TYPE LookUp(Array<TYPE,1>& tabeval, TYPE x, TYPE coefx2t)
        {
            // evaluate via linear interpolation on table
            unsigned long int t;
            TYPE a, tmp;
            tmp = (x-minArg)*coefx2t;
            t = (int)tmp;
            if(t==tabeval.rows()-1)
                t = tabeval.rows()-1;
            a = tmp-t;
            return (1-a)*tabeval(t)+a*tabeval(t+1);
        }

    public:

        LookUpFunction1D(Function<TYPE> *fct, double acc=1e-6) : Fct(fct), Acc(acc), minArg(fct->minArg), maxArg(fct->maxArg)
        {
            // compute look-up table for function
            TYPE maxDer2 = 600.54; // 6.54 > max(|f''|) (to be checked)
            ComputeTable(&Function<TYPE>::Eval,TabEval,maxDer2,CoefEval);
            // compute look-up table for diff
            TYPE maxDer3 = 1000.0;// 10 > max(|f'''|) maybe wrong, just a try
            ComputeTable(&Function<TYPE>::Diff,TabDiff,maxDer3,CoefDiff);
            Period = maxArg-minArg;
        }

        void Write(ofstream &f)
        {
            Fct->Write(f);
        }

        TYPE Eval(TYPE x)
        {
            return LookUp(TabEval,x,CoefEval);
        }

        TYPE Diff(TYPE x)
        {
            return LookUp(TabDiff,x,CoefDiff);
        }

        TYPE Diff2(TYPE z)
        {
            assert(0);
            return z-z;
        }

};

#endif  // LOOKUPFUNCTION1D

