#ifndef ANGLEKERNEL
#define ANGLEKERNEL

// class for kernels on S^n which are expressed as functions of the angle between points

#include "ScalarKernel.h"

template < typename TYPE, int DIM, typename VECT >
class AngleKernel : public ScalarKernel<TYPE,DIM,VECT>
{
        typedef TinyVector<TYPE,DIM> VectPoint;
        using ScalarKernel<TYPE,DIM,VECT>::Fct;

    public:

        AngleKernel(Function<TYPE> *fct = NULL) : ScalarKernel<TYPE,DIM,VECT>(fct) { }

        void Write(ofstream &f)
        {
            f << "Angle,function=" << endl;
            Fct->Write(f);
        }

        TYPE EvalScal(VectPoint &z)
        {
          return acos(max(-1.0,min(1.0,1.0-.5*sum(z*z))));
        }

        TYPE DiffScal(VectPoint &z, VectPoint &beta)
        {
            TYPE s=sum(z*z);
            if(abs(s)>0.0 && abs(s)<0.99999)
            	return -sum(beta*z)/sqrt(s-.25*pow2(s));
            else
            	return 0.0;
        }

        VectPoint GradScal(VectPoint &z)
        {
            VectPoint G;
            TYPE s=sum(z*z);
            if(abs(s)>0.0 && abs(s)<0.99999)
            	G=(-1.0/sqrt(s-.25*pow2(s)))*z;
            else
            	G = 0.0;
            return G;
        }

	VectPoint GradDiffScal(VectPoint &z, VectPoint &beta)
	{
		VectPoint G;
		TYPE s = sum(z*z);
		if(abs(s)>0.0 && abs(s)<0.99999)
		{
			TYPE coef = s-.25*pow2(s);
			G = (1.0/sqrt(coef)) * ( beta - (sum(z*beta)*(1.0-2.0*s)/coef) * z );
		}	 
		else
			G = 0.0;	
		return G;
	}

};

#endif // ANGLEKERNEL
