#ifndef SQDISTSCALARKERNEL
#define SQDISTSCALARKERNEL

#include "ScalarKernel.h"
#include "Function.h"

template < typename TYPE , int DIMPOINT, typename VECT >
class SqDistScalarKernel : public ScalarKernel<TYPE,DIMPOINT,VECT>
{
	private :
	
		//using ScalarKernel<TYPE,DIMPOINT,VECT>::Fct;
		//using ScalarKernel<TYPE,DIMPOINT,VECT>::testdelete;
	
		typedef TinyVector<TYPE,DIMPOINT> VectPoint;
		typedef Array<VectPoint,1> ArrVectPoint;
		typedef Array<VECT,1> ArrVect;
	
	public:
	
		SqDistScalarKernel(Function<TYPE> *fct = NULL) : ScalarKernel<TYPE,DIMPOINT,VECT>(fct) {}
	
		SqDistScalarKernel(ifstream &f) : ScalarKernel<TYPE,DIMPOINT,VECT>(f) {}
	
		void Write(ofstream &f)
		{
			f << "SqDistScalar,function=" << endl;
			this->Fct->Write(f);
		}
	
		TYPE EvalScal(VectPoint &z)
		{
			return sum(z*z);
		}
	
		TYPE DiffScal(VectPoint &z, VectPoint &beta)
		{
			return 2*sum(z*beta);
		}
	
		VectPoint GradScal(VectPoint& z)
		{
			return 2*z;
		}
	
		VectPoint GradDiffScal(VectPoint& z, VectPoint &beta)
		{
			return 2*beta;
		}

		ArrVectPoint GradDiffConv(ArrVectPoint &x, ArrVect &beta, ArrVectPoint &eta)
		{
			int nx = x.extent(firstDim);
			ArrVectPoint gamma(Range(1,nx));
			VectPoint ximxj, etaimetaj;
			TYPE r2;
			for(int i=1; i<nx+1; i++)
			{
				gamma(i) = 0.0;
				for(int j=1; j<nx+1; j++)
				{
					ximxj = x(i)-x(j);
					r2 = sum(ximxj*ximxj);
					etaimetaj = eta(i) - eta(j);
					gamma(i) += sum(beta(i)*beta(j)) * ( (2 * this->Fct->Diff2(r2) * sum(ximxj*etaimetaj)) * ximxj + this->Fct->Diff(r2) * etaimetaj );
				}
			}
			gamma *= 4;
			return gamma;
		}

};


#endif // SQDISTSCALARKERNEL
