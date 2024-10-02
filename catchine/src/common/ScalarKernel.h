#ifndef SCALARKERNEL
#define SCALARKERNEL

#include "Kernel.h"
#include "Function.h"
#include "ReadFunction.h"

template < typename TYPE , int DIMPOINT, typename VECT >
class ScalarKernel : public Kernel< TinyVector<TYPE,DIMPOINT> , VECT , TinyVector<TYPE,DIMPOINT> >
{
	typedef TinyVector<TYPE,DIMPOINT> VectPoint;
	typedef Array<VECT,1> ArrVect;
	typedef Array<VectPoint,1> ArrVectPoint;
	
	virtual TYPE EvalScal(VectPoint&) = 0;
	virtual TYPE DiffScal(VectPoint&, VectPoint&) = 0;
	virtual VectPoint GradScal(VectPoint&) = 0;
	virtual VectPoint GradDiffScal(VectPoint&, VectPoint&) = 0;
	
    protected :
	
	Function<TYPE> *Fct;
	bool testdelete;
	
    public :
	
	ScalarKernel() { testdelete = false; }
	
	void Init(Function<TYPE> *fct)
	{ 
		testdelete = false;
		Fct = fct;
	}
	
	ScalarKernel(Function<TYPE> *fct) { Init(fct); }
	
	void Init(ifstream &f)
	{
		testdelete = true;
		Fct = ReadFunction<TYPE>(f);
	}
	
	ScalarKernel(ifstream &f) { Init(f); }
	
	~ScalarKernel()
	{
		if(testdelete)
			delete(Fct);
	}
	
	VECT Eval(VectPoint &x, VectPoint &y, VECT &beta)
	{
		VECT gamma;
		VectPoint xmy = x-y;
		gamma = Fct->Eval(EvalScal(xmy)) * beta;
		return gamma;
	}
	
	VECT Diff1(VectPoint &x, VectPoint &y, VECT &beta, VectPoint &eta)
	{
		VECT gamma;
		VectPoint xmy = x-y;
		gamma = (Fct->Diff(EvalScal(xmy)) * DiffScal(xmy,eta)) * beta;
		return gamma;
	}
	
	VECT Diff2(VectPoint &x, VectPoint &y, VECT &beta, VectPoint &eta)
	{
		return -Diff1(x,y,beta,eta);
	}
	
	VectPoint Grad1(VECT &alpha, VectPoint &x, VectPoint &y, VECT &beta)
	{
		VectPoint gamma;
		VectPoint xmy = x-y;
		gamma = sum(alpha*beta) * Fct->Diff(EvalScal(xmy)) * GradScal(xmy);
		return gamma;
	}
	
	VectPoint Grad1Diff1(VECT &alpha, VectPoint &x, VectPoint &y, VECT &beta, VectPoint &eta)
	{
		VectPoint gamma;
		VectPoint xmy = x-y;
	    gamma = ((Fct->Diff2(EvalScal(xmy)) * DiffScal(xmy,eta)) * GradScal(xmy) + (Fct->Diff(EvalScal(xmy)) * GradDiffScal(xmy,eta))) * sum(alpha*beta);
 	    return gamma;
	}
	
	VectPoint Grad2Diff1(VECT &alpha, VectPoint &x, VectPoint &y, VECT &beta, VectPoint &eta)
	{
 	    return -Grad1Diff1(alpha,x,y,beta,eta);
	}

	ArrVect DiffConv(ArrVectPoint &x, ArrVect &beta, ArrVectPoint &eta)
	{
		int nx = x.extent(firstDim);
	    ArrVect gamma(Range(1,nx));
		VectPoint ximxj, eimej;
		for(int i=1; i<nx+1; i++)
		{
			gamma(i) = 0.0;
			for(int j=1; j<nx+1; j++)
			{
				ximxj = x(i)-x(j);
				eimej = eta(i)-eta(j);
				gamma(i) += ( (Fct->Diff(EvalScal(ximxj)) * DiffScal(ximxj,eimej)) ) * beta(j);
			}
		}
		return gamma;
	}
		
};


#endif // SCALARKERNEL
