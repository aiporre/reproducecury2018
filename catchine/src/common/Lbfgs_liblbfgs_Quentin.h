#ifndef LBFGS_LIBLBFGS_QUENTIN
#define LBFGS_LIBLBFGS_QUENTIN

#include "Optim.h"
#include <setulb.cpp>


// wrapper class for liblbfgs 


template < typename TYPE >
class Lbfgs_liblbfgs_Quentin : public Optim<TYPE>
{
	typedef Array<TYPE,1> Vect;
	typedef double (*Functional)(void*, Vect);
	typedef Vect (*Gradient)(void*, Vect);
	
	using Optim<TYPE>::FunctValue;
	using Optim<TYPE>::VerboseMode;
	
	Functional Fct;
	Gradient Grad;
	void *PO;
	int M;
	
	void InitEvalGrad(void *po, Functional F, Gradient G)
	{
		Fct = F;
		Grad = G;
		PO = po;
	}
	
	TYPE WrapperEvalGrad(uvector &statevar, uvector &gradstatevar, const int n)
	{
		Vect g(Range(1,n));
		Vect x(Range(1,n));
		for(int i=0; i<n; i++)
			x(i+1) = statevar[i];
		g = Grad(PO,x);
		for(int i=0; i<n; i++)
			gradstatevar[i] = g(i+1);
		return Fct(PO,x);
	}		
		
	public :
	Lbfgs_liblbfgs_Quentin(ifstream &f, int verbosemode = 1)
	{
		VerboseMode = verbosemode;
		if(VerboseMode) cout << "Reading LBFGS (liblbfgs wrapper) parameters" << endl;
		
		string Tag;
		for(int i=0; i<1; i++)
		{
			f >> Tag;
			if(!Tag.compare("M="))
				f >> M;
			else	
			{
				cout << "error reading LBFGS (liblbfgs wrapper) parameters." << endl;
				throw -1;
			}
		}
		
		if(VerboseMode) cout << "done." << endl;
		
	}
	
	void Write(ofstream &f)
	{
		f << "Lbfgs_liblbfgs_Quentin" << endl;
		f << "M=" << endl << M << endl;
	}
	
	void Run(void *po, Functional F, Gradient G, Vect x)
	{
		InitEvalGrad(po,F,G);
		Run(x);
	}
	
	void Run(Vect x)
	{
		if(VerboseMode) cout << "Starting LBFGS (liblbfgs_Quentin wrapper)" << endl;
		int ret, n=x.rows();
		
		uvector statevar(n,0.0);
		uvector gradstatevar(n,0.0);
		TYPE *initsv = x.data();
		for(int i=0; i<n; i++)
			statevar[i] = initsv[i];
					
		TYPE fx;
		
		Lbfgs lb(n,M);    
		int k = 0;
		while(1)
		{
			int r = lb.iterate(statevar,fx,gradstatevar);
			if(r == LBFGS_FG)
				fx = WrapperEvalGrad(statevar,gradstatevar,n);
			else if(r == LBFGS_NEW_X)
			{
				k++;
				if(VerboseMode == 2)
				{
					printf("Iteration %d:\n", k);
					printf("  fx = %f\n", fx);
					printf("\n");
				}
			}
			else
				break;
	 	}
	}
};

#endif // LBFGS_LIBLBFGS_QUENTIN
