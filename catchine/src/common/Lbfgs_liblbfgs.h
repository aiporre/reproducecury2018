#ifndef LBFGS_LIBLBFGS
#define LBFGS_LIBLBFGS

#include "Optim.h"
#include <lbfgsfloat.c>
#include <lbfgsdouble.c>


// wrapper class for liblbfgs


template < typename TYPE >
class Lbfgs_liblbfgs : public Optim<TYPE>
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
	string LineSearch;
	
	
	
	void InitEvalGrad(void *po, Functional F, Gradient G)
	{
		Fct = F;
		Grad = G;
		PO = po;
	}
	
	static TYPE WrapperEvalGrad(void *po, const TYPE *statevar, TYPE *gradstatevar, const int n, const TYPE step)
	{
		Lbfgs_liblbfgs* pa = (Lbfgs_liblbfgs*) po;
		Vect g(Range(1,n));
		Vect x((TYPE*)statevar, shape(n), neverDeleteData);
		x.reindexSelf(1);
		if(gradstatevar!=NULL)
		{
			g = pa->Grad(pa->PO,x);
			for(int i=0; i<n; i++)
				gradstatevar[i] = g(i+1);
		}
		return pa->Fct(pa->PO,x);
	}		
	
	static int progress(
						void *instance,
						const TYPE *x,
						const TYPE *g,
						const TYPE fx,
						const TYPE xnorm,
						const TYPE gnorm,
						const TYPE step,
						int n,
						int k,
						int ls
						)
	{
		printf("Iteration %d:\n", k);
		//printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
		printf("  fx = %f\n", fx);
		printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
		printf("\n");
		return 0;
	}
	
	public :
	Lbfgs_liblbfgs(ifstream &f, int verbosemode = 1)
	{
		VerboseMode = verbosemode;
		if(VerboseMode) cout << "Reading LBFGS (liblbfgs wrapper) parameters" << endl;
		
		string Tag;
		for(int i=0; i<2; i++)
		{
			f >> Tag;
			if(!Tag.compare("M="))
				f >> M;
			else
				if(!Tag.compare("LineSearch="))
					f >> LineSearch;
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
		f << "Lbfgs_liblbfgs" << endl;
		f << "M=" << endl << M << endl;
		f << "LineSearch=" << endl << LineSearch << endl;
	}
	
	
	void Run(void *po, Functional F, Gradient G, Vect x)
	{
		InitEvalGrad(po,F,G);
		Run(x);
	}
	
	void Run(Vect x)
	{
		if(VerboseMode) cout << "Starting LBFGS (liblbfgs wrapper)" << endl;
		int ret, n=x.rows();
		TYPE *statevar = (TYPE*) malloc(x.rows()*sizeof(TYPE));
		TYPE *initsv = x.data();
		TYPE fx;
		lbfgs_parameter_t param;
		lbfgs_parameter_init(&param);
		if(LineSearch=="Backtracking")
			param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
		else
			if(LineSearch=="Backtracking_Armijo")
				param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
			else
				if(LineSearch=="Backtracking_Strong_Wolfe")
					param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
				else 
					if(LineSearch=="Backtracking_Wolfe")
						param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
					else 
						if(LineSearch=="MoreThuente")
							param.linesearch = LBFGS_LINESEARCH_MORETHUENTE;
						else 
							if(LineSearch=="Default")
								param.linesearch = LBFGS_LINESEARCH_DEFAULT;
							else
							{
								cout << "Invalid line search for LBFGS (liblbfgs wrapper)" << endl;
								throw -1;
							}
		param.m = M;
		for(int i=0; i<n; i++)
			statevar[i] = initsv[i];
		if(VerboseMode==2)
			ret = lbfgs(n, statevar, &fx, WrapperEvalGrad, progress, this, &param);
		else
			ret = lbfgs(n, statevar, &fx, WrapperEvalGrad, NULL, this, &param);
		if(VerboseMode) cout << "LBFGS terminated with return value " << ret << "." << endl;
		free(statevar);
		
	}
};

#endif // LBFGS_LIBLBFGS
