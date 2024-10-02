#ifndef OPTIM_DLIB
#define OPTIM_DLIB

#include <dlib/optimization.h>
using namespace dlib;

#include "Optim.h"



// wrapper class for Dlib optimizaton


template < typename TYPE >
class Optim_Dlib : public Optim<TYPE>
{
	typedef Array<TYPE,1> Vect;
	typedef double (*Functional)(void*, Vect);
	typedef Vect (*Gradient)(void*, Vect);
	typedef matrix<double,0,1> DlibVect;
	
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
	
	Vect DlibVect2Vect(DlibVect& dv);
	{
		int n = dv.nc();
		Vect v(Range(1,n));
		for(int i=0; i<n; i++)
			v(i+1) = dv(i);
		return v;
	}
	
	DlibVect Vect2DlibVect(Vect& v);
	{
		int n = v.rows();
		DlibVect dv;
		dv.set_size(n);
		for(int i=0,j=x.base(0); i<n; i++, j++)
			dv(i) = v(j);
		return dv;
	}
	
	class WrapperEval
{
	Functional Fct;
	public :
	WrapperEval(Functional F) { Fct = F; }
	double operator() ( const DlibVect& statevar) const
	{
		Vect x = DlibVect2Vect(statevar);
		return Fct(x);
	}
};

class WrapperGrad
{
	Gradient Grad;
	public :
	WrapperGrad(Grad G) { Grad = G; }
	const DlibVect operator() ( const DlibVect& statevar) const
	{
		int n = statevar.nc();
		Vect x = DlibVect2Vect(statevar);
		Vect g(Range(1,n));
		g = Grad(x);
		DlibVect gradstatevar = Vect2DlibVect(g);
		return gradstatevar;
	}
};

	
	

	
	public :
	Optim_Dlib(ifstream &f, int verbosemode = 1)
	{
		VerboseMode = verbosemode;
		if(VerboseMode) cout << "Reading Optim_Dlib parameters" << endl;
/*		
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
*/	
	void Write(ofstream &f)
	{
		f << "Optim_Dlib" << endl;
/*
		f << "M=" << endl << M << endl;
		f << "LineSearch=" << endl << LineSearch << endl;
*/
	}
	
	
	void Run(void *po, Functional F, Gradient G, Vect x)
	{
		DlibVect statevar = Vect2DlibVect(x);
        find_min(lbfgs_search_strategy(10), objective_delta_stop_strategy(1e-7).be_verbose(), WrapperEval(F), WrapperGrad(G), statevar, -1);
	}

};

#endif // LBFGS_LIBLBFGS
