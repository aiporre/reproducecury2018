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
	TYPE BreakRatio;
	
	void InitEvalGrad(void *po, Functional F, Gradient G)
	{
		Fct = F;
		Grad = G;
		PO = po;
	}
	
	static void DlibVect2Vect(Vect& v, const DlibVect& dv)
	{
		int n = dv.nr();
		v.resize(Range(1,n));
		for(int i=0; i<n; i++)
			v(i+1) = dv(i);
	}
	
	static void Vect2DlibVect(DlibVect& dv, const Vect& v)
	{
		int n = v.rows();
		dv.set_size(n);
		for(int i=0,j=v.base(0); i<n; i++, j++)
			dv(i) = v(j);
	}
	
	class WrapperEval
{
	Functional Fct;
	void *PO;
	Vect x;
	public :
	WrapperEval(void *po, Functional F) 
	{ 
		PO = po; 
		Fct = F; 
	}
	double operator() ( const DlibVect& statevar) const
	{
		DlibVect2Vect((Vect&)x,statevar);
		return Fct(PO,x);
	}
};

class WrapperGrad
{
	Gradient Grad;
	void *PO;
	Vect g, x;
	public :
	WrapperGrad(void *po, Gradient G) 
	{ 
		PO = po; 
		Grad = G;
	}
	const DlibVect operator() ( const DlibVect& statevar) const
	{
		Vect &gc = (Vect&)g;
		gc.resize(Range(1,statevar.nr()));
		DlibVect2Vect((Vect&)x,statevar);
		gc = Grad(PO,x);
		DlibVect gradstatevar;
		Vect2DlibVect(gradstatevar,gc);
		return gradstatevar;
	}
};

class Method
{
	public:
	virtual void Write(ofstream&) = 0;
	virtual void Run(void *po, Functional, Gradient, DlibVect&, TYPE) = 0;
};

class LBFGS : public Method
{
	int M, VerboseMode;
	public:
	LBFGS(ifstream& f, int verbosemode = 1)
	{
		VerboseMode = verbosemode;
		if(VerboseMode) cout << "Reading Dlib LBFGS parameters" << endl;
		string Tag;
		f >> Tag;
		if(!Tag.compare("M="))
			f >> M;
		else
		{
			cout << "error reading Dlib LBFGS parameters : LBFGS method should specify parameter M." << endl;
			throw -1;
		}
		if(VerboseMode) cout << "done" << endl;
	}
	void Write(ofstream& f)
	{
		f << "LBFGS" << endl;
		f << "M=" << endl << M << endl;
	}
	void Run(void *po, Functional F, Gradient G, DlibVect& statevar, TYPE BreakRatio)
	{
		objective_delta_stop_strategy Odss(BreakRatio);
		if(VerboseMode>=2) 
			Odss = Odss.be_verbose();
		find_min(lbfgs_search_strategy(M), Odss, WrapperEval(po,F), WrapperGrad(po,G), statevar, -1);
	}
};

class CG : public Method 
{
	int VerboseMode;
	public:
	CG(ifstream &f, int verbosemode = 1) { VerboseMode = verbosemode; };
	void Write(ofstream& f)
	{
		f << "CG" << endl;
	}
	void Run(void *po, Functional F, Gradient G, DlibVect& statevar, TYPE BreakRatio)
	{
		objective_delta_stop_strategy Odss(BreakRatio);
		if(VerboseMode>=2) 
			Odss = Odss.be_verbose();
		find_min(cg_search_strategy(), Odss, WrapperEval(po,F), WrapperGrad(po,G), statevar, -1);
	}
};

class BFGS : public Method 
{
	int VerboseMode;
	public:
	BFGS(ifstream &f, int verbosemode = 1) { VerboseMode = verbosemode; };
	void Write(ofstream& f)
	{
		f << "BFGS" << endl;
	}
	void Run(void *po, Functional F, Gradient G, DlibVect& statevar, TYPE BreakRatio)
	{
		//objective_delta_stop_strategy Odss(BreakRatio);
		gradient_norm_stop_strategy Odss(BreakRatio);
		if(VerboseMode>=2) 
			Odss = Odss.be_verbose();
		find_min(bfgs_search_strategy(), Odss, WrapperEval(po,F), WrapperGrad(po,G), statevar, -1);
	}
};

Method* ReadMethod(ifstream& f, int verbosemode)
{
	Method *Meth;
	string Tag;
	f >> Tag;
	if(!Tag.compare("LBFGS"))
		Meth = new LBFGS(f,verbosemode);
	else if(!Tag.compare("BFGS"))
		Meth = new BFGS(f,verbosemode);
	else if(!Tag.compare("CG"))
		Meth = new CG(f,verbosemode);
	return Meth;
}

	Method *Meth;

	public :
	Optim_Dlib(ifstream &f, int verbosemode = 1)
	{
		VerboseMode = verbosemode;
		if(VerboseMode) cout << "Reading Dlib Optimization parameters" << endl;
		
		string Tag;
		for(int i=0; i<2; i++)
		{
			f >> Tag;
			if(!Tag.compare("Method="))
				Meth = ReadMethod(f,VerboseMode);
			else if(!Tag.compare("BreakRatio="))
				f >> BreakRatio;
			else	
			{
				cout << "error reading Dlib Optimization parameters." << endl;
				throw -1;
			}
		}
		
		if(VerboseMode) cout << "done." << endl;
		
	}
    
	void Write(ofstream &f)
	{
		f << "Optim_Dlib" << endl;
		f << "Method=" << endl;
		Meth->Write(f);
		f << "BreakRatio=" << endl << BreakRatio << endl;
	}

	void Run(void *po, Functional F, Gradient G, Vect x)
	{
		DlibVect statevar;
		Vect2DlibVect(statevar,x);
		Meth->Run(po,F,G,statevar,BreakRatio);
	}
	
};

#endif // OPTIM_DLIB
