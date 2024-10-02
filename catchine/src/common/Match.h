#ifndef MATCH
#define MATCH

#include "Target.h"
#include "Evol.h"






template < typename TYPE, int DIM >
class Match
{
	typedef TinyVector<TYPE,DIM> Vect;
	typedef Array<Vect,1> ArrVect;
	typedef Array<Target<TYPE,DIM>*,1> ArrPntTarget;
	Evol<TYPE,DIM> *Ev;
	ArrPntTarget *Tg;
	Optim<TYPE> *Opt;
	typedef Array<TYPE,1> OptimVect;
	
	int NumTargets;
	double RegWeight;
	int VerboseMode, testdelete;
	
	ArrVect EndGrad;
	
	void Init(Evol<TYPE,DIM> *ev, ArrPntTarget *tg, Optim<TYPE>* opt, double regweight, int verbosemode)
	{
		Ev = ev;
		EndGrad.resize(Range(1,Ev->GetNumPoints()));
		Tg = tg;
		NumTargets = Tg->rows();
		RegWeight = regweight;
		VerboseMode = verbosemode;
		for(int n=1; n<NumTargets+1; n++)
			(*Tg)(n)->SetPhi(Ev);
		Opt = opt;
	}
	
public:
	
	Array<double,1> FunctValue;
	
	Match() {}
	
	Match(Evol<TYPE,DIM> *ev, ArrPntTarget *tg, double regweight, int verbosemode)
	{
		Init(ev, tg, regweight, verbosemode);
		testdelete = false;
	}
	
	Match(ifstream &f, double regweight, int verbosemode) { Init(f, regweight, verbosemode); }
	
	void Init(ifstream &f, double regweight, int verbosemode)
	{
		testdelete = true;
		string Tag;
		for(int i=0; i<3; i++)
		{
			f >> Tag;
			if(!Tag.compare("Evol"))
				Ev = ReadEvol<TYPE,__Dim__>(f, verbosemode);
			else if(!Tag.compare("Optim"))
				Opt = ReadOptim<TYPE>(f, verbosemode);
			else if(!Tag.compare("Targets"))
				Tg = ReadArrTargets<TYPE,__Dim__>(f, verbosemode);
			else
			{
				if(verbosemode) cout << "Error reading file; it should contain Deformation, Optim and Targets data. Exiting." << endl;
				throw -1;
			}
		}
		Init(Ev, Tg, Opt, regweight, verbosemode);
	}
	
	void Write(ofstream &fout)
	{	
	    fout << "Evol" << endl;
		Ev->Write(fout);
		fout << "Targets" << endl;
		WriteArrTargets(Tg,fout);
		fout << "Functional" << endl;
		WriteArr(FunctValue,fout);
	}
	
	~Match()
	{
		if(testdelete)
		{
			delete(Ev);
			for(int n=1; n<Tg->rows()+1; n++)
				delete((*Tg)(n));
			delete(Tg);
			delete(Opt);
		}
	}
	
	double Functional()
	{
		double J = 0;
		if(RegWeight)
			J = RegWeight * Ev->Energy();
		for(int n=1; n<NumTargets+1; n++)
			J += (*Tg)(n)->TargetWeight * (*Tg)(n)->Eval();
		return J;
	}
	
	OptimVect CompGradient()
	{
		EndGrad = 0.0;
		for(int n=1; n<NumTargets+1; n++)
			EndGrad((*Tg)(n)->RX) += (*Tg)(n)->TargetWeight * (*Tg)(n)->Gradient();
		return Ev->CompGradient(EndGrad,RegWeight);
	}
	
	static double WrapperOptimEval(void *po, OptimVect x)
	{
		Match* pa = (Match*) po;
		pa->Ev->SetStateVar(x);
		return pa->Functional();
	}	
	
	static OptimVect WrapperOptimGrad(void *po, OptimVect x)
	{
		Match* pa = (Match*) po;
		pa->Ev->SetStateVar(x);
		return pa->CompGradient();
	}	

	void DoOptimize()
	{
		OptimVect x(Range(1,Ev->NumVars()));
		Ev->GetStateVar(x);
		Opt->Run(this,WrapperOptimEval,WrapperOptimGrad,x);
	}
	
};





#endif // MATCH

