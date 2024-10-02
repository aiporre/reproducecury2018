#ifndef ADAPTGRADDESCENT
#define ADAPTGRADDESCENT

#include "Optim.h"

// Gradient descent with adaptive stepsize

template < typename TYPE >
class AdaptGradDescent : public Optim<TYPE>
{
	typedef Array<TYPE,1> Vect;
	typedef double (*Functional)(void*, Vect);
	typedef Vect (*Gradient)(void*, Vect);

	using Optim<TYPE>::FunctValue;
	using Optim<TYPE>::VerboseMode;

	int Niters, LoopBreak;
	TYPE StepSize;
	double BreakRatio, Malus, Bonus;

	void Init()
	{
		FunctValue.resize(Range(0,Niters));
		Malus = 0.5;
		Bonus = 1.2;
	}	

	public :
	AdaptGradDescent(ifstream &f, int verbosemode = 1)
	{
		VerboseMode = verbosemode;
		if(VerboseMode) cout << "Reading Adaptive Gradient Descent parameters" << endl;
		string Tag;
		for(int i=0; i<4; i++)
		{
	                f >> Tag;
	                if(!Tag.compare("Niters="))
				f >> Niters;
	                else if(!Tag.compare("StepSize="))
		                f >> StepSize;
	                else if(!Tag.compare("LoopBreak="))
		                f >> LoopBreak;
	                else if(!Tag.compare("BreakRatio="))
		                f >> BreakRatio;
        	        else	
			{
                		cout << "error reading Adaptive Gradient Descent parameters." << endl;
				throw -1;
			}
            	}
		if(VerboseMode) cout << "done." << endl;
		Init();
        }

	void Write(ofstream &f)
	{
		f << "AdaptGradDescent" << endl;
		f << "Niters=" << endl << Niters << endl;
		f << "StepSize" << endl << StepSize << endl;
		f << "LoopBreak" << endl << LoopBreak << endl;
		f << "BreakRatio" << endl << BreakRatio << endl;
	}

	void Run(void *po, Functional F, Gradient G, Vect x)
	{
		if(VerboseMode) cout << "Starting Gradient Descent (adaptive stepsize)" << endl;
	double MalusVar;
	double Jp, J0;
	int i;
		Vect grad(Range(x.base(0),x.base(0)+x.extent(0)));
	Vect xcurr(Range(x.base(0),x.base(0)+x.extent(0)));
	FunctValue(0) = F(po,x);
	J0 = FunctValue(0);
	if(VerboseMode) cout << "starting with J=" << FunctValue(0) << endl;
	for(i=0; i<Niters; i++)
	{
		xcurr = x;
		Jp = FunctValue(i);
		grad = G(po,x);
		MalusVar = -Malus;
		int loopiter;
		for(loopiter=0; loopiter<LoopBreak; loopiter++)
		{
			x -= StepSize * grad;
			FunctValue(i+1) = F(po,x);
			if((Jp-FunctValue(i+1))>(J0-Jp)*BreakRatio)
				break;
			StepSize *= MalusVar;
			MalusVar = Malus;
		}
		if(loopiter==LoopBreak)
		{
			if(VerboseMode) cout<<"no significant minimization, stopping"<<endl;
			x = xcurr;
			F(po,x); // to reset the state
			FunctValue.resizeAndPreserve(i+1);
			FunctValue.reindexSelf(0);
			break;
		}
		StepSize = abs(StepSize) * Bonus;
		if(VerboseMode>=2)
			cout << "iteration " << i+1 << ", J=" << FunctValue(i+1) << ", StepSize=" << StepSize << endl;
	}

	}
};

#endif // ADAPTGRADDESCENT
