#ifndef FIXEDGRADDESCENT
#define FIXEDGRADDESCENT

#include "Optim.h"

// Gradient descent with fixed stepsize

template < typename TYPE >
class FixedGradDescent : public Optim<TYPE>
{
	typedef Array<TYPE,1> Vect;
	typedef double (*Functional)(void*, Vect);
	typedef Vect (*Gradient)(void*, Vect);

	using Optim<TYPE>::FunctValue;
	using Optim<TYPE>::VerboseMode;

	int Niters;
	TYPE StepSize;

	void Init()
	{
		FunctValue.resize(Range(0,Niters));
	}	

	public :
	FixedGradDescent(ifstream &f, int verbosemode = 1)
	{
		VerboseMode = verbosemode;
		if(VerboseMode) cout << "Reading Fixed Gradient Descent parameters" << endl;
		string Tag;
		for(int i=0; i<2; i++)
		{
	                f >> Tag;
	                if(!Tag.compare("Niters="))
				f >> Niters;
	                else if(!Tag.compare("StepSize="))
		                f >> StepSize;
        	        else	
			{
                		cout << "error reading Fixed Gradient Descent parameters." << endl;
				throw -1;
			}
            	}
		if(VerboseMode) cout << "done." << endl;
		Init();
        }

	void Write(ofstream &f)
	{
		f << "FixedGradDescent" << endl;
		f << "Niters=" << endl << Niters << endl;
		f << "StepSize" << endl << StepSize << endl;
	}

	void Run(void *po, Functional F, Gradient G, Vect x)
	{
		if(VerboseMode) cout << "Starting Gradient Descent (fixed stepsize)" << endl;
		Vect grad(Range(x.base(0),x.base(0)+x.extent(0)));
		FunctValue(0) = F(po,x);
		if(VerboseMode) cout << "starting with J=" << FunctValue(0) << endl;
		for(int i=0; i<Niters; i++)
		{
			x = x - StepSize * G(po,x);
			FunctValue(i+1) = F(po,x);
			if(VerboseMode>=2)
				cout << "iteration " << i+1 << ", J=" << FunctValue(i+1) << ", StepSize=" << StepSize << endl;
		}
		if(VerboseMode) cout<<"ending at iteration "<< Niters << ", J=" << FunctValue(Niters) << endl;
	}
};

#endif // FIXEDGRADDESCENT
