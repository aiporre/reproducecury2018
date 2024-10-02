#ifndef FREEEVOL
#define FREEEVOL

#include "Evol.h"

template < typename TYPE, int DIM >
class FreeEvol : public Evol<TYPE,DIM>
{
	
private:
	
	using Evol<TYPE,DIM>::NumPoints;
	using Evol<TYPE,DIM>::Phi;
	
#ifdef DLIB_CATCHINE 
	typedef matrix<double,0,1> column_vector;
#endif
	typedef TinyVector<TYPE,DIM> Vect;
	typedef Array<Vect,1> ArrVect;
	
	ArrVect Position;
	ArrVect Momentum;
	ArrVect Grad;
	
	void ComputePhi()
	{
		Phi = Position + Momentum;
	}
	
	void Init(int npoints)
	{
		NumPoints = npoints;
		Phi.resize(Range(1,NumPoints));
		Grad.resize(Range(1,NumPoints));
	}
	
public:
	
	FreeEvol(ifstream &f, int verbosemode = 1)
	{
		if(verbosemode) cout << "Building FreeEvol from ifstream" << endl;
		string Tag;
		for(int i=0; i<2; i++)
		{
			f >> Tag;
			if(!Tag.compare("Dimension="))
			{
				int dim;
				f >> dim;
				if(dim!=DIM)
				{
					cout << "Error constructing FreeEvol from file; Dimension should be " << DIM << "." << endl;
					throw -1;
				}
			}
			else if(!Tag.compare("NumPoints="))
				f >> NumPoints;
			else
			{
				cout << "Error constructing FreeEvol from file; Dimension and NumPoints should be first inputs." << endl;
				throw -1;
			}
		}
		Init(NumPoints);
		for(int i=0; i<1; i++)
		{
			f >> Tag;
			if(!Tag.compare("Position="))
			{
				ReadArr(Position,f);
				Momentum.resize(Range(1,NumPoints));
				Momentum = 0.0;
				Phi = Position;
			}
			else if(!Tag.compare("Position,Momentum="))
			{
				ReadArr(Position,f);
				ReadArr(Momentum,f);
				Phi.resize(Range(1,NumPoints));
				ComputePhi();
			}
			else if(!Tag.compare("Position,Momentum,Phi="))
			{
				ReadArr(Position,f);
				ReadArr(Momentum,f);
				ReadArr(Phi,f);
			}
			else
			{
				cout << "Error constructing FreeEvol from file; data should be given." << endl;
				throw -1;
			}
		}
		if(verbosemode) cout << "done." << endl;
	}
	
	void Write(ofstream &f)
	{
		f << "FreeEvol" << endl;
		f << "Dimension=" << endl << DIM << endl;
		f << "NumPoints=" << endl << NumPoints << endl;
		f << "Position,Momentum,Phi=" << endl;
		WriteArr(Position,f);
		WriteArr(Momentum,f);
		WriteArr(Phi,f);
	}
	
	FreeEvol(int npoints)
	{
		Init(npoints);
	}
	
	FreeEvol(ArrVect &x, ArrVect &mom)
	{
		Init(x.rows());
		Position.reference(x);
		Momentum.reference(mom);
		ComputePhi();
	}
	
	FreeEvol(ArrVect &x)
	{
		Init(x.rows());
		Position.reference(x);
		Momentum.resize(Range(1,NumPoints));
		Momentum = 0.0;
		Phi = Position;
	}
	
	int GetNumPoints()
	{
		return NumPoints;
	}
	
	void SetStateVar(const TYPE *statevar)
	{
		for(int i=0; i<DIM*NumPoints; i++)
			Momentum(1).data()[i] = statevar[i];
		ComputePhi();
	}
	
#ifdef DLIB_CATCHINE 
	void SetStateVar(const column_vector &statevar)
	{
		for(int i=0; i<DIM*NumPoints; i++)
			Momentum(1).data()[i] = statevar(i);
		ComputePhi();
	}
#endif
	
	
	TYPE* GetStateVarPointer()
	{
		return Momentum(1).data();
	}		
	
	int NumVars()
	{
		return DIM*NumPoints;
	}
	
	TYPE* CompGradient(ArrVect &EndGrad, double regweight)
	{
		(void) regweight;
		Grad = EndGrad;
		return Grad(1).data();
	}
	
	void OneStepDescent(double StepSize)
	{
		Momentum = Momentum - StepSize * Grad;
		ComputePhi();
	}
	
	double Energy()
	{
		return 0.0;
	}
	
	void Flow(ArrVect &X, ArrVect &PhiX)
	{
		PhiX =  X;
	}
	
	void Flow(ArrVect &X, ArrVect &PhiX, int dummy1, int dummy2)
	{
		// this function is for compatibility with Largedef.h
		(void) dummy1;
		(void) dummy2;
		PhiX =  X;
	}
};

#endif // FREEEVOL

