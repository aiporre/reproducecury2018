#ifndef SMALLDEF
#define SMALLDEF

#include "Evol.h"
#include "Kernel.h"

template < typename TYPE, int DIM >
class SmallDef : public Evol<TYPE,DIM>
{
	
private:
	
	using Evol<TYPE,DIM>::NumPoints;
	using Evol<TYPE,DIM>::Phi;
	
#ifdef DLIB_CATCHINE 
	typedef matrix<double,0,1> column_vector;
#endif
	typedef TinyVector<TYPE,DIM> Vect;
	typedef Array<Vect,1> ArrVect;
	typedef Kernel<Vect,Vect,Vect> SmallDefKernel;
	typedef Array<TYPE,1> OptimVect;
	
	ArrVect Position;
	ArrVect Momentum;
	ArrVect Grad;
	
	SmallDefKernel *Ker;
	
	void ComputePhi()
	{
		Phi = Position + Ker->EvalConv(Position,Position,Momentum);
	}
	
	void Init(int npoints)
	{
		NumPoints = npoints;
		Phi.resize(Range(1,NumPoints));
		Grad.resize(Range(1,NumPoints));
	}
	
public:
	
	SmallDef(ifstream &f, int verbosemode=1)
	{
		if(verbosemode) cout << "Building SmallDef from ifstream" << endl;
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
					cout << "Error constructing SmallDef from file; Dimension should be " << DIM << "." << endl;
					throw -1;
				}
			}
			else if(!Tag.compare("NumPoints="))
				f >> NumPoints;
			else
			{
				cout << "Error constructing SmallDef from file; Dimension and NumPoints should be first inputs." << endl;
				throw -1;
			}
		}
		Init(NumPoints);
		for(int i=0; i<2; i++)
		{
			f >> Tag;
			if(!Tag.compare("Kernel="))
				Ker = ReadKernel<TYPE,DIM,DIM>(f);
			else if(!Tag.compare("Position="))
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
				cout << "Error constructing SmallDef from file; data and kernels should be given." << endl;
				throw -1;
			}
		}
		if(verbosemode) cout << "done." << endl;
	}
	
	void Write(ofstream &f)
	{
		f << "SmallDef" << endl;
		f << "Dimension=" << endl << DIM << endl;
		f << "NumPoints=" << endl << NumPoints << endl;
		f << "Kernel=" << endl;
		Ker->Write(f);
		f << "Position,Momentum,Phi=" << endl;
		WriteArr(Position,f);
		WriteArr(Momentum,f);
		WriteArr(Phi,f);
	}
	
	
	SmallDef(int npoints, SmallDefKernel &ker)
	{
		Init(npoints);
		Ker = &ker;
	}
	
	SmallDef(ArrVect &x, ArrVect &mom, SmallDefKernel &ker)
	{
		Init(x.rows());
		Ker = &ker;
		Position.reference(x);
		Momentum.reference(mom);
		ComputePhi();
	}
	
	SmallDef(ArrVect &x, SmallDefKernel &ker)
	{
		Init(x.rows());
		Ker = &ker;
		Position.reference(x);
		Momentum.resize(Range(1,NumPoints));
		Momentum = 0.0;
		Phi = Position;
	}
	
	~SmallDef()
	{
		delete(Ker);
	}
	
	int GetNumPoints()
	{
		return NumPoints;
	}
		
	void SetStateVar(OptimVect statevar)
	{
		for(int i=0; i<DIM*NumPoints; i++)
			Momentum(1).data()[i] = statevar.data()[i];
		ComputePhi();
	}
	
	void GetStateVar(OptimVect statevar)
	{
		for(int i=0; i<NumVars(); i++)
			statevar.data()[i] = Momentum(1).data()[i];
	}		
	
	int NumVars()
	{
		return DIM*NumPoints;
	}
	
	OptimVect CompGradient(ArrVect &EndGrad, double regweight)
	{
		Grad = 2 * regweight * Momentum + EndGrad;
		OptimVect grad(Grad(1).data(), shape(DIM*NumPoints), neverDeleteData);
   		return grad;
	}
	
	void OneStepDescent(double StepSize)
	{
		Momentum = Momentum - StepSize * Grad;
		ComputePhi();
	}
	
	double Energy()
	{
		ArrVect temp(Range(1,NumPoints));
		temp = Momentum*(Phi-Position);
		return fullsum(temp);
	}
	
	void Flow(ArrVect &X, ArrVect &PhiX)
	{
		PhiX =  X + Ker->EvalConv(Position,X,Momentum);
	}
	
	void Flow(ArrVect &X, ArrVect &PhiX, int dummy1, int dummy2)
	{
		// this function is for compatibility with Largedef.h
		(void) dummy1;
		(void) dummy2;
		PhiX =  X + Ker->EvalConv(Position,X,Momentum);
	}
};

#endif // SMALLDEF

