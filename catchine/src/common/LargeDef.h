
// Large deformations class : implements large deformations parametrized via time-dependant momentum vectors

#ifndef LARGEDEF
#define LARGEDEF

#include "LargeDef_abstract.h"
#include "Kernel.h"
#include "Utils.h"

template < typename TYPE, int DIM >
class LargeDef : public LargeDef_abstract<TYPE,DIM>
{
	
private:
	
	using Evol<TYPE,DIM>::NumPoints;
	using Evol<TYPE,DIM>::Phi;
	
#ifdef DLIB_CATCHINE 
	typedef matrix<double,0,1> column_vector;
#endif
	typedef TinyVector<TYPE,DIM> Vect;
	typedef Array<Vect,1> ArrVect;
	typedef Array<ArrVect,1> ArrArrVect;
	typedef Kernel<Vect,Vect,Vect> LargeDefKernel;
	typedef Array<LargeDefKernel*,1> ArrKernelPtr;
	typedef Array<TYPE,1> OptimVect;

	int NumTimeSteps;
	double TimeInc;
	double RegWeight;
	
	bool testdeleteker, testsameker;
	
	ArrKernelPtr Ker;
	
	ArrArrVect Position;
	ArrArrVect Momentum;
	ArrArrVect Grad;
	
	void Integrate(ArrArrVect& A, ArrVect(LargeDef<TYPE,DIM>::*EquationA)(LargeDef<TYPE,DIM>::ArrArrVect&,int),
				   int tdeb, int tfin)
	{
		int Inc = sign(tfin-tdeb);
		ArrVect tempA(Range(1,A(1).rows()));
		for(int t=tdeb; t!=tfin; t+=Inc)
		{
			tempA = (this->*EquationA)(A,t);
			A(t+Inc) = A(t) + Inc * TimeInc * tempA;
			tempA += (this->*EquationA)(A,t+Inc);
			A(t+Inc) = A(t) + Inc * TimeInc * .5 * tempA;
		}
	}
	
	void Integrate(ArrArrVect& A, ArrVect(LargeDef<TYPE,DIM>::*EquationA)(LargeDef<TYPE,DIM>::ArrArrVect&,int),
				   ArrArrVect& B, ArrVect(LargeDef<TYPE,DIM>::*EquationB)(LargeDef<TYPE,DIM>::ArrArrVect&,int),
				   int tdeb, int tfin)
	{
		int Inc = sign(tfin-tdeb);
		ArrVect tempA(Range(1,A(1).rows())), tempB(Range(1,B(1).rows()));
		for(int t=tdeb; t!=tfin; t+=Inc)
		{
			tempA = (this->*EquationA)(A,t);
			tempB = (this->*EquationB)(B,t);
			A(t+Inc) = A(t) + Inc * TimeInc * tempA;
			B(t+Inc) = B(t) + Inc * TimeInc * tempB;
			A(t+Inc) = A(t) + Inc * TimeInc * .5 * ((this->*EquationA)(A,t+Inc)+tempA);
			B(t+Inc) = B(t) + Inc * TimeInc * .5 * ((this->*EquationB)(B,t+Inc)+tempB);
		}
	}
	
	ArrVect MomentumEq(ArrArrVect &Mom, int t)
	{
		ArrVect dMom(Range(1,Mom(t).rows()));
		dMom = - Ker(t)->Grad1Conv(Mom(t),Position(t),Position(t),Mom(t));
		return dMom;
	}
	
	ArrVect PositionEq(ArrArrVect &Pos, int t)
	{
		ArrVect dPos(Range(1,Pos(t).rows()));
		dPos = Ker(t)->EvalConv(Pos(t),Pos(t),Momentum(t));
		return dPos;
	}
	
	ArrVect FlowEq(ArrArrVect &Points, int t)
	{
		ArrVect dPoints(Range(1,Points(t).rows()));
		dPoints = Ker(t)->EvalConv(Points(t),Position(t),Momentum(t));
		return dPoints;
	}
	
	ArrVect GradEq(ArrArrVect &G, int t)
	{
		ArrVect beta(Range(1,NumPoints));
		beta = G(t) + RegWeight * Momentum(t);
		beta = -Ker(t)->GradConv(beta,Position(t), Momentum(t));
		return beta;
	}
	
	void ComputePosition()
	{
		Integrate(Position, &LargeDef<TYPE,DIM>::PositionEq,1,NumTimeSteps);
	}
	
	
	void Init(int ntimesteps, int npoints)
	{
		NumTimeSteps = ntimesteps;
		TimeInc = 1.0/(NumTimeSteps-1);
		NumPoints = npoints;
		AllocateArrArr(Position,NumTimeSteps,NumPoints);
		AllocateArrArr(Momentum,NumTimeSteps,NumPoints);
		AllocateArrArr(Grad,NumTimeSteps,NumPoints);
		Phi.reference(Position(NumTimeSteps));
	}
	
public:
	
	LargeDef(ifstream& f, int verbosemode=1)
	{
		if(verbosemode) cout << "Building LargeDef from ifstream" << endl;
		string Tag;
		for(int i=0; i<3; i++)
		{
			f >> Tag;
			if(!Tag.compare("Dimension="))
			{
				int dim;
				f >> dim;
				if(dim!=DIM)
				{
					cout << "Error constructing LargeDef from file; Dimension should be " << DIM << "." << endl;
					throw -1;
				}
			}
			else if(!Tag.compare("NumTimeSteps="))
				f >> NumTimeSteps;
			else if(!Tag.compare("NumPoints="))
				f >> NumPoints;
			else
			{
				cout << "Error constructing LargeDef from file; Dimension, NumTimeSteps, NumPoints should be first inputs." << endl;
				throw -1;
			}
		}
		Init(NumTimeSteps,NumPoints);
		Ker.resize(Range(1,NumTimeSteps));
		for(int i=0; i<2; i++)
		{
			f >> Tag;
			if(!Tag.compare("Kernel=same"))
			{
				testsameker = true;
				Ker(1) = ReadKernel<TYPE,DIM,DIM>(f);
				for(int t=2; t<NumTimeSteps+1; t++)
					Ker(t) = Ker(1);
			}
			else if(!Tag.compare("Kernel="))
			{
				testsameker = false;
				for(int t=1; t<NumTimeSteps+1; t++)
					Ker(t) = ReadKernel<TYPE,DIM,DIM>(f);
			}
			else if(!Tag.compare("Position(1)="))
			{
				ReadArr(Position(1),f);
				for(int t=1; t<NumTimeSteps+1; t++)
				{
					if(t!=1)
						Position(t) = Position(1);
					Momentum(t) = 0.0;
				}
			}
			else if(!Tag.compare("Position(1),Momentum(1)="))
			{
				ReadArr(Position(1),f);
				ReadArr(Momentum(1),f);
				Shooting();
			}
			else if(!Tag.compare("Position,Momentum="))
			{
				ReadArr(Position,f);
				ReadArr(Momentum,f);
			}
			else if(!Tag.compare("Position(1),Momentum="))
			{
				ReadArr(Position(1),f);
				ReadArr(Momentum,f);
				ComputePosition();
			}
			else
			{
				cout << "Error constructing LargeDef from file; Data and kernels should be given." << endl;
				throw -1;
			}
		}
		if(verbosemode) cout << "done." << endl;
		testdeleteker = true;
	}
	
	void Write(ofstream &f)
	{
		f << "LargeDef" << endl;
		f << "Dimension=" << endl << DIM << endl;
		f << "NumTimeSteps=" << endl << NumTimeSteps << endl;
		f << "NumPoints=" << endl << NumPoints << endl;
		f << "Kernel=" << endl;
		for(int t=1; t<NumTimeSteps+1; t++)
			Ker(t)->Write(f);
		f << "Position,Momentum=" << endl;
		WriteArr(Position,f);
		WriteArr(Momentum,f);
	}
	
	
	LargeDef(int ntimesteps, int npoints, ArrKernelPtr &ker)
	{
		Init(ntimesteps, npoints);
		Ker.reference(ker);
		testdeleteker = false;
	}
	
	LargeDef(int ntimesteps, ArrVect &x, ArrVect &mom, ArrKernelPtr &ker)
	{
		Init(ntimesteps, x.rows());
		Ker.reference(ker);
		Position(1) = x;
		Momentum(1) = mom;
		Shooting();
		testdeleteker = false;
	}
	
	LargeDef(int ntimesteps, ArrVect &x, ArrKernelPtr &ker)
	{
		Init(ntimesteps, x.rows());
		Ker.reference(ker);
		testdeleteker = false;
		for(int t=1; t<NumTimeSteps+1; t++)
		{
			Position(t) = x;
			Momentum(t) = 0.0;
		}
	}
	
	LargeDef(int ntimesteps, ArrArrVect &X, ArrArrVect &Mom, ArrKernelPtr &ker)
	{
		Init(ntimesteps, X(1).rows());
		testdeleteker = false;
		Ker.reference(ker);
		for(int t=1; t<NumTimeSteps+1; t++)
		{
			Position(t) = X(t);
			Momentum(t) = Mom(t);
		}
	}
	
	LargeDef(int ntimesteps, ifstream &f, ArrKernelPtr &ker)
	{
		ArrVect x;
		ReadArr(x,f);
		Init(ntimesteps, x.rows());
		Ker.reference(ker);
		testdeleteker = false;
		for(int t=1; t<NumTimeSteps+1; t++)
		{
			Position(t) = x;
			Momentum(t) = 0.0;
		}
	}
	
	LargeDef(ifstream &f, ArrKernelPtr &ker)
	{
		ReadArr(Position,f);
		ReadArr(Momentum,f);
		NumTimeSteps = Position.rows();
		TimeInc = 1.0/(NumTimeSteps-1);
		NumPoints = Position(1).rows();
		AllocateArrArr(Grad,NumTimeSteps,NumPoints);
		Ker.reference(ker);
		testdeleteker = false;
	}
	
	~LargeDef()
	{
		if(testdeleteker)
            {
            if(testsameker)
				delete(Ker(Ker.base(firstDim)));
			else
				for(int i=Ker.base(firstDim); i<Ker.base(firstDim)+Ker.rows(); i++)
					delete(Ker(i));
            }
	}
	
	
	int GetNumPoints()
	{
		return NumPoints;
	}
	
	int GetNumTimeSteps()
	{
		return NumTimeSteps;
	}
	
	void Shooting()
	{
		Integrate(Position, &LargeDef<TYPE,DIM>::PositionEq, Momentum, &LargeDef<TYPE,DIM>::MomentumEq,1,NumTimeSteps);
	}
	
	OptimVect CompGradient(ArrVect &EndGrad, double regweight)
	{
		RegWeight = regweight;
		Grad(NumTimeSteps) = EndGrad;
		Integrate(Grad, &LargeDef<TYPE,DIM>::GradEq,NumTimeSteps,1);
		ArrAddMult(Grad, 2*RegWeight, Momentum, Grad);
		OptimVect grad(Grad(1)(1).data(), shape(DIM*NumPoints*NumTimeSteps), neverDeleteData);
   		return grad;
	}
	
	void OneStepDescent(double StepSize)
	{
		ArrAddMult(Momentum, -StepSize, Grad, Momentum);
		ComputePosition();
	}
	
	double Energy()
	{
		static ArrVect temp(Range(1,NumPoints));
		double E = 0;
		for(int t=1; t<NumTimeSteps; t++)
		{
			temp = (Position(t+1)-Position(t))*(Momentum(t)+Momentum(t+1));
			E += fullsum(temp);
		}
		return E/2;
	}
	
	void Flow(ArrVect &X, ArrVect &EndPhiX)
	{
		Flow(X, EndPhiX, 1, NumTimeSteps);
	}
	
	void Flow(ArrArrVect &PhiX, int tdeb, int tfin)
	{
		Integrate(PhiX, &LargeDef<TYPE,DIM>::FlowEq,tdeb,tfin);
	}
	
	void Flow(ArrVect &X, ArrVect &EndPhiX, int tdeb, int tfin)
	{
		if(tdeb==-1)
			tdeb = 1;
		if(tfin==-1)
			tfin = NumTimeSteps;
		ArrArrVect PhiX;
		AllocateArrArr(PhiX,NumTimeSteps,X.rows());
		PhiX(tdeb) = X;
		Flow(PhiX,tdeb,tfin);
		EndPhiX = PhiX(tfin);
	}
	
	
	void AllFlow(ArrVect &X, ArrArrVect &PhiX)
	{
		PhiX(1) = X;
		Flow(PhiX,1,NumTimeSteps);
	}
	
	void SetStateVar(OptimVect statevar)
	{
		for(int i=0; i<NumVars(); i++)
			Momentum(1)(1).data()[i] = statevar.data()[i];
		ComputePosition();
	}
	
	void GetStateVar(OptimVect statevar)
	{
		for(int i=0; i<NumVars(); i++)
			statevar.data()[i] = Momentum(1)(1).data()[i];
	}
	
	int NumVars()
	{
		return DIM*NumPoints*NumTimeSteps;
	}
	
};

#endif // LARGEDEF

