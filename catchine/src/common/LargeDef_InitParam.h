#ifndef LARGEDEF_INITPARAM
#define LARGEDEF_INITPARAM

#include "LargeDef_abstract.h"
#include "Kernel.h"
#include "Utils.h"
#include "Solvers.h"

template < typename TYPE, int DIM >
class LargeDef_InitParam : public LargeDef_abstract<TYPE,DIM>
{
private:
	
	using Evol<TYPE,DIM>::NumPoints;
	using Evol<TYPE,DIM>::Phi;
	
	typedef TinyVector<TYPE,DIM> Vect;
	typedef Array<Vect,1> ArrVect;
	typedef Array<TYPE,1> OptimVect;
	typedef Kernel<Vect,Vect,Vect> LargeDefKernel;
	
	int NumTimeSteps;
	double TimeInc;
	double RegWeight;
	Range IndPos, IndMom, IndPosMom;
	LargeDefKernel* Ker;
	
	ArrVect FinalPosMom, FinalPos, FinalMom, InitPosMom, InitPos, InitMom, Grad;
	Array<ArrVect,1> AllPosMom;
	
	ArrVect GeodEq(const ArrVect &PosMom)
	{
		ArrVect Pos, Mom, dPosMom(IndPosMom);
		Pos.reference(PosMom(IndPos));
		Mom.reference(PosMom(IndMom));
		dPosMom(IndPos) = Ker->EvalConv(Pos,Pos,Mom);
		dPosMom(IndMom) = - Ker->Grad1Conv(Mom,Pos,Pos,Mom);
		return dPosMom;
	}

	ArrVect FlowEq(const ArrVect &PosMomPoints)
	{
		ArrVect Points, PosMom, Pos, Mom, dPosMomPoints(Range(1,PosMomPoints.rows()));
		Range IndPoints = Range(2*NumPoints+1,PosMomPoints.rows());
		PosMom.reference(PosMomPoints(IndPosMom));
		Pos.reference(PosMomPoints(IndPos));
		Mom.reference(PosMomPoints(IndMom));
		Points.reference(PosMomPoints(IndPoints));
		dPosMomPoints(IndPosMom) = GeodEq(PosMom);
		dPosMomPoints(IndPoints) = Ker->EvalConv(Points,Pos,Mom);
		return dPosMomPoints;
	}

	ArrVect GradEq(const ArrVect &PosMomBeta)
	{
		ArrVect BetaPos, BetaMom, PosMom, Pos, Mom, dPosMomBeta(Range(1,PosMomBeta.rows()));
		Range IndBetaPos = Range(2*NumPoints+1,3*NumPoints);
		Range IndBetaMom = Range(3*NumPoints+1,4*NumPoints);
		PosMom.reference(PosMomBeta(IndPosMom));
		Pos.reference(PosMomBeta(IndPos));
		Mom.reference(PosMomBeta(IndMom));
		BetaPos.reference(PosMomBeta(IndBetaPos));
		BetaMom.reference(PosMomBeta(IndBetaMom));
		dPosMomBeta(IndPosMom) = GeodEq(PosMom);
		dPosMomBeta(IndBetaMom) = Ker->DiffConv(Pos,Mom,BetaMom) - Ker->EvalConv(Pos,Pos,BetaPos) ;
		dPosMomBeta(IndBetaPos) = .5 * Ker->GradDiffConv(Pos,Mom,BetaMom) - Ker->GradConv(BetaPos,Pos,Mom);
		return dPosMomBeta;
	}

	class FunctorEq
	{
		LargeDef_InitParam* Parent;  
		typedef ArrVect (LargeDef_InitParam::*equation)(const ArrVect &);
		equation Eq;
	public:
		FunctorEq(LargeDef_InitParam* parent, equation eq) : Parent(parent), Eq(eq) { }
		void operator() (const ArrVect &X, ArrVect &dX, const double t)
		{
			dX = (Parent->*Eq)(X);
		}
	};

	Solver<FunctorEq,ArrVect> *Solv;
	
	void Init(int ntimesteps, int npoints)
	{
		NumTimeSteps = ntimesteps;
		TimeInc = 1.0/(NumTimeSteps-1);
		NumPoints = npoints;
		IndPosMom = Range(1,2*npoints);
		InitPosMom.resize(IndPosMom);
		FinalPosMom.resize(IndPosMom);
		IndPos = Range(1,npoints);
		IndMom = Range(npoints+1,2*npoints);
		InitPos.reference(InitPosMom(IndPos));
		InitMom.reference(InitPosMom(IndMom));
		FinalPos.reference(FinalPosMom(IndPos));
		FinalMom.reference(FinalPosMom(IndMom));
		Phi.reference(FinalPos);
		Grad.resize(Range(1,npoints));
	}
	
public:
	
	LargeDef_InitParam(ifstream& f, int verbosemode=1) //: WrapperGeodEq(this), WrapperFlowEq(this), WrapperGradEq(this)
	{
		if(verbosemode) cout << "Building LargeDef_InitParam from ifstream" << endl;
		string Tag;
		int npoints, ntimesteps;
		for(int i=0; i<4; i++)
		{
			f >> Tag;
			if(!Tag.compare("Dimension="))
			{
				int dim;
				f >> dim;
				if(dim!=DIM)
				{
					cout << "Error constructing LargeDef_InitParam from file; Dimension should be " << DIM << "." << endl;
					throw -1;
				}
			}
			else if(!Tag.compare("NumTimeSteps="))
				f >> ntimesteps;
			else if(!Tag.compare("NumPoints="))
				f >> npoints;
			else if(!Tag.compare("Solver="))
				Solv = ReadSolver<FunctorEq,ArrVect>(f);
			else
			{
				cout << "Error constructing LargeDef_InitParam from file; Dimension, NumPoints, NumTimeSteps, Solver should be first inputs." << endl;
				throw -1;
			}
		}
		Init(ntimesteps,npoints);
		for(int i=0; i<2; i++)
		{
			f >> Tag;			
			if(!Tag.compare("Kernel="))
				Ker = ReadKernel<TYPE,DIM,DIM>(f);
			else if(!Tag.compare("InitPos="))
			{
				ReadArr(InitPos,f);
				InitMom = 0.0;
			}
			else if(!Tag.compare("InitPos,InitMom="))
			{
				ReadArr(InitPos,f);
				ReadArr(InitMom,f);
			}
			else
			{
				cout << "Error constructing LargeDef_InitParam from file; Data and kernels should be given." << endl;
				throw -1;
			}
		}
		if(verbosemode) cout << "done." << endl;
	}
	
	void Write(ofstream &f)
	{
		f << "LargeDef_InitParam" << endl;
		f << "Dimension=" << endl << DIM << endl;
		f << "NumTimeSteps=" << endl << NumTimeSteps << endl;
		f << "NumPoints=" << endl << NumPoints << endl;
		f << "Solver=" << endl;
		Solv->Write(f);
		f << "Kernel=" << endl;
		Ker->Write(f);
		f << "InitPos,InitMom=" << endl;
		WriteArr(InitPos,f);
		WriteArr(InitMom,f);            
		f << "FinalPos=" << endl;
		Shooting();
		WriteArr(Phi,f);
		f << "distIdPhi=" << endl << sqrt(Energy()) << endl;
		AllShooting();
		if(AllPosMom.rows())
		{
			f << "AllPosMom=" << endl;
			WriteArr(AllPosMom,f);
		}
		else
			f << "AllPosMom=NotComputed" << endl;
	}
	
	int GetNumPoints()
	{
		return NumPoints;
	}
	
	void SetStateVar(OptimVect statevar)
	{
		for(int i=0; i<NumVars(); i++)
			InitMom(1).data()[i] = statevar.data()[i];
		Shooting();
	}
	
	void GetStateVar(OptimVect statevar)
	{
		for(int i=0; i<NumVars(); i++)
			statevar.data()[i] = InitMom(1).data()[i];
	}
	
	int NumVars()
	{
		return DIM*NumPoints;
	}
	
	void Shooting()
	{
		FinalPosMom = InitPosMom;
		Solv->Integrate(FunctorEq(this,&LargeDef_InitParam::GeodEq),FinalPosMom,0.0,1.0,NumTimeSteps);
	}
	
	void AllShooting()
	{
		FinalPosMom = InitPosMom;
		AllPosMom.resize(1); // just to have non-zero size as a flag for Integrate method; then will be resized inside
		Solv->Integrate(FunctorEq(this,&LargeDef_InitParam::GeodEq),FinalPosMom,0.0,1.0,NumTimeSteps,AllPosMom);
	}
	
	OptimVect CompGradient(ArrVect &EndGrad, double regweight)
	{
		Range IndBetaPos = Range(1+2*NumPoints,3*NumPoints);
		Range IndBetaMom = Range(1+3*NumPoints,4*NumPoints);
		RegWeight = regweight;
		Shooting();
		ArrVect PosMomBeta(Range(1,4*NumPoints));
		PosMomBeta(IndPosMom) = FinalPosMom;
		PosMomBeta(IndBetaMom) = 0.0;
		PosMomBeta(IndBetaPos) = EndGrad;
		Solv->Integrate(FunctorEq(this,&LargeDef_InitParam::GradEq),PosMomBeta,1.0,0.0,NumTimeSteps);
		Grad = PosMomBeta(IndBetaMom) + (2.0*RegWeight)*Ker->EvalConv(InitPos,InitPos,InitMom);
		OptimVect grad(Grad(1).data(), shape(DIM*NumPoints*NumTimeSteps), neverDeleteData);
		return grad;
	}
	
	double Energy()
	{
		ArrVect temp(Range(1,NumPoints));
		temp = InitMom * Ker->EvalConv(InitPos,InitPos,InitMom);
		return fullsum(temp);
	}
	
    void Flow(ArrVect &Points, ArrVect &PhiPoints)
    {
        Flow(Points, PhiPoints, 1, NumTimeSteps);
    }
              
	void Flow(ArrVect &Points, ArrVect &PhiPoints, int tdeb, int tfin)
	{
        double timestart = ((double)tdeb-1.0)/((double)NumTimeSteps-1);
        double timeend = ((double)tfin-1.0)/((double)NumTimeSteps-1);
		int a = 2*NumPoints;
		int b = Points.rows();
		ArrVect PosMomPoints(Range(1,a+b));
		PosMomPoints(IndPosMom) = InitPosMom;
		PosMomPoints(Range(a+1,a+b)) = Points;
		Solv->Integrate(FunctorEq(this,&LargeDef_InitParam::FlowEq),PosMomPoints,timestart,timeend,NumTimeSteps);
		PhiPoints = PosMomPoints(Range(a+1,a+b));
	}
	
	void AllFlow(ArrVect &Points, Array<ArrVect,1> &PhiPoints)
	{
		int a = 2*NumPoints;
		int b = Points.rows();
		ArrVect PosMomPoints(Range(1,a+b));
		PosMomPoints(IndPosMom) = InitPosMom;
		PosMomPoints(Range(a+1,a+b)) = Points;
		PhiPoints.resize(1); // just to have non-zero size as a flag for Integrate method; then will be resized inside
		Solv->Integrate(FunctorEq(this,&LargeDef_InitParam::FlowEq),PosMomPoints,0.0,1.0,NumTimeSteps,PhiPoints);
		ArrVect temp(Range(1,b));
		for(int i=PhiPoints.base(0); i<PhiPoints.base(0)+PhiPoints.rows(); i++)
		{
			temp = PhiPoints(i)(Range(a+1,a+b));
			PhiPoints(i).resize(Range(1,b));
			PhiPoints(i) = temp;
		}
	}
	
};

#endif // LARGEDEF_INITPARAM

