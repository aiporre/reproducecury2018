
// Implementation of Translation and Rotation Invariant (TRI) kernels
// K(x,y) = ftilde(|x-y|)(x-y)(x-y)^T + fortho(|x-y|)Id


#ifndef TRIKERNEL
#define TRIKERNEL

#include "Kernel.h"
#include "Function.h"

template < typename TYPE, int DIMPOINT, int DIMVECT >
class TriKernel : public Kernel< TinyVector<TYPE,DIMPOINT> , TinyVector<TYPE,DIMVECT> , TinyVector<TYPE,DIMPOINT> >
{
	typedef TinyVector<TYPE,DIMPOINT> Point;
	typedef TinyVector<TYPE,DIMVECT> Vect;
	typedef TinyVector<TYPE,DIMPOINT> VectPoint;
	
protected:
	
	Function<TYPE> *FctTilde, *FctOrtho;
	bool testdelete;
	
public:

	TriKernel() { testdelete = false; }
	
	void Init(Function<TYPE> *fcttilde, Function<TYPE> *fctortho)
	{
		testdelete = false;
		FctTilde = fcttilde;
	    FctOrtho = fctortho;
	}	
	
	TriKernel(Function<TYPE> *fcttilde, Function<TYPE> *fctortho) { Init(fcttilde, fctortho); }
	
	void Init(ifstream &f)
	{
		testdelete = true;
		FctTilde = ReadFunction<TYPE>(f);
		string tag;
		f >> tag;
		if(tag!="functionOrtho=")
		{
			cout << "Error reading TriKernel from file; functionOrtho should be given." << endl;
			throw -1;
		}
		FctOrtho = ReadFunction<TYPE>(f);
	}
	
	TriKernel(ifstream &f) { Init(f); }
	
	~TriKernel()
	{
		if(testdelete)
		{
			delete(FctTilde);
			delete(FctOrtho);
		}
	}

	void Write(ofstream &f)
	{
		f << "Tri,functionTilde=" << endl;
		FctTilde->Write(f);
		f << "functionOrtho=" << endl;
		FctOrtho->Write(f);			
	}
	
	Vect Eval(Point &x, Point &y, Vect &beta)
	{
		Vect gamma;
		VectPoint xmy = x-y;
		TYPE r2 = sum(xmy*xmy);
		gamma = (FctTilde->Eval(r2) * sum(beta*xmy)) * xmy + FctOrtho->Eval(r2) * beta;
		return gamma;
	}
	
	Vect Diff1(Point &x, Point &y, Vect &beta, VectPoint &eta)
	{
		Vect gamma;
		VectPoint xmy = x-y;
		TYPE r2 = sum(xmy*xmy);
		TYPE xmyeta = sum(xmy*eta);
		TYPE xmybeta = sum(xmy*beta);
	    TYPE ftilde = FctTilde->Eval(r2);
		gamma = (FctTilde->Diff(r2) * 2.0 * xmyeta * xmybeta + ftilde * sum(eta*beta)) * xmy
		+ (ftilde * xmybeta) * eta + (FctOrtho->Diff(r2) * 2.0 * xmyeta) * beta;
	    return gamma;
	}
	
	Vect Diff2(Point &x, Point &y, Vect &beta, VectPoint &eta)
	{
	    return -Diff1(x,y,beta,eta);
	}
	
	VectPoint Grad1(Vect &alpha, Point &x, Point &y, Vect &beta)
	{
		VectPoint gamma;
		VectPoint xmy = x-y;
		TYPE r2 = sum(xmy*xmy);
		TYPE xmyalpha = sum(xmy*alpha);
		TYPE xmybeta = sum(xmy*beta);
		TYPE ftilde = FctTilde->Eval(r2);
		gamma = (FctTilde->Diff(r2) * 2.0 * xmyalpha * xmybeta + FctOrtho->Diff(r2) * 2.0 * sum(alpha*beta)) * xmy
		+ (ftilde * xmyalpha) * beta + (ftilde * xmybeta) * alpha;
		return gamma;
	}
	
	VectPoint Grad1Diff1(Vect &alpha, Point &x, Point &y, Vect &beta, VectPoint &eta)
	{
		VectPoint gamma;
		VectPoint xmy = x-y;
		TYPE r2 = sum(xmy*xmy);
		TYPE xmyalpha = sum(xmy*alpha);
		TYPE xmybeta = sum(xmy*beta);
		TYPE xmyeta = sum(xmy*eta);
		TYPE alphaeta = sum(alpha*eta);
		TYPE betaeta = sum(eta*beta);
		TYPE alphabeta = sum(alpha*beta);
		TYPE ftilde = FctTilde->Eval(r2);
		TYPE ftildep = FctTilde->Diff(r2);
		gamma = (4 * FctTilde->Diff2(r2) * xmyeta * xmybeta * xmyalpha 
				 + 2 * ftildep * (betaeta*xmyalpha+xmybeta*alphaeta)
				 + 4 * FctOrtho->Diff2(r2) * xmyeta * alphabeta) * xmy
		+ (2 * ftildep * xmybeta * xmyalpha + 2 * FctOrtho->Diff(r2) * alphabeta) * eta
		+ (2 * ftildep * xmyeta * xmyalpha + ftilde * alphaeta) * beta
		+ (2 * ftildep * xmyeta * xmybeta + ftilde * betaeta) * alpha;
		return gamma;
	}		
	
	VectPoint Grad2Diff1(Vect &alpha, Point &x, Point &y, Vect &beta, VectPoint &eta)
	{
		return -Grad1Diff1(alpha,x,y,beta,eta);
	}		
	
	
};

#endif // TRIKERNEL


