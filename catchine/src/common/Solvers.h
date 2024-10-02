
// ODE Solvers
// Implements Euler, Trapezoidal rule with Euler as predictor, and 4th-order Runge-Kutta

#ifndef SOLVERS
#define SOLVERS

#include "Utils.h"

template < class EQUATION, class BLITZARRAY >
class Solver
{
	public:
	virtual void Integrate(EQUATION, BLITZARRAY, double, double, int, Array<BLITZARRAY,1>& X = *(new Array<BLITZARRAY,1>(0))) = 0;
	virtual void Write(ofstream&) = 0;
};

template < class EQUATION, class BLITZARRAY >
class RungeKutta4 : public Solver<EQUATION,BLITZARRAY>
{
	public:
	void Write(ofstream& f)
	{
		f << "RungeKutta4" << endl;
	}
	void Integrate(EQUATION equation, BLITZARRAY x, double TimeStart, double TimeEnd, int NumTimeSteps, Array<BLITZARRAY,1>& X = *(new Array<BLITZARRAY,1>(0)))
	{
		Range indx = Range(x.base(0),x.base(0)+x.rows()-1);
		BLITZARRAY k1(indx), k2(indx), k3(indx), k4(indx), f(indx), dxdt(indx);
		double TimeInc = (TimeEnd-TimeStart)/(NumTimeSteps-1);
		double t;
		if(X.rows())
		{
			X.resize(Range(1,NumTimeSteps));
			X(1).resize(indx);
			X(1) = x;
		}
		for(int i=0, t=TimeStart; i<NumTimeSteps-1; i++, t+=TimeInc)
		{
			equation(x,dxdt,t);
			k1 = TimeInc * dxdt;
			f = x + .5*k1;
			equation(f,dxdt,t+.5*TimeInc);
			k2 = TimeInc * dxdt;
			f = x + .5*k2;
			equation(f,dxdt,t+.5*TimeInc);
			k3 = TimeInc * dxdt;
			f = x + .5*k3;
			equation(f,dxdt,t+TimeInc);
			k4 = TimeInc * dxdt;
			x += (1.0/6.0)*k1 + (1.0/3.0)*k2 + (1.0/3.0)*k3 + (1.0/6.0)*k4;
			if(X.rows())
			{
				X(i+2).resize(indx);
				X(i+2) = x;
			}
		}
	}
};

template < class EQUATION, class BLITZARRAY >
class EulerTrapezoidal : public Solver<EQUATION,BLITZARRAY>
{
	public:
	void Write(ofstream& f)
	{
		f << "EulerTrapezoidal" << endl;
	}
	void Integrate(EQUATION equation, BLITZARRAY x, double TimeStart, double TimeEnd, int NumTimeSteps, Array<BLITZARRAY,1>& X = *(new Array<BLITZARRAY,1>(0)))
	{
		Range indx = Range(x.base(0),x.base(0)+x.rows()-1);
		BLITZARRAY temp1(indx), temp2(indx), xtmp(indx);
		xtmp = x;
		double TimeInc = (TimeEnd-TimeStart)/(NumTimeSteps-1);
		double t;
		if(X.rows())
		{
			X.resize(Range(1,NumTimeSteps));
			X(1).resize(indx);
			X(1) = x;
		}
		for(int i=0, t=TimeStart; i<NumTimeSteps-1; i++, t+=TimeInc)
		{
			equation(x,temp1,t);
			xtmp += TimeInc * temp1;
			equation(xtmp,temp2,t+TimeInc);
			temp1 += temp2;
			x += (TimeInc*.5) * temp1;
			if(X.rows())
			{
				X(i+2).resize(indx);
				X(i+2) = x;
			}
		}
	}
};

template < class EQUATION, class BLITZARRAY >
class Euler : public Solver<EQUATION,BLITZARRAY>
{
	public:
	void Write(ofstream& f)
	{
		f << "Euler" << endl;
	}
	void Integrate(EQUATION equation, BLITZARRAY x, double TimeStart, double TimeEnd, int NumTimeSteps, Array<BLITZARRAY,1>& X = *(new Array<BLITZARRAY,1>(0)))
	{
		Range indx = Range(x.base(0),x.base(0)+x.rows()-1);
		BLITZARRAY temp(indx);
		double TimeInc = (TimeEnd-TimeStart)/(NumTimeSteps-1);
		double t;
		if(X.rows())
		{
			X.resize(Range(1,NumTimeSteps));
			X(1).resize(indx);
			X(1) = x;
		}
		for(int i=0, t=TimeStart; i<NumTimeSteps-1; i++, t+=TimeInc)
		{
			equation(x,temp,t);
			x += TimeInc * temp;
			if(X.rows())
			{
				X(i+2).resize(indx);
				X(i+2) = x;
			}
		}
	}
};

template < class EQUATION, class BLITZARRAY >
Solver<EQUATION,BLITZARRAY>* ReadSolver(ifstream& f)
{
	Solver<EQUATION,BLITZARRAY> *Solv = NULL;
	string Tag;
	f >> Tag;
	if(!Tag.compare("Euler"))
        	Solv = new Euler<EQUATION,BLITZARRAY>();
	else if(!Tag.compare("EulerTrapezoidal"))
		Solv = new EulerTrapezoidal<EQUATION,BLITZARRAY>();
	else if(!Tag.compare("RungeKutta4"))
		Solv = new RungeKutta4<EQUATION,BLITZARRAY>();
	else
	{
        	cout << "Error reading file: unknown Solver method '" << Tag << "'." << endl;
		throw -1;
    	}
    	return Solv;
}


#endif // SOLVERS
