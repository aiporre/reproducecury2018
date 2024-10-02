#ifndef EVOL
#define EVOL

template < typename TYPE, int DIM >
class Evol
{
    protected:
        typedef TinyVector<TYPE,DIM> Vect;
        typedef Array<Vect,1> ArrVect;
	typedef Array<TYPE,1> OptimVect;
        int NumPoints;
    public:
        ArrVect Phi;
        virtual void Write(ofstream &f) = 0;
        double RegWeight;
        virtual int GetNumPoints() = 0;
        virtual double Energy() = 0;
	virtual void SetStateVar(OptimVect) = 0;
	virtual void GetStateVar(OptimVect) = 0;
	virtual int NumVars() = 0;
        virtual OptimVect CompGradient(ArrVect &, double) = 0;
        virtual void Flow(ArrVect &, ArrVect &) = 0;
        virtual void Flow(ArrVect &, ArrVect &, int, int) = 0;
        virtual ~Evol() { };
};

#endif // EVOL
