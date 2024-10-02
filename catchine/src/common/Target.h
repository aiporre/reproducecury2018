#ifndef TARGET
#define TARGET

#include "Evol.h"

template < typename TYPE, int DIM >
class Target
{

        typedef TinyVector<TYPE,DIM> Vect;
        typedef Array<Vect,1> ArrVect;

	protected:

		ArrVect Phi;

    public:

        virtual void Write(ofstream &f) = 0;
        virtual double Eval() = 0;
        virtual ArrVect Gradient() = 0;
		void SetPhi(Evol<TYPE,DIM>* ev) { Phi.reference(ev->Phi(RX)); }
        Range RX;
        double TargetWeight;
        virtual ~Target() {};
};

template < typename TYPE, int DIM >
void WriteArrTargets(Array<Target<TYPE,DIM>*,1>* Tgt, ofstream &f)
{
    f << Tgt->rows() << endl;
    for(int i=1; i<Tgt->rows()+1; i++)
        (*Tgt)(i)->Write(f);
}


#endif // TARGET

