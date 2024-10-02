#ifndef LARGEDEF_ABSTRACT
#define LARGEDEF_ABSTRACT

#include "Evol.h"

template < typename TYPE, int DIM >
class LargeDef_abstract : public Evol<TYPE,DIM>
{
	typedef TinyVector<TYPE,DIM> Vect;
	typedef Array<Vect,1> ArrVect;
public:
	virtual void AllFlow(ArrVect &, Array<ArrVect,1> &) = 0;
};

#endif // LARGEDEF_ABSTRACT
