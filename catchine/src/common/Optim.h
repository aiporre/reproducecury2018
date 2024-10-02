#ifndef OPTIM
#define OPTIM

template < typename TYPE >
class Optim
{
	typedef Array<TYPE,1> Vect;
	typedef double (*Functional)(void*, Vect);
	typedef Vect (*Gradient)(void*, Vect);

	protected:
	int VerboseMode;	
	Array<double,1> FunctValue;

	public :
	virtual void Run(void*, Functional, Gradient, Vect) = 0;
    virtual ~Optim() { };
    
};

#endif // OPTIM







