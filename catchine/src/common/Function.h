#ifndef FUNCTION
#define FUNCTION

// class for functions taking and returning scalars of same type (possibly complex)

template < typename TYPE >
class Function
{
    public:
        TYPE minArg, maxArg; // minimal and maximal argument values (used for real types only)
        virtual TYPE Eval(TYPE r2) = 0;
        virtual TYPE Diff(TYPE r2) = 0;
        virtual TYPE Diff2(TYPE r2) = 0;
	virtual void Write(ofstream &f) = 0;    
	virtual ~Function() {};
};

template < typename TYPE >
Function<TYPE>* ReadFunction(ifstream &f);

template < typename TYPE >
class SumFunction : public Function<TYPE>
{

	int NumFuns;
	Array<Function<TYPE>*,1> Funs;
	Array<TYPE,1> Weights;

    public:

    void Init(int numfuns)
    {
        NumFuns = numfuns;
        Funs.resize(Range(1,NumFuns));
        Weights.resize(Range(1,NumFuns));
    }
    
    SumFunction(ifstream &f)
    {
        int numfuns;
        f >> numfuns;
        Init(numfuns);
        for(int i=1; i<NumFuns+1; i++)
        {
            f >> Weights(i);
            Funs(i) = ReadFunction<TYPE>(f);
        }
    }
    
    SumFunction(Array<Function<TYPE>*,1>& funs, Array<TYPE,1>& weights)
    {
        Init(funs.rows());
        for(int i=1; i<NumFuns+1; i++)
        {
            Weights(i) = weights(i);
            Funs(i) = funs(i);
        }
    }
    
        TYPE Eval(TYPE r2)
		{
			TYPE res = 0.0;
			for(int i=1; i<NumFuns+1; i++)
				res += Weights(i) * Funs(i)->Eval(r2);
			return res;
		}

        TYPE Diff(TYPE r2)
		{
			TYPE res = 0.0;
			for(int i=1; i<NumFuns+1; i++)
				res += Weights(i) * Funs(i)->Diff(r2);
			return res;
		}

        TYPE Diff2(TYPE r2)
		{
			TYPE res = 0.0;
			for(int i=1; i<NumFuns+1; i++)
				res += Weights(i) * Funs(i)->Diff2(r2);
			return res;
		}

		void Write(ofstream &f)
		{
			f << "Sum" << endl << NumFuns << endl;
			for(int i=1; i<NumFuns+1; i++)
			{
				f << Weights(i) << endl;
				Funs(i)->Write(f);
			}
		}
    
		~SumFunction()
		{
			for(int i=1; i<NumFuns+1; i++)
				delete(Funs(i));
		};
};

#endif // FUNCTION

