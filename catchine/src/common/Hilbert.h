#ifndef HILBERT
#define HILBERT

template < typename VECTOR >
class Hilbert
{

    public:

        virtual double ScalProd(VECTOR&, VECTOR&) = 0;

        double SqNorm(VECTOR& X)
        {
            return ScalProd(X,X);
        }

        double Norm(VECTOR& X)
        {
            return sqrt(SqNorm(X));
        }

        double SqNormDiff(VECTOR& X, VECTOR& Y)
        {
            return ScalProd(X,X) + ScalProd(Y,Y) - 2*ScalProd(X,Y);
        }

        virtual ~Hilbert() {};
};

#endif // HILBERT

