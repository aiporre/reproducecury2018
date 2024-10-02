#ifndef VECTDIRACSMEASURE
#define VECTDIRACSMEASURE

template < typename TYPE, typename POINT, typename VECT >
class VectDiracsMeasure
{
        typedef Array<POINT,1> ArrPoint;
        typedef Array<VECT,1> ArrVect;

    public:

        ArrPoint Points;
        ArrVect Vectors;

        VectDiracsMeasure() { }

        VectDiracsMeasure(int n)
        {
            Init(n);
        }

        void Init(int n)
        {
            Points.resize(Range(1,n));
            Vectors.resize(Range(1,n));
        }

        VectDiracsMeasure(ArrPoint &points)
        {
            Reference(points);
        }

        VectDiracsMeasure(ArrPoint &points, ArrVect &vectors)
        {
            Reference(points,vectors);
        }

        VectDiracsMeasure(VectDiracsMeasure &dm) : Points(Range(1,dm.Points.rows())), Vectors(Range(1,dm.Vectors.rows()))
        {
            Points = dm.Points.copy();
            Points.reindexSelf(1);
            Vectors = dm.Vectors.copy();
            Vectors.reindexSelf(1);
        }

        void operator=(VectDiracsMeasure &dm)
        {
            Init(dm.Points.rows());
            Points = dm.Points.copy();
            Points.reindexSelf(1);
            Vectors = dm.Vectors.copy();
            Vectors.reindexSelf(1);
        }

/*
// à effacer : c'est une copie pas une référence...
        void Reference(VectDiracsMeasure &dm)
        {
 	    Points.resize(Range(1,dm.Points.rows()));
            Vectors.resize(Range(1,dm.Vectors.rows()));
            Points = dm.Points.copy();
            Points.reindexSelf(1);
            Vectors = dm.Vectors.copy();
            Vectors.reindexSelf(1);
        }
*/

        void Reference(VectDiracsMeasure &mu)
        {
 	    Reference(mu.Points,mu.Vectors);
        }

        void Reference(ArrPoint &points)
        {
            Points.reference(points);
            Points.reindexSelf(1);
        }

        void Reference(ArrPoint &points, ArrVect &vectors)
        {
            Points.reference(points);
            Points.reindexSelf(1);
            Vectors.reference(vectors);
            Vectors.reindexSelf(1);
        }
        
        void Reference(VectDiracsMeasure& mu, Range r)
        {
            //Reference(mu.Points,mu.Vectors);
            Points.reference(mu.Points(r));
            Vectors.reference(mu.Vectors(r));
        }

        int NumPoints()
        {
            return Points.rows();
        }

        virtual ~VectDiracsMeasure() {};
};

template < typename TYPE, typename POINT, typename VECT >
void lincomb(VectDiracsMeasure<TYPE,POINT,VECT>& Z, VectDiracsMeasure<TYPE,POINT,VECT>& X, VectDiracsMeasure<TYPE,POINT,VECT>& Y, TYPE a, TYPE b)
{
    if(&Z==&X)
    {
        VectDiracsMeasure<TYPE,POINT,VECT> Xc(X);
        if(&Z==&Y)
        {
            VectDiracsMeasure<TYPE,POINT,VECT> Yc(Y);
            lincomb(Z,Xc,Yc,a,b);
        }
        else
            lincomb(Z,Xc,Y,a,b);
    }
    else if(&Z==&Y)
    {
        VectDiracsMeasure<TYPE,POINT,VECT> Yc(Y);
        lincomb(Z,X,Yc,a,b);
    }
    else
    {
        Z.Init(X.Points.rows()+Y.Points.rows());
        Z.Points(Range(1,X.Points.rows())) = X.Points;
        Z.Points(Range(X.Points.rows()+1,X.Points.rows()+Y.Points.rows())) = Y.Points;
        Z.Vectors(Range(1,X.Points.rows())) = a*X.Vectors;
        Z.Vectors(Range(X.Points.rows()+1,X.Points.rows()+Y.Points.rows())) = b*Y.Vectors;
    }
}

template < typename TYPE, typename POINT, typename VECT >
VectDiracsMeasure<TYPE,POINT,VECT>& operator+(VectDiracsMeasure<TYPE,POINT,VECT>& X, VectDiracsMeasure<TYPE,POINT,VECT>& Y)
{
    VectDiracsMeasure<TYPE,POINT,VECT>* out = new VectDiracsMeasure<TYPE,POINT,VECT>(X.Points.rows()+Y.Points.rows());
    lincomb(*out,X,Y,(TYPE)1.0,(TYPE)1.0);
    return *out;
}

template < typename TYPE, typename POINT, typename VECT >
VectDiracsMeasure<TYPE,POINT,VECT>& operator-(VectDiracsMeasure<TYPE,POINT,VECT>& X, VectDiracsMeasure<TYPE,POINT,VECT>& Y)
{
    VectDiracsMeasure<TYPE,POINT,VECT>* out = new VectDiracsMeasure<TYPE,POINT,VECT>(X.Points.rows()+Y.Points.rows());
    lincomb(*out,X,Y,(TYPE)1.0,(TYPE)(-1.0));
    return *out;
}

template < typename TYPE, typename POINT, typename VECT >
VectDiracsMeasure<TYPE,POINT,VECT>& operator*(TYPE a, VectDiracsMeasure<TYPE,POINT,VECT>& X)
{
    VectDiracsMeasure<TYPE,POINT,VECT>* out = new VectDiracsMeasure<TYPE,POINT,VECT>(X);
    out->Vectors *= a;
    return *out;
}







#endif // VECTDIRACSMEASURE
