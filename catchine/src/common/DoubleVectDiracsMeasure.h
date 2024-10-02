#ifndef DOUBLEVECTDIRACSMEASURE
#define DOUBLEVECTDIRACSMEASURE


template < typename TYPE, typename POINT1, typename POINT2, typename VECT >
class DoubleVectDiracsMeasure
{
        typedef Array<POINT1,1> ArrPoint1;
        typedef Array<POINT2,1> ArrPoint2;
        typedef Array<VECT,1> ArrVect;

    public:

        ArrPoint1 Points1;
        ArrPoint2 Points2;
        ArrVect Vectors;

        DoubleVectDiracsMeasure() { }

        DoubleVectDiracsMeasure(int n)
        {
            Init(n);
        }

        DoubleVectDiracsMeasure(ArrPoint1 &points1, ArrPoint2 &points2, ArrVect &vectors) : Points1(points1), Points2(points2), Vectors(vectors) {}

        void Init(int n)
        {
            Points1.resize(Range(1,n));
            Points2.resize(Range(1,n));
            Vectors.resize(Range(1,n));
        }

        DoubleVectDiracsMeasure(ArrPoint1 &points1, ArrPoint2 &points2)
        {
            Reference(points1,points2);
        }

        DoubleVectDiracsMeasure(DoubleVectDiracsMeasure &dm) : Points1(Range(1,dm.Points1.rows())), Points2(Range(1,dm.Points1.rows())), Vectors(Range(1,dm.Vectors.rows()))
        {
            Points1 = dm.Points1.copy();
            Points1.reindexSelf(1);
            Points2 = dm.Points2.copy();
            Points2.reindexSelf(1);
            Vectors = dm.Vectors.copy();
            Vectors.reindexSelf(1);
        }

        void operator=(DoubleVectDiracsMeasure &dm)
        {
            Init(dm.Points1.rows());
            Points1 = dm.Points1.copy();
            Points1.reindexSelf(1);
            Points2 = dm.Points2.copy();
            Points2.reindexSelf(1);
            Vectors = dm.Vectors.copy();
            Vectors.reindexSelf(1);
        }

        void Reference(ArrPoint1 &points1, ArrPoint2 &points2)
        {
            Points1.reference(points1);
            Points1.reindexSelf(1);
            Points2.reference(points2);
            Points2.reindexSelf(1);
        }

        void Reference(ArrPoint1 &points1, ArrPoint2 &points2, ArrVect &vectors)
        {
            Reference(points1,points2);
            Vectors.reference(vectors);
            Vectors.reindexSelf(1);
        }

        int NumPoints()
        {
            return Points1.rows();
        }

        virtual ~DoubleVectDiracsMeasure() {};
};

template < typename TYPE, typename POINT1, typename POINT2, typename VECT >
void lincomb(DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>& Z, DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>& X,
             DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>& Y, TYPE a, TYPE b)
{
    if(&Z==&X)
    {
        DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT> Xc(X);
        if(&Z==&Y)
        {
            DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT> Yc(Y);
            lincomb(Z,Xc,Yc,a,b);
        }
        else
            lincomb(Z,Xc,Y,a,b);
    }
    else if(&Z==&Y)
    {
        DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT> Yc(Y);
        lincomb(Z,X,Yc,a,b);
    }
    else
    {
        Z.Init(X.Points1.rows()+Y.Points1.rows());
        Z.Points1(Range(1,X.Points1.rows())) = X.Points1;
        Z.Points2(Range(1,X.Points1.rows())) = X.Points2;
        Z.Points1(Range(X.Points1.rows()+1,X.Points1.rows()+Y.Points1.rows())) = Y.Points1;
        Z.Points2(Range(X.Points1.rows()+1,X.Points1.rows()+Y.Points1.rows())) = Y.Points2;
        Z.Vectors(Range(1,X.Points1.rows())) = a*X.Vectors;
        Z.Vectors(Range(X.Points1.rows()+1,X.Points1.rows()+Y.Points1.rows())) = b*Y.Vectors;
    }
}

template < typename TYPE, typename POINT1, typename POINT2, typename VECT >
DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>& operator+(DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>& X, DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>& Y)
{
    DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>* out = new DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>(X.Points1.rows()+Y.Points1.rows());
    lincomb(*out,X,Y,(TYPE)1.0,(TYPE)1.0);
    return *out;
}

template < typename TYPE, typename POINT1, typename POINT2, typename VECT >
DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>& operator-(DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>& X, DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>& Y)
{
    DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>* out = new DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>(X.Points1.rows()+Y.Points1.rows());
    lincomb(*out,X,Y,(TYPE)1.0,(TYPE)(-1.0));
    return *out;
}

template < typename TYPE, typename POINT1, typename POINT2, typename VECT >
DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>& operator*(TYPE a, DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>& X)
{
    DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>* out = new DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT>(X);
    out->Vectors *= a;
    return *out;
}





#endif // DOUBLEVECTDIRACSMEASURE
