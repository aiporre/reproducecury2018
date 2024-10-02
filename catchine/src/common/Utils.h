
#ifndef UTILS
#define UTILS

#include <random/normal.h>
using namespace ranlib;

int sign(int n)
{
    return (n>0)*2-1;
}

/*
template < typename TYPE, int DIM >
class Type
{
	typedef TinyVector<TYPE,DIM> Vect;
	typedef Array<Vect,1> Vect2;
	typedef Array<Vect2,1> Vect3;
};
*/

/*
VectVect operator*(double h, const VectVect& v)
{
	VectVect hv(v.size());
	for(int i=0; i<v.size(); i++)
		hv(i) = h*v(i);
	return hv;
}
*/

/*
template < typename TYPE >
TinyVector<TYPE,3> cross(TinyVector<TYPE,3> &u, TinyVector<TYPE,3> &v)
{
    int bu = u.base(0), bv = v.base(0);
    TinyVector<TYPE,3> w(Range(1,3));
    w(1) = u(bu+1)*v(bv+2) - u(bu+2)*v(bv+1);
    w(2) = u(bu+2)*v(bv) - u(bu)*v(bv+1);
    w(3) = u(bu)*v(bv+1) - u(bu+1)*v(bv);
}
*/


double quotfact(Array<int,1> p, Array<int,1> a, Array<int,1> q, Array<int,1> b)
{
    // computes prod(p(i)!(a(i)))/prod(q(j)!(b(j)))
    double c=1.0;
	firstIndex i;
    p = max(p(i),1);
    q = max(q(i),1);
    for(int k=0; any(p>1) || any(q>1); k++)
    {
    	c *= (double)product(pow(p+0.0,1-k%a))/(double)product(pow(q+0.0,1-k%b));
    	p = max(p(i)-1,1);
    	q = max(q(i)-1,1);
    }
    return c;
}


template < typename TYPE >
TYPE det(TinyVector<TYPE,3> &u, TinyVector<TYPE,3> &v, TinyVector<TYPE,3> &w)
{
    return u(0)*v(1)*w(2) + u(1)*v(2)*w(0) + u(2)*v(0)*w(1) - u(2)*v(1)*w(0) - u(1)*v(0)*w(2) - u(0)*v(2)*w(1);
}

template < typename TYPE, int DIM >
TYPE sqnorm(TinyVector<TYPE,DIM> &u)
{
    return sum(u*u);
}

template < typename TYPE, int DIM >
Array<TinyVector<TYPE,DIM>,1> Abs(Array<TinyVector<TYPE,DIM>,1> &U)
{
	Array<TinyVector<TYPE,DIM>,1> absU(Range(U.base(0),U.base(0)+U.extent(0)-1));
    for(int t=U.base(0); t<U.base(0)+U.extent(0); t++)
		for(int i=0; i<DIM; i++)
			absU(t)(i) = abs(U(t)(i)); 
    return absU;
}

template < typename TYPE, int DIM >
Array<TinyVector<TYPE,DIM>,1> Sign(Array<TinyVector<TYPE,DIM>,1> &U)
{
	Array<TinyVector<TYPE,DIM>,1> signU(Range(U.base(0),U.base(0)+U.extent(0)-1));
    for(int t=U.base(0); t<U.base(0)+U.extent(0); t++)
		for(int i=0; i<DIM; i++)
			signU(t)(i) = (U(t)(i)>0)?1.0:-1.0; 
    return signU;
}

template < typename TYPE >
void ArrAdd(Array<TYPE,1> &A, Array<TYPE,1> &B, Array<TYPE,1> &C)
{
    for(int t=A.base(0); t<A.base(0)+A.extent(0); t++)
        A(t) = B(t) + C(t);
}

template < typename TYPE >
void ArrDiff(Array<TYPE,1> &A, Array<TYPE,1> &B, Array<TYPE,1> &C)
{
    for(int t=A.base(0); t<A.base(0)+A.extent(0); t++)
        A(t) = B(t) - C(t);
}

template < typename TYPE >
void ArrAddMult(Array<TYPE,1> &A, double m, Array<TYPE,1> &B, Array<TYPE,1> &C)
{
    for(int t=A.base(0); t<A.base(0)+A.extent(0); t++)
        A(t) = m * B(t) + C(t);
}

template < typename TYPE >
void ArrMult(Array<TYPE,1> &A, Array<TYPE,1> &B, Array<TYPE,1> &C)
{
    for(int t=A.base(0); t<A.base(0)+A.extent(0); t++)
        A(t) = B(t) * C(t);
}

template < typename TYPE, int DIM >
TYPE *contig(Array<TinyVector<TYPE,DIM>,1> &X)
{
    int nx = X.extent(firstDim);
    Array<TYPE,1> *xcont = new Array<TYPE,1>(nx*DIM);
    for(int i=0; i<nx; i++)
        for(int d=0; d<DIM; d++)
            (*xcont)(DIM*i+d) = X(i+X.base(0))(d);
    return xcont->data();
}

template < typename TYPE, int DIM >
bool iscontig(Array<TinyVector<TYPE,DIM>,1> &X)
{
    int nx = X.extent(firstDim);
	int bx = X.base(firstDim);
    bool test = 1;
    for(int i=bx; i<bx+nx; i++)
        for(int d=0; d<DIM; d++)
			test = test && ((TYPE*)X(i).data()+d == (TYPE*)X.data()+(i-bx)*DIM+d);
    return test;
}

template < typename TYPEIN, typename TYPEOUT, int DIM >
TYPEOUT *contig(Array<TinyVector<TYPEIN,DIM>,1> &X, TYPEOUT dum)
{
    (void)dum;
    int nx = X.extent(firstDim);
    Array<TYPEOUT,1> *xcont = new Array<TYPEOUT,1>(nx*DIM);
    for(int i=0; i<nx; i++)
        for(int d=0; d<DIM; d++)
            (*xcont)(DIM*i+d) = X(i+X.base(0))(d);
    return xcont->data();
}

template < typename TYPE, int DIM >
int indmax(TinyVector<TYPE,DIM> &x)
{
    int imax = 0;
    TYPE max = x(0);
    for(int i=1; i<DIM; i++)
        if(x(i)>max)
        {
            max = x(i);
            imax = i;
        }
    return imax;
}

template < typename TYPE >
int indmax(Array<TYPE,1> &x)
{
    int n = x.extent(0), imax = 1;
    TYPE max = x(1);
    for(int i=2; i<n+1; i++)
        if(x(i)>max)
        {
            max = x(i);
            imax = i;
        }
    return imax;
}

template < typename TYPE, int DIM>
Array<TinyVector<TYPE,DIM>,1> minMaxPerDim(Array<TinyVector<TYPE,DIM>,1> &X)
{
    int nx = X.extent(firstDim);
    Array<TinyVector<TYPE,DIM>,1> R(2);
    R(firstDim) = X(X.base(0));
    R(secondDim) = X(X.base(0));
    for(int i=1; i<nx; i++)
    {
        for(int d=0; d<DIM; d++)
        {
            if ( R(firstDim)(d) > X(i+X.base(0))(d) )
            {
                R(firstDim)(d) = X(i+X.base(0))(d);
            };
            if ( R(secondDim)(d) < X(i+X.base(0))(d) )
            {
                R(secondDim)(d) = X(i+X.base(0))(d);
            }
        }
    }
    return R;
}

/*
template < typename TYPEFLOAT, typename TYPE >
TYPEFLOAT fullsum(Array<TYPE,1> &A)
{
	TYPEFLOAT S;
	for(int i=A.base(0); i<A.base(0)+A.extent(0); i++)
		S += fullsum(A(i));
	return S;
}
*/

template < typename TYPE, int DIM >
double fullsum(Array<Array<TinyVector<TYPE,DIM>,1>,1> &A)
{
    double S = 0;
    for(int i=A.base(0); i<A.base(0)+A.extent(0); i++)
        for(int j=A(i).base(0); j<A(i).base(0)+A(i).extent(0); j++)
            S += sum(A(i)(j));
    return S;
}

template < typename TYPE, int DIM >
double fullsum(Array<TinyVector<TYPE,DIM>,1> &A, int start=-1, int end=-1)
{
    if(start==-1)
        start = A.base(0);
    if(end==-1)
        end = A.base(0)+A.extent(0)-1;
    double S = 0;
    for(int i=start; i<end+1; i++)
        S += sum(A(i));
    return S;
}

/*
template < typename TYPE, int DIM >
TYPE fullsum(TinyVector<TYPE,DIM> &V)
{
	return sum(V);
}
*/

template < typename TYPE, int DIM >
TYPE norminf(Array<TinyVector<TYPE,DIM>,1> &A)
{
    TYPE M = -1.0, v;
    for(int i=A.base(0); i<A.base(0)+A.extent(0); i++)
        for(int j=0; j<DIM; j++)
        {
            v = abs(A(i)(j));
            M = v>M?v:M;
	}

    return M;
}


template < typename TYPE, int DIM >
Array<TYPE,1> innersum(Array<TinyVector<TYPE,DIM>,1> &A)
{
    int base = A.base(0);
    int end = base+A.rows()-1;
    Array<TYPE,1> innS(Range(base,end));
    for(int i=base; i<end+1; i++)
        innS(i) = sum(A(i));
    return innS;
}


template < typename TYPE >
void AllocateArrArr(Array<Array<TYPE,1>,1> &Arr, int p, int q)
{
    Array<TYPE,2> tmp(Range(1,p),Range(1,q));
    Arr.resize(Range(1,p));
    for(int i=1; i<p+1; i++)
	Arr(i).reference(tmp(i,Range::all()));
        //Arr(i).resize(Range(1,q));
}

/// IO of arrays

template < typename TYPE >
int Rank(TYPE &A)
{
    (void) A;
    return 0;
}

template < typename TYPE, int DIM>
int Rank(TinyVector<TYPE,DIM> &A)
{
    (void) A;
    return 1;
}

template < typename TYPE, int DIM>
int Rank(Array<TYPE,DIM> &A)
{
    return Rank(A(A.base())) + DIM;
}

template < typename TYPE >
void WriteArrDims(TYPE &A, ofstream &fout)
{
    (void) A;
    (void) fout;
}

template < typename TYPE, int DIM>
void WriteArrDims(TinyVector<TYPE,DIM> &A, ofstream &fout)
{
    (void) A;
    (void) fout;
    fout << DIM << " ";
}

template < typename TYPE>
void WriteArrDims(Array<TYPE,3> &A, ofstream &fout)
{
    WriteArrDims(A(A.base()),fout);
    fout << A.rows() << " " << A.columns() << " " << A.depth() << " ";
}

template < typename TYPE>
void WriteArrDims(Array<TYPE,1> &A, ofstream &fout)
{
    WriteArrDims(A(A.base(0)),fout);
    fout << A.rows() << " ";
}

template < typename TYPE>
void WriteArrData(TYPE &A, ofstream &fout)
{
    fout << A << " ";
}

template < typename TYPE, int DIM>
void WriteArrData(TinyVector<TYPE,DIM> &A, ofstream &fout)
{
    for(int i=0; i<DIM; i++)
        fout << A(i) << " ";
}

template < typename TYPE >
void WriteArrData(Array<TYPE,1> &A, ofstream &fout)
{
    for(int i=A.base(0); i<A.base(0)+A.extent(0); i++)
        WriteArrData(A(i),fout);
}

template < typename TYPE>
void WriteArr(Array<TYPE,1> &A, ofstream &fout)
{
    fout << Rank(A) << endl;
    WriteArrDims(A,fout);
    fout << endl;
    WriteArrData(A,fout);
    fout << endl;
}



template < typename TYPE >
void WriteArrData(Array<TYPE,3> &A, ofstream &fout)
{
    for(int k=A.base(2); k<A.base(2)+A.extent(2); k++)
        for(int j=A.base(1); j<A.base(1)+A.extent(1); j++)
            for(int i=A.base(0); i<A.base(0)+A.extent(0); i++)
                WriteArrData(A(i,j,k),fout);
}

template < typename TYPE>
void WriteArr(Array<TYPE,3> &A, ofstream &fout)
{
    fout << Rank(A) << endl;
    WriteArrDims(A,fout);
    fout << endl;
    WriteArrData(A,fout);
    fout << endl;
}




template < typename TYPE>
void ReadArrData(TYPE &A, ifstream &fin, Array<int,1> &Dims)
{
    (void) Dims;
    fin >> A;
}


template < typename TYPE, int DIM>
void ReadArrData(TinyVector<TYPE,DIM> &A, ifstream &fin, Array<int,1> &Dims)
{
    (void) Dims;
    for(int i=0; i<DIM; i++)
        fin >> A(i);
}

template < typename TYPE >
void ReadArrData(Array<TYPE,1> &A, ifstream &fin, Array<int,1> &Dims)
{
    if(Dims(Dims.rows())!=A.rows())
        A.resize(Range(1,Dims(Dims.rows())));
    A.reindexSelf(1);
    Array<int,1> SubDims(Range(1,Dims.rows()-1));
    SubDims = Dims(Range(1,Dims.rows()-1));
    for(int i=1; i<Dims(Dims.rows())+1; i++)
        ReadArrData(A(i),fin,SubDims);
}

template < typename TYPE>
void ReadArr(Array<TYPE,1> &A, ifstream &fin)
{
    int Rank;
    fin >> Rank;
    Array<int,1> Dims(Range(1,Rank));
    for(int r=1; r<Rank+1; r++)
        fin >> Dims(r);
    ReadArrData(A,fin,Dims);
}


template < typename TYPE >
void ReadArrData(Array<TYPE,3> &A, ifstream &fin, Array<int,1> &Dims)
{
    A.resize(Dims(Dims.rows()-2),Dims(Dims.rows()-1),Dims(Dims.rows()));
    TinyVector<int,3> ones;
    ones = 1,1,1;
    A.reindexSelf(ones);
    Array<int,1> SubDims(Range(1,Dims.rows()-3));
    SubDims = Dims(Range(1,Dims.rows()-3));
    for(int k=1; k<Dims(Dims.rows())+1; k++)
        for(int j=1; j<Dims(Dims.rows()-1)+1; j++)
            for(int i=1; i<Dims(Dims.rows()-2)+1; i++)
                ReadArrData(A(i,j,k),fin,SubDims);
}

template < typename TYPE >
void ReadArr(Array<TYPE,3> &A, ifstream &fin)
{
    int Rank;
    fin >> Rank;
    Array<int,1> Dims(Range(1,Rank));
    for(int r=1; r<Rank+1; r++)
        fin >> Dims(r);
    ReadArrData(A,fin,Dims);
}

template < typename TYPE >
int partitionner(Array<TYPE,1> &tableau, int p, int r, Array<int,1> &tabind) 
{
    TYPE pivot = tableau(p);
    int i = p-1, j = r+1;
    TYPE temp;
	int tempind;
    while (1) 
	{
        do
            j--;
        while (tableau(j) > pivot);
        do
            i++;
        while (tableau(i) < pivot);
        if (i < j) 
		{
            temp = tableau(i);
            tableau(i) = tableau(j);
            tableau(j) = temp;
			tempind = tabind(i);
			tabind(i) = tabind(j);
			tabind(j) = tempind;
        }
        else
            return j;
    }
}

template < typename TYPE >
void quickSort(Array<TYPE,1> &tableau, int p, int r, Array<int,1> &tabind) 
{
    int q;
    if (p < r) 
	{
        q = partitionner(tableau, p, r, tabind);
        quickSort(tableau, p, q, tabind);
        quickSort(tableau, q+1, r, tabind);
    }
}

Array<int,1>& RandPerm(int n)
{
	Array<int,1>* randind = new Array<int,1>(Range(1,n));
	Array<double,1> tmp(Range(1,n));
	firstIndex k;
	*randind = k;
	NormalUnit<double> normunit;
	normunit.seed((unsigned int)time(0));
	for(int k=1; k<n+1; k++)
		tmp(k) = normunit.random();
	quickSort(tmp,1,n,*randind);
	return *randind;
}

#endif //UTILS


