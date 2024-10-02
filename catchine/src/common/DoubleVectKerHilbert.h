#ifndef DOUBLEVECTKERHILBERT
#define DOUBLEVECTKERHILBERT

#include "Hilbert.h"
#include "DoubleVectDiracsMeasure.h"
#include "DoubleKernel.h"
#include "Utils.h"


template < typename TYPE, typename POINT1, typename POINT2, typename VECT, typename VECTPOINT1, typename VECTPOINT2 >
class DoubleVectKerHilbert : public Hilbert<DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT> >
{
        typedef Array<POINT1,1> ArrPoint1;
        typedef Array<POINT2,1> ArrPoint2;
        typedef Array<VECT,1> ArrVect;
        typedef Array<VECTPOINT1,1> ArrVectPoint1;
        typedef Array<VECTPOINT2,1> ArrVectPoint2;
        typedef DoubleVectDiracsMeasure<TYPE,POINT1,POINT2,VECT> VectMeas;

    public:

        DoubleKernel<POINT1,POINT2,VECT,VECTPOINT1,VECTPOINT2>* Ker;

        DoubleVectKerHilbert() {}

        DoubleVectKerHilbert(DoubleKernel<POINT1,POINT2,VECT,VECTPOINT1,VECTPOINT2>* ker)
        {
            Ker = ker;
        }

        double ScalProd(VectMeas& X, VectMeas& Y)
        {
            ArrVect temp(Range(1,X.NumPoints()));
            temp = X.Vectors * Ker->EvalConv(X.Points1,X.Points2,Y.Points1,Y.Points2,Y.Vectors);     
            return fullsum(temp);
        }

        VectMeas& GradScalProd(VectMeas& X, VectMeas& Y)
        {
            ArrVectPoint1 Grad1(Range(1,X.NumPoints()));
            Grad1 = Ker->Grad1Conv(X.Vectors,X.Points1,X.Points2,Y.Points1,Y.Points2,Y.Vectors);
            ArrVectPoint2 Grad2(Range(1,X.NumPoints()));
            Grad2 = Ker->Grad2Conv(X.Vectors,X.Points1,X.Points2,Y.Points1,Y.Points2,Y.Vectors);
            ArrVect GradV(Range(1,X.NumPoints()));
            GradV = Ker->EvalConv(X.Points1,X.Points2,Y.Points1,Y.Points2,Y.Vectors);
            VectMeas* G = new VectMeas(X.NumPoints());
            G->Reference(Grad1,Grad2,GradV);
            return *G;
        }

};



#endif // DOUBLEVECTKERHILBERT
