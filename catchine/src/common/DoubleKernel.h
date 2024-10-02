#ifndef DOUBLEKERNEL
#define DOUBLEKERNEL

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

using namespace blitz;

template < typename POINT1, typename POINT2, typename VECT, typename VECTPOINT1, typename VECTPOINT2 >
class DoubleKernel
{
        typedef Array<VECT,1> ArrVect;
        typedef Array<POINT1,1> ArrPoint1;
        typedef Array<POINT2,1> ArrPoint2;
        typedef Array<VECTPOINT1,1> ArrVectPoint1;
        typedef Array<VECTPOINT2,1> ArrVectPoint2;

    public:

        virtual void Write(ofstream &f) = 0;

        virtual VECT Eval(POINT1 &x1, POINT2 &x2, POINT1 &y1, POINT2 &y2, VECT &beta) = 0;

        virtual VECTPOINT1 Grad1(VECT &alpha, POINT1 &x1, POINT2 &x2, POINT1 &y1, POINT2 &y2, VECT &beta) = 0;

        virtual VECTPOINT2 Grad2(VECT &alpha, POINT1 &x1, POINT2 &x2, POINT1 &y1, POINT2 &y2, VECT &beta) = 0;

        virtual ArrVect EvalConv(ArrPoint1 &x1, ArrPoint2 &x2, ArrPoint1 &y1, ArrPoint2 &y2, ArrVect &beta)
        {
            int nx = x1.extent(firstDim);
            int ny = y1.extent(firstDim);
            ArrVect gamma(Range(1,nx));
            for(int i=1; i<nx+1; i++)
            {
                gamma(i) = 0;
                for(int j=1; j<ny+1; j++)
                    gamma(i) += Eval(x1(i),x2(i),y1(j),y2(j),beta(j));
            }
            return gamma;
        }

        virtual ArrVectPoint1 Grad1Conv(ArrVect &alpha, ArrPoint1 &x1, ArrPoint2 &x2, ArrPoint1 &y1, ArrPoint2 &y2, ArrVect &beta)
        {
            int nx = x1.extent(firstDim);
            int ny = y1.extent(firstDim);
            ArrVectPoint1 gamma(Range(1,nx));
            for(int i=1; i<nx+1; i++)
            {
                gamma(i) = 0;
                for(int j=1; j<ny+1; j++)
                    gamma(i) += Grad1(alpha(i),x1(i),x2(i),y1(j),y2(j),beta(j));
            }
            return gamma;
        }

        virtual ArrVectPoint2 Grad2Conv(ArrVect &alpha, ArrPoint1 &x1, ArrPoint2 &x2, ArrPoint1 &y1, ArrPoint2 &y2, ArrVect &beta)
        {
            int nx = x1.extent(firstDim);
            int ny = y1.extent(firstDim);
            ArrVectPoint2 gamma(Range(1,nx));
            for(int i=1; i<nx+1; i++)
            {
                gamma(i) = 0;
                for(int j=1; j<ny+1; j++)
                    gamma(i) += Grad2(alpha(i),x1(i),x2(i),y1(j),y2(j),beta(j));
            }
            return gamma;
        }

        virtual ~DoubleKernel() {};
};


#endif // DOUBLEKERNEL
