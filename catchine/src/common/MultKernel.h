
#ifndef MULTKERNEL
#define MULTKERNEL

#include "DoubleKernel.h"

template < typename TYPE, typename POINT1, typename POINT2, typename VECT, typename VECTPOINT1, typename VECTPOINT2 >
class MultKernel : public DoubleKernel< POINT1, POINT2, VECT, VECTPOINT1, VECTPOINT2 >
{
        typedef Kernel<POINT1,VECT,VECTPOINT1> KerType1;
        typedef Kernel<POINT2,VECT,VECTPOINT2> KerType2;

    protected:

        KerType1 *Ker1;
        KerType2 *Ker2;

    public:

        MultKernel(KerType1 *ker1, KerType2 *ker2) : Ker1(ker1), Ker2(ker2) { }

        void Write(ofstream &f)
        {
            f << "Mult,kernel1=" << endl;
            Ker1->Write(f);
            f << "kernel2=" << endl;
            Ker2->Write(f);
        }

        VECT Eval(POINT1 &x1, POINT2 &x2, POINT1 &y1, POINT2 &y2, VECT &beta)
        {
            VECT tmp = Ker2->Eval(x2,y2,beta);
            return Ker1->Eval(x1,y1,tmp);
        }

        VECTPOINT1 Grad1(VECT &alpha, POINT1 &x1, POINT2 &x2, POINT1 &y1, POINT2 &y2, VECT &beta)
        {
            VECT tmp = Ker2->Eval(x2,y2,beta);
            return Ker1->Grad1(alpha,x1,y1,tmp);
        }

        VECTPOINT2 Grad2(VECT &alpha, POINT1 &x1, POINT2 &x2, POINT1 &y1, POINT2 &y2, VECT &beta)
        {
            VECT tmp = Ker1->Eval(y1,x1,alpha);
            return Ker2->Grad1(tmp,x2,y2,beta);
        }

};


#endif // MULTKERNEL
