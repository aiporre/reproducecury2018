
// virtual class for general kernels of RKHS spaces of the form
// gamma=K(x,y)beta where x,y are points (POINT type) and alpha,gamma vectors (VECTOR type) 
//
// the following operations are virtual :
//		Eval(x,y,beta) : computes K(x,y)beta 
//		Grad1(alpha,x,y,beta) : computes the gradient of <alpha,K(x,y)beta)> with respect to x
//		Diff1(x,y,beta,eta) : computes the partial differential of K(x,y)beta wrt x, applied to eta
//		Diff2(x,y,beta,eta) : computes the partial differential of K(x,y)beta wrt y, applied to eta
//		Grad1Diff1(alpha,x,y,beta,eta) : computes the gradient in x of <alpha,Diff1(x,y,beta,eta)>
//		Grad2Diff1(alpha,x,y,beta,eta) : computes the gradient in y of <alpha,Diff1(x,y,beta,eta)>

// the convolution operations are defined here but can be overloaded for faster implementations :
//		EvalConv(x,y,beta) : computes gamma where gamma(i) = sum_j K(x(i),y(j))beta(j)
//		Grad1Conv(alpha,x,y,beta) : computes the gradient of sum_i sum_j <alpha(i),K(x(i),y(j))beta(j))> wrt x
//		GradConv(alpha,x,beta) : computes the gradient of sum_i sum_j <alpha(i),K(x(i),x(j))beta(j))> wrt x
//		DiffConv(x,beta,eta) : computes gamma where gamma(i) = differential of sum_j K(x(i),x(j))beta(j) wrt x, applied to eta
//		GradDiffConv(x,beta,eta) : computes the gradient wrt x of the differential of sum_i sum_j <beta(i),K(x(i),x(j))beta(j))> wrt x applied to eta

#ifndef KERNEL
#define KERNEL

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include "Utils.h"

using namespace blitz;

template < typename POINT, typename VECT, typename VECTPOINT >
class Kernel
{
        typedef Array<VECT,1> ArrVect;
        typedef Array<POINT,1> ArrPoint;
        typedef Array<VECTPOINT,1> ArrVectPoint; // VECTPOINT is an element of a tangent space at a POINT

    public:

        virtual void Write(ofstream &f) = 0;

        virtual VECT Eval(POINT &x, POINT &y, VECT &beta) = 0;

        virtual VECTPOINT Grad1(VECT &alpha, POINT &x, POINT &y, VECT &beta) = 0;	
		
	//VECTPOINT Grad2(VECT &alpha, POINT &x, POINT &y, VECT &beta) { return Grad1(beta,y,x,alpha); }

        virtual VECT Diff1(POINT &x, POINT &y, VECT &beta, VECTPOINT &eta) = 0;

        virtual VECT Diff2(POINT &x, POINT &y, VECT &beta, VECTPOINT &eta) = 0;

        virtual VECTPOINT Grad1Diff1(VECT &alpha, POINT &x, POINT &y, VECT &beta, VECTPOINT &eta) = 0;

        virtual VECTPOINT Grad2Diff1(VECT &alpha, POINT &x, POINT &y, VECT &beta, VECTPOINT &eta) = 0;

        virtual ArrVect EvalConv(ArrPoint &x, ArrPoint &y, ArrVect &beta)
        {
            int nx = x.extent(firstDim);
            int ny = y.extent(firstDim);
            ArrVect gamma(Range(1,nx));
            for(int i=1; i<nx+1; i++)
            {
                gamma(i) = 0;
                for(int j=1; j<ny+1; j++)
                    gamma(i) += Eval(x(i),y(j),beta(j));
            }
            return gamma;
        }

        virtual ArrVectPoint Grad1Conv(ArrVect &alpha, ArrPoint &x, ArrPoint &y, ArrVect &beta)
        {
            int nx = x.extent(firstDim);
            int ny = y.extent(firstDim);
            ArrVectPoint gamma(Range(1,nx));
            for(int i=1; i<nx+1; i++)
            {
                gamma(i) = 0;
                for(int j=1; j<ny+1; j++)
                    gamma(i) += Grad1(alpha(i),x(i),y(j),beta(j));
            }
            return gamma;
        }
        
        virtual ArrVectPoint GradConv(ArrVect &alpha, ArrPoint &x, ArrVect &beta)
        {
            int nx = x.extent(firstDim);
            ArrVectPoint gamma(Range(1,nx));
            for(int i=1; i<nx+1; i++)
            {
                gamma(i) = 0;
                for(int j=1; j<nx+1; j++)
                    gamma(i) += Grad1(alpha(i),x(i),x(j),beta(j)) + Grad1(beta(i),x(i),x(j),alpha(j));
            }
            return gamma;
        }
        
        virtual ArrVect DiffConv(ArrPoint &x, ArrVect &beta, ArrVectPoint &eta)
        {
		int nx = x.extent(firstDim);
	    ArrVect gamma(Range(1,nx));
            for(int i=1; i<nx+1; i++)
            {
                gamma(i) = 0.0;
                for(int j=1; j<nx+1; j++)
                    gamma(i) += Diff1(x(i),x(j),beta(j),eta(i)) + Diff2(x(i),x(j),beta(j),eta(j));
            }
            return gamma;
        }
        		
        virtual ArrVectPoint GradDiffConv(ArrPoint &x, ArrVect &beta, ArrVectPoint &eta)
        {
			int nx = x.extent(firstDim);
			ArrVectPoint gamma(Range(1,nx));
            for(int i=1; i<nx+1; i++)
            {
                gamma(i) = 0.0;
                for(int j=1; j<nx+1; j++)
                    gamma(i) += Grad1Diff1(beta(i),x(i),x(j),beta(j),eta(i)) + Grad2Diff1(beta(j),x(j),x(i),beta(i),eta(j));
            }
			gamma *= 2;
            return gamma;
        }
        		
        ArrVect DeConv(ArrPoint &x, ArrVect &gamma, ArrVect &alpha)
        {
            // inversion using basic conjugate gradient scheme. input alpha is used as initial guess and gives solution at the end
            int nx = x.extent(firstDim);
            ArrVect r(Range(1,nx)), w(Range(1,nx)), rw(Range(1,nx)), wz(Range(1,nx)), rz(rw), rr(rw);
            r = Kernel::EvalConv(x,x,alpha);//EvalConv(x,x,alpha);//
            r = gamma - r;
            w = -r;
            ArrVect z(Kernel::EvalConv(x,x,w));//(EvalConv(x,x,w));//
            rw = r*w;
            wz = w*z;
            double a = fullsum(rw)/fullsum(wz);
            alpha += a*w;
            double B = 0;         
            for(int i=1; norminf(r)>1e-10; i++)//for(int i=1; i<nx+1; i++)
            {
SHOW(norminf(r))
            	r -= a*z;
            	//if(norminf(r)<1e-10)
            	//    break;
            	rz = r*z;
            	wz = w*z;
            	B = fullsum(rz)/fullsum(wz);
            	w = -r + B*w;
            	z = Kernel::EvalConv(x,x,w);//EvalConv(x,x,w);//
            	rw = r*w;
            	wz = w*z;
            	a = fullsum(rw)/fullsum(wz);
            	alpha += a*w;
            }
            return alpha;
        }
        
        virtual ~Kernel() { };
};




#endif // KERNEL


