#ifndef DIRECTKERNEL
#define DIRECTKERNEL

#include "Kernel.h"

template < typename TYPE, int DIMPOINT, int DIMVECT >
class DirectKernel : public Kernel<TYPE,DIMPOINT,DIMVECT>
{	
	typedef TinyVector<TYPE,DIMVECT> Vect;
	typedef Array<Vect,1> ArrVect;
	typedef TinyVector<TYPE,DIMPOINT> Point;
	typedef Array<Point,1> ArrPoint;
	typedef TinyVector<TYPE,DIMPOINT> VectPoint;
	typedef Array<VectPoint,1> ArrVectPoint;
	typedef TinyMatrix<TYPE,DIMVECT,DIMVECT> Matrix; 
		
    public:
		
	virtual void Write(ofstream &f) = 0;
	virtual Vect Eval(Point &x, Point &y, Vect &alpha) = 0;
        virtual Vect Diff(Point &x, VectPoint &beta, Point &y, Vect &alpha) = 0;
        virtual VectPoint Grad(Vect &gamma, Point &x, Point &y, Vect &alpha) = 0;
	VectPoint SpecDiffTKer(Point &x, Point &y, Vect &alphax, Vect &alphay, Vect &etax, Vect &etay)
        {
            return Grad(etay,x,y,alphax) + Grad(alphay,x,y,etax);
        }
		
        ArrVect ConvKer(ArrVectPoint &x, ArrVectPoint &y, ArrVect &alpha)
        {
            int nx = x.extent(firstDim);
            int ny = y.extent(firstDim);
            ArrVect gamma(Range(1,ny));
            for(int j=1;j<ny+1;j++)
            {
                gamma(j) = 0;
                for(int i=1;i<nx+1;i++)
                    gamma(j) += Eval(x(i),y(j),alpha(i));
            }
            return gamma;
        }
       
        ArrVectPoint ConvGradKer(ArrVect &gamma, ArrVectPoint &x, ArrVectPoint &y, ArrVect &alpha)
        {
            int nx = x.extent(firstDim);
            int ny = y.extent(firstDim);
            ArrVectPoint beta(Range(1,nx));
            for(int i=1;i<nx+1;i++)
            {
                beta(i) = 0;
                for(int j=1;j<ny+1;j++)
                    beta(i) += Grad(gamma(j),x(i),y(j),alpha(i));
            }
            return beta;
        }

        ArrVectPoint ConvSpecDiffTKer(ArrPoint &x, ArrVect &alpha, ArrVect &eta)
        {
        	int nx = x.extent(firstDim);
			ArrVectPoint gamma(Range(1,nx));
            
			for(int i=1;i<nx+1;i++)
				{
        	        gamma(i) = 0.0;
        	        for(int j=1;j<nx+1;j++)
						gamma(i) += SpecDiffTKer(x(i),x(j),alpha(i),alpha(j),eta(i),eta(j));
				}
			return gamma;
        }
        
        virtual ~DirectKernel() {};
};


#endif  // DIRECTKERNEL

        
