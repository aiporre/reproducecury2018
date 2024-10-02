#ifndef VECTKERHILBERT
#define VECTKERHILBERT

#include "Hilbert.h"
#include "VectDiracsMeasure.h"
#include "Kernel.h"
#include "Utils.h"

template < typename TYPE, typename POINT, typename VECT, typename VECTPOINT >
class VectKerHilbert : public Hilbert<VectDiracsMeasure<TYPE,POINT,VECT> >
{
        typedef Array<POINT,1> ArrPoint;
        typedef Array<VECT,1> ArrVect;
        typedef Array<VECTPOINT,1> ArrVectPoint;
        typedef VectDiracsMeasure<TYPE,POINT,VECT> VectMeas;
        
    public:

        Kernel<POINT,VECT,VECTPOINT>* Ker;

        VectKerHilbert() {}

        VectKerHilbert(Kernel<POINT,VECT,VECTPOINT>* ker)
        {
            Ker = ker;
        }

        double ScalProd(VectMeas& X, VectMeas& Y)
        {
            ArrVect temp(Range(1,Y.NumPoints()));
            temp = Y.Vectors * Ker->EvalConv(Y.Points,X.Points,X.Vectors);
            return fullsum(temp);
        }

        ArrVectPoint GradScalProdPoints(VectMeas& X, VectMeas& Y) // gradient wrt X.Points
        {
            ArrVectPoint Grad(Range(1,X.NumPoints()));
            Grad = Ker->Grad1Conv(X.Vectors,X.Points,Y.Points,Y.Vectors);
            return Grad;
        }

        ArrVect GradScalProdVectors(VectMeas& X, VectMeas& Y) // gradient wrt X.Vectors
        {
            ArrVect Grad(Range(1,X.NumPoints()));
            Grad = Ker->EvalConv(X.Points,Y.Points,Y.Vectors);
            return Grad;
        }

        ArrVect EvalVectorField(VectMeas& X, ArrPoint& y)
        {
            ArrVect gamma(Range(1,y.rows()));
            gamma = Ker->EvalConv(y,X.Points,X.Vectors);
            return gamma;
        }

        ArrVect EvalVectorFieldDirect(VectMeas& X, ArrPoint& y)
        {
            ArrVect gamma(Range(1,y.rows()));
            gamma = Ker->Kernel<POINT,VECT,VECTPOINT>::EvalConv(y,X.Points,X.Vectors);
            return gamma;
        }

/*        VectMeas& ApproxMP(VectMeas& mu, TYPE epsilon)
        {
            TYPE valuestart, residue, residue2, invk0 = 1.0; // works only with scalar kernel s.t. k(0)=1 !!!
            ArrVect Gamma(EvalVectorField(mu,mu.Points));
            int ind, indprev, i, n = mu.NumPoints();
            ArrVect Gamma2(Range(1,n)), tmp_arr(Range(1,n)), Alpha(Range(1,n));
            Alpha = mu.Vectors;
            VectMeas* nu = new VectMeas(n);
            nu->Vectors = 0.0;
            Array<TYPE,1> sqnrmGamma(Range(1,n));
            VectMeas tmp(1);
            i = 1;
		    tmp_arr = Gamma*mu.Vectors;
		    valuestart = fullsum(tmp_arr);
		    residue = valuestart;
		    residue2 = valuestart;
int cnt = 0;
            while(1)
            {
cnt++;
                Gamma2 = Gamma*Gamma;
            	sqnrmGamma = innersum(Gamma2);
            	ind = indmax(sqnrmGamma);
SHOW(residue)
SHOW(residue2)
            	if(residue<epsilon*valuestart)
            	{
            	    i--;
            	    break;
            	}
            	//if(tab(ind))
            	//    iw = tab(ind);
            	//else
            	//{
            	    iw = i;
            	    i++;
            	//}
            	//tab(ind) = iw;
            	
            	nu->Points(iw) = mu.Points(ind(i));
            	//nu->Vectors(iw) += Gamma(ind)*invk0;
            	nu->Vectors = GaussSeidel(tmpmat,Gamma0(ind)));
            	//tmp.Points(1) = mu.Points(ind);
            	//tmp.Vectors(1) = Gamma(ind)*invk0;
            	Gamma -= EvalVectorFieldDirect(tmp,mu.Points);
            	Alpha(ind) -= Gamma(ind)*invk0;
            	tmp_arr = Gamma*Alpha;
            	residue = fullsum(tmp_arr);
            	residue2 -= invk0*sqnrmGamma(ind);
            	indprev = ind;
            }
            nu->Points.resizeAndPreserve(i);
            nu->Points.reindexSelf(1);
            nu->Vectors.resizeAndPreserve(i);
            nu->Vectors.reindexSelf(1);
SHOW(n)
SHOW(nu->NumPoints())
SHOW(cnt)
            return *nu;
        }
*/

        VectMeas& ApproxMP(VectMeas& mu, TYPE epsilon)
        {
        	int n = mu.NumPoints();
SHOW(n)
        	VectMeas nu0(n), *nu = new VectMeas;
        	ArrVect tmp(Range(1,n));
        	ArrVect Gamma0(EvalVectorField(mu,mu.Points));
        	ArrVect Gamma(Range(1,n));
        	Gamma = Gamma0;
        	ArrVect Gamma2(Range(1,n));
        	ArrVect Gammanu0(Range(1,n)), Gammanu;
        	tmp = Gamma*mu.Vectors;
        	double valinit = sqrt(fullsum(tmp));
        	Array<TYPE,1> sqnrmGamma(Range(1,n));
		int i = 1;
		int ind;
		while(1)
		{
SHOW(i)
			Gamma2 = Gamma*Gamma;
            		sqnrmGamma = innersum(Gamma2);
            		ind = indmax(sqnrmGamma);
			nu0.Points(i) = mu.Points(ind);
			nu->Reference(nu0,Range(1,i));
			Gammanu0(i) = Gamma0(ind);
			Gammanu.reference(Gammanu0(Range(1,i)));
			Ker->DeConv(nu->Points,Gammanu,nu->Vectors);
			Gamma = Gamma0 - EvalVectorField(*nu,mu.Points);
			tmp(Range(1,i)) = Gammanu*nu->Vectors;
SHOW((valinit-sqrt(fullsum(tmp,1,i)))/valinit)
			if(valinit-sqrt(fullsum(tmp,1,i))<epsilon*valinit)
				break;			
			i++;
		}
		return *nu;
        }
	
	VectMeas& ApproxRandom(VectMeas& mu, TYPE epsilon)
	{
		int n = mu.NumPoints();
		SHOW(n)
		VectMeas munu(n), mushuffle(n), *nu = new VectMeas;
		ArrVect tmp(Range(1,n));
		ArrVect Gammamu(EvalVectorField(mu,mu.Points));
		//ArrVect Gamma(Range(1,n));
		//Gamma = Gamma0;
		//ArrVect Gamma2(Range(1,n));
		ArrVect Gammamunu(Range(1,n)), Gammamushuffle(Range(1,n)), Gammanu;
		tmp = Gammamu*mu.Vectors;
		double valinit = sqrt(fullsum(tmp));
		SHOW(valinit)
		//Array<TYPE,1> sqnrmGamma(Range(1,n));
		int i = 100;
		//int ind;
		Array<int,1> randind(RandPerm(n));
		for(int k=1; k<n+1; k++)
		{
			randind(k) = k;
			mushuffle.Points(k) = mu.Points(randind(k));
			mushuffle.Vectors(k) = mu.Vectors(randind(k));
			Gammamushuffle(k) = Gammamu(randind(k));
		}
		munu.Points = mushuffle.Points;
		munu.Vectors = mushuffle.Vectors;
		while(i<150)
		{
			SHOW(i)
			//Gamma2 = Gamma*Gamma;
			//sqnrmGamma = innersum(Gamma2);
			//ind = indmax(sqnrmGamma);		
			nu->Points.reference(mushuffle.Points(Range(1,i)));
			nu->Vectors.resize(Range(1,i));
			Gammanu.reference(Gammamushuffle(Range(1,i)));
			Ker->DeConv(nu->Points,Gammanu,nu->Vectors);
ArrVect ess(Range(1,i));
ess = Gammamushuffle(Range(1,i)) - EvalVectorField(*nu,nu->Points);
SHOW(ess(Range(1,i)))
SHOW(Gammamushuffle(Range(1,i)))			
			munu.Vectors(Range(1,i)) = mushuffle.Vectors(Range(1,i)) - nu->Vectors;
			Gammamunu = EvalVectorField(munu,munu.Points);
SHOW(Gammamunu(Range(1,100)))
			tmp = Gammamunu*munu.Vectors;
			SHOW((sqrt(fullsum(tmp)))/valinit)
			if(sqrt(fullsum(tmp))<epsilon*valinit)
				break;			
			i+=100;
			if(i>n)
				i = n;
		}
		return *nu;
	}

};

#endif // VECTKERHILBERT
