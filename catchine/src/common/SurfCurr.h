#ifndef SURFCURR
#define SURFCURR

#include "VectKerHilbert.h"
#include "VectDiracsMeasure.h"
#include "Target.h"

template < typename TYPE >
class SurfCurr : public Target<TYPE,3>
{

		using Target<TYPE,3>::Phi;
		using Target<TYPE,3>::RX;
		using Target<TYPE,3>::TargetWeight;
		
        typedef TinyVector<TYPE,3> Point;
        typedef Array<Point,1> ArrPoint;
        typedef TinyVector<TYPE,3> Vect;
        typedef Array<Vect,1> ArrVect;
        typedef TinyVector<TYPE,3> VectPoint;
        typedef Array<VectPoint,1> ArrVectPoint;
        typedef TinyVector<int,3> Face;
        typedef Array<Face,1> ArrFace;
        typedef Array<TYPE,1> ArrWeight;

        VectDiracsMeasure<TYPE,Point,Vect> muY, muPhiX;
        ArrVect VY, CPhiX, NPhiX;
        ArrFace FY, FX;
        ArrWeight WX, WY;
        VectKerHilbert<TYPE,Point,Vect,VectPoint> H;
        double SqNormY, MatchingPursuitEpsilon;

        void Set_muY_and_SqNormY()
        {
            ComputeCenters(muY.Points,VY,FY);
            ComputeNormals(muY.Vectors,VY,FY,WY);
            SqNormY = H.ScalProd(muY,muY);
        }

        void ComputeCenters(ArrVect& C, ArrVect &V, ArrFace &F)
        {
            TYPE oot = (TYPE)1.0/(TYPE)3.0;
            for(int i=1; i<F.rows()+1; i++)
                C(i) = oot*(V(F(i)(0))+V(F(i)(1))+V(F(i)(2)));
        }

        void ComputeNormals(ArrVect &N, ArrVect &V, ArrFace &F, ArrWeight &W)
        {
            Vect u,v,w;
            for(int i=1; i<F.rows()+1; i++)
            {
                u = V(F(i)(1))-V(F(i)(0));
                v = V(F(i)(2))-V(F(i)(0));
                w = cross(u,v);
                N(i) = (W(i) * 0.5) * w;
            }
        }

        ArrVect TransDiffCenters(ArrVect &eC, ArrVect &V, ArrFace &F)
        {
            // compute transpose of derivative of 'vertices to centers' function
            TYPE oot = (TYPE)1.0/(TYPE)3.0;
            ArrVect eV(Range(1,V.rows()));
            eV = 0.0;
            for(int i=1; i<F.rows()+1; i++)
            {
                eV(F(i)(0)) += oot*eC(i);
                eV(F(i)(1)) += oot*eC(i);
                eV(F(i)(2)) += oot*eC(i);
            }
            return eV;
        }

        ArrVect TransDiffNormals(ArrVect &eN, ArrVect &V, ArrFace &F, ArrWeight &W)
        {
            // compute transpose of derivative of 'vertices to normals' function
            ArrVect eV(Range(1,V.rows()));
            eV = 0.0;
            Vect u,v,w;
            for(int i=1; i<F.rows()+1; i++)
            {
                u = V(F(i)(1))-V(F(i)(0));
                u = .5*W(i)*cross(u,eN(i));
                v = V(F(i)(2))-V(F(i)(0));
                v = .5*W(i)*cross(v,eN(i));
                w = u-v;
                eV(F(i)(0)) += w;
                eV(F(i)(1)) += v;
                eV(F(i)(2)) -= u;
            }
            return eV;
        }

    public:

        SurfCurr(ifstream &f, int verbosemode=1)
        {
            if(verbosemode) cout << "Reading Surface Current Target" << endl;
            string Tag;
            for(int i=0; i<8; i++)
            {
                f >> Tag;
                if(!Tag.compare("Range="))
                {
                    int debut, fin;
                    f >> debut >> fin;
                    RX = Range(debut,fin);
                }
                else if(!Tag.compare("Weight="))
                    f >> TargetWeight;
                else if(!Tag.compare("Kernel="))
                    H.Ker = ReadKernel<TYPE,3,3>(f);
                else if(!Tag.compare("VY="))
                    ReadArr(VY,f);
                else if(!Tag.compare("FY="))
                    ReadArr(FY,f);
                else if(!Tag.compare("FX="))
                    ReadArr(FX,f);
                else if(!Tag.compare("WX,WY="))
                {
                    ReadArr(WX,f);
                    ReadArr(WY,f);
                }
                else if(!Tag.compare("MatchingPursuitEpsilon="))
                    f >> MatchingPursuitEpsilon;
                else
				{
                    cout << "error reading file." << endl;
					throw -1;
				}

            }
            VY.reindexSelf(1);
            FY.reindexSelf(1);
            FX.reindexSelf(1);
            WX.reindexSelf(1);
            WY.reindexSelf(1);
            CPhiX.resize(Range(1,FX.rows()));
            NPhiX.resize(Range(1,FX.rows()));
            muPhiX.Reference(CPhiX,NPhiX);
            muY.Init(FY.rows());
            Set_muY_and_SqNormY();
            if(MatchingPursuitEpsilon)
            {
                cout << "Using Matching Pursuit Approx of Target" << endl;
            	muY.Reference(H.ApproxRandom(muY,MatchingPursuitEpsilon));
             	ComputeVerticesFacesAndWeights(muY.Points,muY.Vectors,VY,FY,WY);
            	Set_muY_and_SqNormY();
            	SqNormY = H.ScalProd(muY,muY);
            }
			if(verbosemode) cout << "done." << endl;
            	
        }
        
        void ComputeVerticesFacesAndWeights(ArrPoint &C, ArrVect &N, ArrPoint &V, ArrFace &F, ArrWeight &W)
        {
            int n = C.rows();
            TYPE nNi;
            Vect Nin, absNin, v1, v2, v3, e1, e2, tmp;
            TYPE lambda;
            V.resize(Range(1,3*n));
            F.resize(Range(1,n));
            W.resize(Range(1,n));
            // compute min distance between all centers
            Vect Cij = C(1)-C(2);
            TYPE mindist = n==1?0:sqnorm(Cij);
            for(int i=1; i<n+1; i++)
            	for(int j=i+1; j<n+1; j++)
            	{
            	    Cij = C(i)-C(j);
            	    mindist = min(mindist,sqnorm(Cij));
            	}
            mindist = sqrt(mindist);
            lambda = mindist*.25;
            TYPE sqrt3o2 = sqrt(3.0)*.5;
            TYPE coef = 4.0/(3.0*sqrt(3.0)*lambda*lambda);
            for(int i=1; i<n+1; i++)
            {
            	nNi = norm(N(i));
            	Nin = N(i)/nNi;
            	tmp = Nin;
            	absNin = abs(Nin);
            	tmp(indmax(absNin)) = 0;
            	tmp = tmp/norm(tmp);
            	e1 = cross(Nin,tmp);
            	e1 = e1/norm(e1);
            	e2 = cross(Nin,e1);
            	e2 = e2/norm(e2);
            	V(3*i-2) = C(i) + lambda * e1;
            	V(3*i-1) = C(i) + lambda * (-.5*e1+sqrt3o2*e2);
            	V(3*i) = C(i) + lambda * (-.5*e1-sqrt3o2*e2);
            	F(i) = 3*i-2,3*i-1,3*i;
            	W(i) = nNi*coef;
            }
        }

        void Write(ofstream &f)
        {
            f << "SurfCurr" << endl;
            f << "Range=" << endl << RX(0) << " " << RX(1)<< endl;
            f << "Weight=" << endl << TargetWeight << endl;
            f << "Kernel=" << endl;
            H.Ker->Write(f);
            f << "VY=" << endl;
            WriteArr(VY,f);
            f << "FX=" << endl;
            WriteArr(FX,f);
            f << "FY=" << endl;
            WriteArr(FY,f);
            f << "WX,WY=" << endl;
            WriteArr(WX,f);
            WriteArr(WY,f);
        }

        double Eval()
        {
            ComputeCenters(CPhiX,Phi,FX);
            ComputeNormals(NPhiX,Phi,FX,WX);
            return H.ScalProd(muPhiX,muPhiX-(TYPE)2.0*muY)+SqNormY;
        }

        ArrVect Gradient()
        {
            ComputeCenters(CPhiX,Phi,FX);
            ComputeNormals(NPhiX,Phi,FX,WX);
            ArrVect temp1(2.0*H.GradScalProdPoints(muPhiX,muPhiX-muY));
            ArrVect temp2(2.0*H.GradScalProdVectors(muPhiX,muPhiX-muY));
            ArrVect Grad(TransDiffCenters(temp1,Phi,FX) + TransDiffNormals(temp2,Phi,FX,WX));
            return Grad;
        }

};


#endif // SURFCURR
