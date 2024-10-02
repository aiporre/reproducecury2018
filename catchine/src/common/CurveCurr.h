#ifndef CURVECURR
#define CURVECURR

#include "VectKerHilbert.h"
#include "Target.h"

template < typename TYPE, int DIM >
class CurveCurr : public Target<TYPE,DIM>
{

		using Target<TYPE,DIM>::Phi;
		using Target<TYPE,DIM>::RX;
		using Target<TYPE,DIM>::TargetWeight;
		
        typedef TinyVector<TYPE,DIM> Point;
        typedef Array<Point,1> ArrPoint;
        typedef TinyVector<TYPE,DIM> Vect;
        typedef Array<Vect,1> ArrVect;
        typedef TinyVector<TYPE,DIM> VectPoint;
        typedef Array<VectPoint,1> ArrVectPoint;
        typedef TinyVector<int,2> Face;
        typedef Array<Face,1> ArrFace;
        typedef Array<TYPE,1> ArrWeight;

        VectDiracsMeasure<TYPE,Point,Vect> muY, muPhiX;
        ArrPoint VY, CPhiX;
        ArrVect TPhiX;
        ArrFace FY, FX;
        ArrWeight WX, WY;
        VectKerHilbert<TYPE,Point,Vect,VectPoint> H;
        double SqNormY, MatchingPursuitEpsilon;

        void Set_muY_and_SqNormY()
        {
            ComputeCenters(muY.Points,VY,FY);
            ComputeTangents(muY.Vectors,VY,FY,WY);
            SqNormY = H.ScalProd(muY,muY);
        }

        void ComputeCenters(ArrPoint& C, ArrPoint &V, ArrFace& F)
        {
            for(int i=1; i<F.rows()+1; i++)
                C(i) = .5*(V(F(i)(0))+V(F(i)(1)));
        }

        void ComputeTangents(ArrVect &T, ArrPoint &V, ArrFace &F, ArrWeight &W)
        {
            for(int i=1; i<F.rows()+1; i++)
                T(i) = W(i) * (V(F(i)(1))-V(F(i)(0)));
        }

        ArrVectPoint TransDiffCenters(ArrVectPoint &eC, ArrPoint &V, ArrFace &F)
        {
            // compute transpose of derivative of 'vertices to centers' function
            ArrVectPoint eV(Range(1,V.rows()));
            eV = 0.0;
            for(int i=1; i<F.rows()+1; i++)
            {
                eV(F(i)(0)) += .5*eC(i);
                eV(F(i)(1)) += .5*eC(i);
            }
            return eV;
        }

        ArrVectPoint TransDiffTangents(ArrVect &eT, ArrPoint &V, ArrFace &F, ArrWeight &W)
        {
            // compute transpose of derivative of 'vertices to tangents' function
            ArrVectPoint eV(Range(1,V.rows()));
            eV = 0.0;
            for(int i=1; i<F.rows()+1; i++)
            {
                eV(F(i)(0)) -= W(i) * eT(i);
                eV(F(i)(1)) += W(i) * eT(i);
            }
            return eV;
        }

    public:

        CurveCurr(ifstream &f, int verbosemode = 1)
        {
            if(verbosemode) cout << "Reading Curve Current Target" << endl;
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
                    H.Ker = ReadKernel<TYPE,DIM,DIM>(f);
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
                    cout << "error reading Curve Current Target." << endl;
					throw -1;
				}

            }
            VY.reindexSelf(1);
            FY.reindexSelf(1);
            FX.reindexSelf(1);
            WX.reindexSelf(1);
            WY.reindexSelf(1);
            CPhiX.resize(Range(1,FX.rows()));
            TPhiX.resize(Range(1,FX.rows()));
            muPhiX.Reference(CPhiX,TPhiX);
            muY.Init(FY.rows());
            Set_muY_and_SqNormY();
            if(MatchingPursuitEpsilon)
            {
                if(verbosemode) cout << "Using Matching Pursuit Approx of Target" << endl;
            	muY.Reference(H.ApproxMP(muY,MatchingPursuitEpsilon));
             	//ComputeVerticesFacesAndWeights(muY.Points,muY.Vectors,VY,FY,WY);
            	//Set_muY_and_SqNormY();
            	SqNormY = H.ScalProd(muY,muY);
            }
			if(verbosemode) cout << "done." << endl;
        }

        void Write(ofstream &f)
        {
            f << "CurveCurr" << endl;
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
            ComputeTangents(TPhiX,Phi,FX,WX);
            return H.ScalProd(muPhiX,muPhiX-(TYPE)2.0*muY)+SqNormY;
        }

        ArrVectPoint Gradient()
        {
            ComputeCenters(CPhiX,Phi,FX);
            ComputeTangents(TPhiX,Phi,FX,WX);
            ArrVectPoint temp1(2.0*H.GradScalProdPoints(muPhiX,muPhiX-muY));
            ArrVect temp2(2.0*H.GradScalProdVectors(muPhiX,muPhiX-muY));
            ArrVectPoint Grad(TransDiffCenters(temp1,Phi,FX) + TransDiffTangents(temp2,Phi,FX,WX));
            return Grad;
        }

};


#endif // CURVECURR
