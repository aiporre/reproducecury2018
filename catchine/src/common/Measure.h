#ifndef MEASURE
#define MEASURE

#include "VectKerHilbert.h"
#include "VectKerHilbert.h"
#include "Target.h"

template < typename TYPE, int DIM >
class Measure : public Target<TYPE,DIM>
{

		using Target<TYPE,DIM>::Phi;
		using Target<TYPE,DIM>::RX;
		using Target<TYPE,DIM>::TargetWeight;
		
        typedef TinyVector<TYPE,DIM> Point;
        typedef Array<Point,1> ArrPoint;
        typedef TinyVector<TYPE,DIM> VectPoint;
        typedef Array<VectPoint,1> ArrVectPoint;
        typedef TinyVector<TYPE,1> Weight;
        typedef Array<Weight,1> ArrWeight;

        VectDiracsMeasure<TYPE,Point,Weight> muY;
        VectDiracsMeasure<TYPE,Point,Weight> muPhiX;
        VectKerHilbert<TYPE,Point,Weight,VectPoint> W;
        double SqNormY;

        void SetSqNormY()
        {
            SqNormY = W.ScalProd(muY,muY);
        }

    public:

        Measure(ifstream &f, int verbosemode=1)
        {
            if(verbosemode) cout << "Reading Measure Target" << endl;
            string Tag;
            for(int i=0; i<5; i++)
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
                    W.Ker = ReadKernel<TYPE,DIM,1>(f);
                else if(!Tag.compare("Y="))
                    ReadArr(muY.Points,f);
                else if(!Tag.compare("WX,WY="))
                {
                    ReadArr(muPhiX.Vectors,f);
                    ReadArr(muY.Vectors,f);
                }
                else
				{
                    cout << "error reading file." << endl;
					throw -1;
				}
            }
            muY.Points.reindexSelf(1);
            muPhiX.Vectors.reindexSelf(1);
            muY.Vectors.reindexSelf(1);
            SetSqNormY();
			if(verbosemode) cout << "done." << endl;
        }

        void Write(ofstream &f)
        {
            f << "Measure" << endl;
            f << "Range=" << endl << RX(0) << " " << RX(0)+muPhiX.NumPoints()-1<< endl;
            f << "Weight=" << endl << TargetWeight << endl;
            f << "Kernel=" << endl;
            W.Ker->Write(f);
            f << "Y=" << endl;
            WriteArr(muY.Points,f);
            f << "WX,WY=" << endl;
            WriteArr(muPhiX.Vectors,f);
            WriteArr(muY.Vectors,f);
        }

        double Eval()
        {
            muPhiX.Reference(Phi);
            return W.ScalProd(muPhiX,muPhiX-(TYPE)2.0*muY)+SqNormY;
        }

        ArrVectPoint Gradient()
        {
            muPhiX.Reference(Phi);
            ArrVectPoint Grad(Range(1,Phi.rows()));
            Grad = 2.0 * W.GradScalProdPoints(muPhiX,muPhiX-muY);
            return Grad;
        }

};


#endif // MEASURE
