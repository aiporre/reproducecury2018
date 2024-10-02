#ifndef LANDMARKS
#define LANDMARKS

#include "Target.h"

template < typename TYPE, int DIM >
class Landmarks : public Target<TYPE,DIM>
{

		using Target<TYPE,DIM>::Phi;
		using Target<TYPE,DIM>::RX;
		using Target<TYPE,DIM>::TargetWeight;
		
        typedef TinyVector<TYPE,DIM> Vect;
        typedef Array<Vect,1> ArrVect;
        typedef Array<TinyVector<TYPE,1>,1> ArrWeight;

        ArrVect Y;
        int NumPoints;

    public:

        Landmarks(ifstream &f, int verbosemode=1)
        {
            if(verbosemode) cout << "Reading Landmarks Target" << endl;
            string Tag;
            for(int i=0; i<3; i++)
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
                else if(!Tag.compare("Y="))
                    ReadArr(Y,f);
                else
				{
                    cout << "error reading file." << endl;
					throw -1;
				}
                NumPoints = Y.rows();
                Y.reindexSelf(1);
            }
			if(verbosemode) cout << "done." << endl;
        }

        void Write(ofstream &f)
        {
            f << "Landmarks" << endl;
            f << "Range=" << endl << RX(0) << " " << RX(0)+NumPoints-1 << endl;
            f << "Weight=" << endl << TargetWeight << endl;
            f << "Y=" << endl;
            WriteArr(Y,f);
        }

        Landmarks(ArrVect &y, Range rx, ArrWeight w, double weight)
        {
            NumPoints = y.rows();
            Y.resize(Range(1,NumPoints));
            Y = y;
            RX = rx;
            TargetWeight = weight;
        }

        double Eval()
        {
            ArrVect PhimY(Range(1,NumPoints));
            PhimY = pow2(Phi-Y);
            return fullsum(PhimY);
        }

        ArrVect Gradient()
        {
            ArrVect Grad(Range(1,NumPoints));
            Grad = 2.0*(Phi-Y);
            return Grad;
        }

};

#endif // LANDMARKS    


