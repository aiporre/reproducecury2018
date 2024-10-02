#ifndef CURVECYCLE
#define CURVECYCLE

#include "DoubleVectKerHilbert.h"
#include "Target.h"
#include "CycleFunctions.h"
#include "CycleFunctions3D.h"
#include "LookUpFunction1D.h"
#include "AngleKernel.h"
#include "SqDistScalarKernel.h"
#include "MultKernel.h"
#include <cmath>

template < typename TYPE, int DIM >
class CurveCycle
{
        typedef TinyVector<TYPE,DIM> Point;
        typedef Array<Point,1> ArrPoint;

        typedef TinyVector<TYPE,DIM> VectPoint;
        typedef Array<VectPoint,1> ArrVectPoint;

        typedef TinyVector<int,2> Face;
        typedef Array<Face,1> ArrFace;

        typedef TinyVector<TYPE,1> Weight;
        typedef Array<Weight,1> ArrWeight;

        typedef TinyVector<TYPE,DIM> Cone;
        typedef Array<Cone,1> ArrCone;

        typedef TinyVector<TYPE,DIM> VectCone;
        typedef Array<VectCone,1> ArrVectCone;

        typedef DoubleVectDiracsMeasure<TYPE,Point,Cone,Weight> Cycle;
        typedef DoubleVectDiracsMeasure<TYPE,VectPoint,VectCone,Weight> VectCycle;

        int NumVertices, NumFaces;
        ArrWeight WeightVertex, WeightFace;
        ArrFace Faces;
        TYPE PI;

    public :

        Cycle CycleVertex, CycleFace;

        CurveCycle() {}

        void Init(ArrFace &face)
        {
            PI = 3.14159265358979323846;
            Faces.reference(face);
            Faces.reindexSelf(1);
            NumFaces = Faces.rows();
            CycleVertex.Init(2*NumFaces);
            CycleFace.Init(NumFaces);
        }

        void Compute(ArrPoint &Vertex)
        {
            int NumVertices = Vertex.rows();
            Point v0, v1, c;
            VectPoint tau;
            Cone taun;
            TYPE l;
            Array<bool,1> IsNode(Range(1,Vertex.rows()));
            IsNode = 0;
            for(int i=1; i<NumFaces+1; i++)
            {
                int i2 = i+NumFaces;
                v0 = Vertex(Faces(i)(0));
                v1 = Vertex(Faces(i)(1));
                c = .5*(v0+v1);
                tau = v1-v0;
                l = sqrt(sum(tau*tau));
                taun = tau/l;
                CycleFace.Points1(i) = c;
                CycleFace.Points2(i) = taun;
                CycleFace.Vectors(i) = l;

                CycleVertex.Points1(i) = v0;
                if(IsNode(Faces(i)(0)))
                {
                    CycleVertex.Points2(i) = taun;
                    CycleVertex.Vectors(i) = -1.0;
                }
                else
                {
                    CycleVertex.Points2(i) = -taun;
                    CycleVertex.Vectors(i) = 1.0;
                    IsNode(Faces(i)(0)) = 1;
                }

                CycleVertex.Points1(i2) = v1;
                if(IsNode(Faces(i)(1)))
                {
                    CycleVertex.Points2(i2) = -taun;
                    CycleVertex.Vectors(i2) = -1.0;
                }
                else
                {
                    CycleVertex.Points2(i2) = taun;
                    CycleVertex.Vectors(i2) = 1.0;
                    IsNode(Faces(i)(1)) = 1;
                }
            }
        }

        ArrVectPoint TransDiff(ArrPoint &Vertex, Cycle &eCV, Cycle &eCE)
        {
            int NumVertices = Vertex.rows();
            // compute transpose of derivative of 'vertices to cycle' function
            Point v0, v1;
            VectPoint tau, etau;
            TYPE l;
            Cone taun;
            VectCone etaun;
            ArrVectPoint eV(Range(1,NumVertices));
            eV = 0.0;
            Array<bool,1> IsNode(Range(1,Vertex.rows()));
            IsNode = 0;
            for(int i=1; i<Faces.rows()+1; i++)
            {
                int i2 = i+NumFaces;
                v0 = Vertex(Faces(i)(0));
                v1 = Vertex(Faces(i)(1));
                tau = v1-v0;
                l = sqrt(sum(tau*tau));
                taun = tau/l;
                eV(Faces(i)(0)) += .5*eCE.Points1(i);
                eV(Faces(i)(1)) += .5*eCE.Points1(i);
                etaun = eCE.Points2(i);
                eV(Faces(i)(0)) += eCV.Points1(i);
                if(IsNode(Faces(i)(0)))
                    etaun += eCV.Points2(i);
                else
                {
                    etaun -= eCV.Points2(i);
                    IsNode(Faces(i)(0)) = 1;
                }
                eV(Faces(i)(1)) += eCV.Points1(i2);
                if(IsNode(Faces(i)(1)))
                    etaun -= eCV.Points2(i2);
                else
                {
                    etaun += eCV.Points2(i2);
                    IsNode(Faces(i)(1)) = 1;
                }
                etau = (etaun-sum(taun*etaun)*taun)/l + (eCE.Vectors(i)(0)/l)*tau;
                eV(Faces(i)(0)) -= etau;
                eV(Faces(i)(1)) += etau;
            }
            return eV;
        }

};

template < typename TYPE, int DIM >
class CurveCycleTarget : public Target<TYPE,DIM>
{
	using Target<TYPE,DIM>::Phi;
	using Target<TYPE,DIM>::RX;
	using Target<TYPE,DIM>::TargetWeight;
        
	typedef TinyVector<TYPE,DIM> Point;
        typedef Array<Point,1> ArrPoint;

        typedef TinyVector<TYPE,DIM> VectPoint;
        typedef Array<VectPoint,1> ArrVectPoint;

        typedef TinyVector<int,2> Face;
        typedef Array<Face,1> ArrFace;

        typedef TinyVector<TYPE,1> Weight;
        typedef Array<Weight,1> ArrWeight;

        typedef TinyVector<TYPE,DIM> Cone;
        typedef TinyVector<TYPE,DIM> VectCone;
        typedef Array<Cone,1> ArrCone;

        typedef DoubleVectDiracsMeasure<TYPE,Point,Cone,Weight> Cycle;
        typedef DoubleVectDiracsMeasure<TYPE,VectPoint,VectCone,Weight> VectCycle;

        Function<TYPE> *PointFun;
        TYPE orderCone, accCone, orderEdge, accEdge;

        CurveCycle<TYPE,DIM> CycleY, CyclePhiX;
        ArrPoint VertexY;
        ArrFace FaceY, FaceX;
        ArrWeight WeightX, WeightY;
        DoubleVectKerHilbert<TYPE,Point,Cone,Weight,VectPoint,VectCone> HFace, HVertex;
        static const TYPE LambdaFace=1.0, LambdaCone=1.0;
        double SqNormY;

        void Set_muY_and_SqNormY()
        {
            CycleY.Compute(VertexY);
	    SqNormY = LambdaFace * HFace.ScalProd(CycleY.CycleFace,CycleY.CycleFace) + LambdaCone * HVertex.ScalProd(CycleY.CycleVertex,CycleY.CycleVertex);
        }

    public:

        CurveCycleTarget(ifstream &f, int verbosemode = 1)
        {
            if(verbosemode) cout << "Reading Curve Cycle Target" << endl;
            string Tag;
            int debut, fin;
            for(int i=0; i<10; i++)
            {
                f >> Tag;
                if(!Tag.compare("Range="))
                {
                    f >> debut >> fin;
                    Target<TYPE,DIM>::RX = Range(debut,fin);
                }
                else if(!Tag.compare("Weight="))
                    f >> TargetWeight;
                else if(!Tag.compare("FunctionPoint="))
                    PointFun = ReadFunction<TYPE>(f);
                else if(!Tag.compare("OrderConeFunction="))
                    f >> orderCone;
                else if(!Tag.compare("AccConeFunction="))
                    f >> accCone;
                else if(!Tag.compare("OrderEdgeFunction="))
                    f >> orderEdge;
                else if(!Tag.compare("AccEdgeFunction="))
                    f >> accEdge;
                else if(!Tag.compare("VY="))
                    ReadArr(VertexY,f);
                else if(!Tag.compare("FY="))
                    ReadArr(FaceY,f);
                else if(!Tag.compare("FX="))
                    ReadArr(FaceX,f);
                else
                    cout << "error reading file." << endl;

            }

            VertexY.reindexSelf(1);
            FaceY.reindexSelf(1);
            FaceX.reindexSelf(1);
            Kernel<Point,Weight,VectPoint> *PointKer = new SqDistScalarKernel<TYPE,DIM,Weight>(PointFun);
            Kernel<Cone,Weight,VectCone> *EdgeKer, *ConeKer;
            switch(DIM)
            {
            	case 2 :
            		EdgeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new EdgeFunction<TYPE>(orderEdge,accEdge)));
            		ConeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new ConeFunction<TYPE>(orderCone,accCone)));
            		break;
            	case 3 :
            		EdgeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new EdgeFunction3D<TYPE>(orderEdge,accEdge)));
            		ConeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new ConeFunction3D<TYPE>(orderCone,accCone)));
            		//EdgeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new EdgeFunction<TYPE>(orderEdge,accEdge)));
            		//ConeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new ConeFunction<TYPE>(orderCone,accCone)));
            		break;
            	default :
            		cout << "error: CurveCycle implemented for dim 2 or 3 only" << endl;
            		assert(0);
            		break;
            }
            HFace.Ker = new MultKernel<TYPE,Point,Cone,Weight,VectPoint,VectCone>(PointKer,EdgeKer);
            HVertex.Ker = new MultKernel<TYPE,Point,Cone,Weight,VectPoint,VectCone>(PointKer,ConeKer);
            CycleY.Init(FaceY);
            CyclePhiX.Init(FaceX);
            Set_muY_and_SqNormY();
        }

        void Write(ofstream &f)
        {
            f << "CurveCycleTarget" << endl;
            f << "Range=" << endl << Target<TYPE,DIM>::RX(0) << " " << Target<TYPE,DIM>::RX(1)<< endl;
            f << "Weight=" << endl << TargetWeight << endl;
            f << "FunctionPoint=" << endl;
            PointFun->Write(f);
            f << "OrderConeFunction=" << endl << orderCone << endl;
            f << "AccConeFunction=" << endl << accCone << endl;
            f << "OrderEdgeFunction=" << endl << orderEdge << endl;
            f << "AccEdgeFunction=" << endl << accEdge << endl;
            f << "VY=" << endl;
            WriteArr(VertexY,f);
            f << "FX=" << endl;
            WriteArr(FaceX,f);
            f << "FY=" << endl;
            WriteArr(FaceY,f);
            f << "WX,WY=" << endl;
            WriteArr(WeightX,f);
            WriteArr(WeightY,f);
        }

        double Eval()
        {
            CyclePhiX.Compute(Phi);
            return LambdaFace * HFace.ScalProd(CyclePhiX.CycleFace,CyclePhiX.CycleFace-(TYPE)2.0*CycleY.CycleFace)
                   + LambdaCone * HVertex.ScalProd(CyclePhiX.CycleVertex,CyclePhiX.CycleVertex-(TYPE)2.0*CycleY.CycleVertex)
                   + SqNormY;
        }

        ArrVectPoint Gradient()
        {
            CyclePhiX.Compute(Phi);
            VectCycle eCycleFace(HFace.GradScalProd(CyclePhiX.CycleFace,CyclePhiX.CycleFace-CycleY.CycleFace));
            eCycleFace.Points1 *= (TYPE)(2.0*LambdaFace);
            eCycleFace.Points2 *= (TYPE)(2.0*LambdaFace);
            eCycleFace.Vectors *= (TYPE)(2.0*LambdaFace);
            VectCycle eCycleVertex(HVertex.GradScalProd(CyclePhiX.CycleVertex,CyclePhiX.CycleVertex-CycleY.CycleVertex));
            eCycleVertex.Points1 *= (TYPE)(2.0*LambdaCone);
            eCycleVertex.Points2 *= (TYPE)(2.0*LambdaCone);
            eCycleVertex.Vectors *= (TYPE)(2.0*LambdaCone);
SHOW(LambdaFace)  
SHOW(LambdaCone)          
SHOW(eCycleFace.Points1)          
SHOW(eCycleFace.Points2)          
SHOW(eCycleFace.Vectors)          
SHOW(eCycleVertex.Points1)          
SHOW(eCycleVertex.Points2)
SHOW(eCycleVertex.Vectors)
            ArrVectPoint Grad(CyclePhiX.TransDiff(Phi,eCycleVertex,eCycleFace));
            return Grad;
        }

};

template < typename TYPE, int DIM >
class CurveCycleScalarProducts
{
        typedef TinyVector<TYPE,DIM> Point;
        typedef Array<Point,1> ArrPoint;

        typedef TinyVector<TYPE,DIM> VectPoint;
        typedef Array<VectPoint,1> ArrVectPoint;

        typedef TinyVector<int,2> Face;
        typedef Array<Face,1> ArrFace;

        typedef TinyVector<TYPE,1> Weight;
        typedef Array<Weight,1> ArrWeight;

        typedef TinyVector<TYPE,DIM> Cone;
        typedef TinyVector<TYPE,DIM> VectCone;
        typedef Array<Cone,1> ArrCone;

        typedef DoubleVectDiracsMeasure<TYPE,Point,Cone,Weight> Cycle;
        typedef DoubleVectDiracsMeasure<TYPE,VectPoint,VectCone,Weight> VectCycle;

        Function<TYPE> *PointFun;
        TYPE orderCone, accCone, orderEdge, accEdge;

        Array<CurveCycle<TYPE,DIM>,1> AllCycle;
        Array<ArrPoint,1> AllVertex;
        Array<ArrFace,1> AllFace;
        DoubleVectKerHilbert<TYPE,Point,Cone,Weight,VectPoint,VectCone> HFace, HVertex;
        static const TYPE LambdaFace=1.0, LambdaCone=1.0;

        int N;
	Array<TYPE,2> TabScalProds;

        double ComputeScalProd(CurveCycle<TYPE,DIM>& CycleX, CurveCycle<TYPE,DIM>& CycleY)
        {
	    double ScalProd = LambdaFace * HFace.ScalProd(CycleX.CycleFace,CycleY.CycleFace) + LambdaCone * HVertex.ScalProd(CycleX.CycleVertex,CycleY.CycleVertex);
            return ScalProd;
        }

    public:

        CurveCycleScalarProducts(ifstream &f)
        {
            cout << "Reading Curve Cycle parameters" << endl;
            string Tag;
            int debut, fin;
            for(int i=0; i<6; i++)
            {
                f >> Tag;
                if(!Tag.compare("FunctionPoint="))
                    PointFun = ReadFunction<TYPE>(f);
                else if(!Tag.compare("OrderConeFunction="))
                    f >> orderCone;
                else if(!Tag.compare("AccConeFunction="))
                    f >> accCone;
                else if(!Tag.compare("OrderEdgeFunction="))
                    f >> orderEdge;
                else if(!Tag.compare("AccEdgeFunction="))
                    f >> accEdge;
		else if(!Tag.compare("N="))
		    f >> N;
                else
                    cout << "error reading file." << endl;
	    }
            TabScalProds.resize(Range(1,N),Range(1,N));
            AllCycle.resize(Range(1,N));
            AllVertex.resize(Range(1,N));
            AllFace.resize(Range(1,N));
            Kernel<Point,Weight,VectPoint> *PointKer = new SqDistScalarKernel<TYPE,DIM,Weight>(PointFun);
            Kernel<Cone,Weight,VectCone> *EdgeKer, *ConeKer;
            switch(DIM)
            {
            	case 2 :
            		EdgeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new EdgeFunction<TYPE>(orderEdge,accEdge)));
            		ConeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new ConeFunction<TYPE>(orderCone,accCone)));
            		break;
            	case 3 :
            		EdgeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new EdgeFunction3D<TYPE>(orderEdge,accEdge)));
            		ConeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new ConeFunction3D<TYPE>(orderCone,accCone)));
            		//EdgeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new EdgeFunction<TYPE>(orderEdge,accEdge)));
            		//ConeKer = new AngleKernel<TYPE,DIM,Weight>(new LookUpFunction1D<TYPE>(new ConeFunction<TYPE>(orderCone,accCone)));
            		break;
            	default :
            		cout << "error: CurveCycle implemented for dim 2 or 3 only" << endl;
            		assert(0);
            		break;
            }
            HFace.Ker = new MultKernel<TYPE,Point,Cone,Weight,VectPoint,VectCone>(PointKer,EdgeKer);
            HVertex.Ker = new MultKernel<TYPE,Point,Cone,Weight,VectPoint,VectCone>(PointKer,ConeKer);

            cout << "Starting reading Curve Cycle data" << endl;    
	    for(int p=1; p<N+1; p++)
            {
                    cout << "p=" << p << endl;
            	    for(int i=0; i<2; i++)
                    {
                        f >> Tag;
                        if(!Tag.compare("V="))
                            ReadArr(AllVertex(p),f);
                        else if(!Tag.compare("F="))
                            ReadArr(AllFace(p),f);
                        else
                            cout << "error reading file." << endl;
                    }
		    AllCycle(p).Init(AllFace(p));
	            AllCycle(p).Compute(AllVertex(p));
	    }
        }

        void ComputeAllScalarProds()
        {
            cout << "Starting computing scalar products" << endl;
            for(int p=1; p<N+1; p++)
                for(int q=1; q<N+1; q++)
                {
                    cout << "p=" << p << ", q=" << q << endl;
                    TabScalProds(p,q) = ComputeScalProd(AllCycle(p),AllCycle(q));
                }
        }

        void Write(ofstream &f)
        {
            f << "CurveCycleScalarProducts" << endl;
            f << TabScalProds;
        }

};

#endif // CURVECYCLE
