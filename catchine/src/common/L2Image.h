#ifndef L2IMAGE
#define L2IMAGE

#include "Target.h"

template < typename TYPE >
class L2Image : public Target<TYPE,3>
{

		using Target<TYPE,3>::Phi;
		using Target<TYPE,3>::RX;
		using Target<TYPE,3>::TargetWeight;
		
        typedef TinyVector<TYPE,3> Vect;
        typedef Array<Vect,1> ArrVect;
        typedef Array<Vect,3> GridVect;
        typedef TinyVector<int,3> Index3;
        typedef Array<TYPE,1> ArrIntens;
        typedef Array<TYPE,3> GridIntens;
        struct RegGrid
        {
            Vect Base;
            Vect VoxSize;
        };

        Index3 indexones;
        GridVect SourceCenters;
        GridIntens ImSource, ImTarget, HitCoef;
        RegGrid RegGridTarget;

        void Init(GridIntens &imsource, GridIntens &imtarget, Vect basetarget, Vect voxsizetarget, Range rx, double weight)
        {
            indexones = 1,1,1;
            ImSource.reference(imsource);
            ImSource.reindexSelf(indexones);
            SourceCenters.resize(ImSource.shape());
            ImTarget.reference(imtarget);
            ImTarget.reindexSelf(indexones);
            HitCoef.resize(ImTarget.shape());
            HitCoef.reindexSelf(indexones);
            RegGridTarget.Base = basetarget;
            RegGridTarget.VoxSize = voxsizetarget;
            RX = rx;
            TargetWeight = weight;
        }

        GridVect CompGridCenters(GridVect &V)
        {
            GridVect C(V.shape()-1);
            C.reindexSelf(indexones);
            for(int i=1; i<V.rows(); i++)
                for(int j=1; j<V.columns(); j++)
                    for(int k=1; k<V.depth(); k++)
                        C(i,j,k) = 0.125 * (V(i,j,k)+V(i,j,k+1)+V(i,j+1,k)+V(i,j+1,k+1)+V(i+1,j,k)+V(i+1,j,k+1)+V(i+1,j+1,k)+V(i+1,j+1,k+1));
            return C;
        }

        void ImageInterp(GridIntens &GridIntensInterp, GridVect &GridPosInterp, GridIntens &Image, RegGrid RegGridImage)
        {
            ArrVect PosInterp(GridPosInterp.data(),shape(GridPosInterp.size()),neverDeleteData);
            PosInterp.reindexSelf(1);
            ArrIntens IntensInterp(GridIntensInterp.data(),PosInterp.shape(),neverDeleteData);
            IntensInterp.reindexSelf(1);
            ImageInterp(IntensInterp, PosInterp, Image, RegGridImage);
        }

        void ImageInterp(ArrIntens &IntensInterp, ArrVect &PosInterp, GridIntens &Image, RegGrid RegGridImage)
        {
            int dohit;
            dohit = HitCoef.rows();
            IntensInterp = 0.0;
            Array<TinyVector<int,8>,3> HitFlags;
            if(dohit)
            {
                HitFlags.resize(Image.shape());
                HitFlags = 1;
                HitFlags.reindexSelf(indexones);
            }
            Index3 otf;
            otf = 1,2,4;
            for(int i=1; i<PosInterp.rows()+1; i++)
            {
                Vect f = (PosInterp(i)-RegGridImage.Base)/RegGridImage.VoxSize + 0.5;
                Index3 index = floor(f);
                f = f-index;
                Index3 e;
                int be;
                for(e(0)=0; e(0)<2; e(0)++)
                    for(e(1)=0; e(1)<2; e(1)++)
                        for(e(2)=0; e(2)<2; e(2)++)
                        {
                            Index3 indpe = index+e;
                            Vect ef = e*f+(1-e)*(1-f);
                            if(indpe(0)>0 && indpe(1)>0 && indpe(2)>0 && indpe(0)<Image.rows()+1 && indpe(1)<Image.columns()+1 && indpe(2)<Image.depth()+1)
                            {
                                IntensInterp(i) += product(ef) * Image(indpe);
                                if(dohit)
                                {
                                    be = sum(otf*e);
                                    HitFlags(indpe)(be) = 0;
                                }
                            }
                        }
            }
            if(dohit)
                for(int i=1; i<HitCoef.rows(); i++)
                    for(int j=1; j<HitCoef.columns(); j++)
                        for(int k=1; k<HitCoef.depth(); k++)
                            HitCoef(i, j, k) = sum(HitFlags(i, j, k));
        }

        void DiffImageInterp(GridVect &GridPosInterp, GridIntens &Image, RegGrid RegGridImage, GridIntens &CoefIn, GridVect &GradOut)
        {
            Index3 i;
            for(i(0)=1; i(0)<GridPosInterp.rows()+1; i(0)++)
                for(i(1)=1; i(1)<GridPosInterp.columns()+1; i(1)++)
                    for(i(2)=1; i(2)<GridPosInterp.depth()+1; i(2)++)
                    {
                        TYPE coef = 0.125 * CoefIn(i);
                        Vect f = (GridPosInterp(i)-RegGridImage.Base)/RegGridImage.VoxSize + 0.5;
                        Index3 index = floor(f);
                        f = f-index;
                        Index3 e;
                        for(e(0)=0; e(0)<2; e(0)++)
                            for(e(1)=0; e(1)<2; e(1)++)
                                for(e(2)=0; e(2)<2; e(2)++)
                                {
                                    Index3 indpe = index+e;
                                    TYPE coef2 = 0;
                                    if(indpe(0)>0 && indpe(1)>0 && indpe(2)>0 && indpe(0)<Image.rows()+1 && indpe(1)<Image.columns()+1 && indpe(2)<Image.depth()+1)
                                        coef2 = coef * Image(indpe);
                                    Vect ef = e*f+(1-e)*(1-f);
                                    Vect grad;
                                    grad(0) = (2*e(0)-1)*ef(1)*ef(2);
                                    grad(1) = (2*e(1)-1)*ef(2)*ef(0);
                                    grad(2) = (2*e(2)-1)*ef(0)*ef(1);
                                    grad *= (coef2/RegGridImage.VoxSize);
                                    Index3 ep;
                                    for(ep(0)=0; ep(0)<2; ep(0)++)
                                        for(ep(1)=0; ep(1)<2; ep(1)++)
                                            for(ep(2)=0; ep(2)<2; ep(2)++)
                                            {
                                                Index3 ipep = i+ep;
                                                GradOut(ipep) += grad;
                                            }
                                }

                    }

        }

        GridIntens CompGridVolume(GridVect &V)
        {
            GridIntens Vol(V.shape()-1);
            Vol.reindexSelf(indexones);
            for(int i=1; i<V.rows(); i++)
                for(int j=1; j<V.columns(); j++)
                    for(int k=1; k<V.depth(); k++)
                    {
                        Vect u = V(i+1,j,k)-V(i,j,k);
                        Vect v = V(i,j+1,k)-V(i,j,k);
                        Vect w = V(i,j,k+1)-V(i,j,k);
                        Vol(i,j,k) = det(u,v,w);
                    }
            return Vol;
        }


        void DiffGridVolume(GridVect &V, GridIntens &CoefIn, GridVect &GradOut)
        {
            GridIntens Vol(V.shape()-1);
            Vol.reindexSelf(indexones);
            for(int i=1; i<V.rows(); i++)
                for(int j=1; j<V.columns(); j++)
                    for(int k=1; k<V.depth(); k++)
                    {
                        Vect u = V(i+1,j,k)-V(i,j,k);
                        Vect v = V(i,j+1,k)-V(i,j,k);
                        Vect w = V(i,j,k+1)-V(i,j,k);
                        Vect cvw = CoefIn(i,j,k) * cross(v,w);
                        GradOut(i+1,j,k) += cvw;
                        Vect cwu = CoefIn(i,j,k) * cross(w,u);
                        GradOut(i,j+1,k) += cwu;
                        Vect cuv = CoefIn(i,j,k) * cross(u,v);
                        GradOut(i,j,k+1) += cuv;
                        GradOut(i,j,k) -= cvw + cwu + cuv;
                    }
        }


    public:

        L2Image(ifstream &f, int verbosemode = 1)
        {
            if(verbosemode) cout << "Reading L2 Image Target" << endl;
            GridIntens imsource, imtarget;
            Vect basetarget, voxsizetarget;
            Range rx;
            double weight;

            string Tag;
            for(int i=0; i<6; i++)
            {
                f >> Tag;
                if(!Tag.compare("Range="))
                {
                    int debut, fin;
                    f >> debut >> fin;
                    rx = Range(debut,fin);
                }
                else if(!Tag.compare("Weight="))
                    f >> weight;
                else if(!Tag.compare("SourceImage="))
                    ReadArr(imsource,f);
                else if(!Tag.compare("TargetImage="))
                    ReadArr(imtarget,f);
                else if(!Tag.compare("TargetGridBase="))
                    f >> basetarget(0) >> basetarget(1) >> basetarget(2);
                else if(!Tag.compare("TargetGridVoxSize="))
                    f >> voxsizetarget(0) >> voxsizetarget(1) >> voxsizetarget(2);
                else
				{
                    cout << "error reading file." << endl;
					throw -1;
				}
            }
            Init(imsource, imtarget, basetarget, voxsizetarget, rx, weight);
			if(verbosemode) cout << "done." << endl;
        }

        void Write(ofstream &f)
        {
            f << "L2Image" << endl;
            f << "Range=" << endl << RX(0) << " " << RX(1)<< endl;
            f << "Weight=" << endl << TargetWeight << endl;
            f << "SourceImage=" << endl;
            WriteArr(ImSource,f);
            f << "TargetImage=" << endl;
            WriteArr(ImTarget,f);
            f << "TargetGridBase=" << endl << RegGridTarget.Base(0) << " " << RegGridTarget.Base(1) << " " << RegGridTarget.Base(2) << " " << endl;
            f << "TargetGridVoxSize=" << endl << RegGridTarget.VoxSize(0) << " " << RegGridTarget.VoxSize(1) << " " << RegGridTarget.VoxSize(2) << endl;
        }

        L2Image(GridIntens &imsource, GridIntens &imtarget, Vect basetarget, Vect voxsizetarget, Range rx, double weight)
        {
            Init(imsource, imtarget, basetarget, voxsizetarget, rx, weight);
        }

        double Eval()
        {
            GridVect PhiGrid(Phi.data(),ImSource.shape()+1,neverDeleteData);
            PhiGrid.reindexSelf(indexones);
            SourceCenters = CompGridCenters(PhiGrid);
            GridIntens ImTargetInterp(SourceCenters.shape());
            ImTargetInterp.reindexSelf(indexones);
            ImageInterp(ImTargetInterp,SourceCenters,ImTarget,RegGridTarget);
            GridIntens VolGridSource = CompGridVolume(PhiGrid);
            return sum(sqr(ImSource-ImTargetInterp)*VolGridSource) + sum(sqr(ImTarget)*HitCoef*(product(RegGridTarget.VoxSize)/8));
        }

        ArrVect Gradient()
        {
            GridVect PhiGrid(Phi.data(),ImSource.shape()+1,neverDeleteData);
            PhiGrid.reindexSelf(indexones);
            SourceCenters = CompGridCenters(PhiGrid);
            GridIntens ImTargetInterp(SourceCenters.shape());
            ImTargetInterp.reindexSelf(indexones);
            ImageInterp(ImTargetInterp,SourceCenters,ImTarget,RegGridTarget);
            GridIntens VolGridSource = CompGridVolume(PhiGrid);
            ArrVect grad(Range(1,Phi.rows()));
            grad = 0.0;
            GridVect gradGrid(grad.data(),ImSource.shape()+1,neverDeleteData);
            gradGrid.reindexSelf(indexones);
            GridIntens CoefIn(ImSource.shape());
            CoefIn.reindexSelf(indexones);
            CoefIn = -2.0*(ImSource-ImTargetInterp)*VolGridSource;
            DiffImageInterp(SourceCenters,ImTarget,RegGridTarget,CoefIn,gradGrid);
            CoefIn = sqr(ImSource-ImTargetInterp);
            DiffGridVolume(PhiGrid,CoefIn,gradGrid);
            return grad;
        }

};

#endif // L2IMAGE

