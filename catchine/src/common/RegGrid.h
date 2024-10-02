#ifndef REGGRID
#define REGGRID

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <fftw3.h>


using namespace blitz;

template < typename TYPE, int DIMPOINT >
class RegGrid
{
        typedef TinyVector<TYPE,DIMPOINT> VectPoint;
        typedef Array<VectPoint,1> ArrVectPoint;
        typedef TinyVector<int,DIMPOINT> Param;
        typedef Array<TYPE,DIMPOINT> ImageND;

        VectPoint Origine;
        Param Long;
        TYPE Pas;
        ImageND Image;
        // char *fftw_wisdom;

    public:
        RegGrid(VectPoint origine, Param leng, TYPE pas, ImageND image) : Origine(origine), Long(leng), Pas(pas), Image(image) {}

        RegGrid(TYPE pas)
        {
            Pas = pas;
            Image = NULL;
            // fftw_wisdom = NULL;
            // int t = fftw_init_threads();
            // cout << "t = " << t << endl;
        }

        TinyVector<TYPE,DIMPOINT> getLength()
        {
            return Long;
        }

        TYPE* getImage()
        {
            return Image.data();
        }

        void changeOrigin(VectPoint origine)
        {
            Origine = origine;
        }

        void resizeGrid(int longX,int longY,int longZ)
        {
            Long(0) = longX;
            Long(1) = longY;
            Long(2) = longZ;

            Image.resize(longX,longY,longZ);
        }

        bool isSameSize(int LongX, int LongY, int LongZ) const
        {
            bool p = (Long(0)==LongX)&&(Long(1)==LongY)&&(Long(2)==LongZ);
            return p;
        }

        void setGridFFTFct(int longX, int longY, int longZ, TYPE sigma)
        {
            Long(0) = longX;
            Long(1) = longY;
            Long(2) = longZ;
            TYPE facteur = pow2(Pas / sigma);

            Image.resize(longX,longY,longZ);
            TYPE *ImagePt = Image.data();

            firstIndex i;
            secondIndex j;
            thirdIndex k;
            Image = exp( -( pow2(i) + pow2(j) + pow2(k) )*facteur );

            fftw_plan forward = fftw_plan_r2r_3d(longX,longY,longZ, ImagePt, ImagePt, FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
            fftw_execute(forward);

            fftw_destroy_plan(forward);
            return;
        }

        void setGridFct(int longX, int longY, int longZ, TYPE sigma)
        {
            Long(0) = longX;
            Long(1) = longY;
            Long(2) = longZ;
            TYPE facteur = pow2(Pas / sigma);

            Image.resize(longX,longY,longZ);
            TYPE *ImagePt = Image.data();

            firstIndex i;
            secondIndex j;
            thirdIndex k;
            Image = exp( -( pow2(i) + pow2(j) + pow2(k) )*facteur );

            return;
        }

        void proj3DLin(ArrVectPoint &x, TYPE *alpha, bool flag)
        {
            // flag = 0 for interpolation, flag = 1 for extrapolation

            TYPE* xcont = contig(x);
            TYPE* ImagePt = Image.data();
            int Nx = x.extent(firstDim);
            int c0X, c0Y,c0Z;
            int nx = Long(0), ny = Long(1), nz = Long(2), N = ny*nz;
            TYPE deltaX, deltaY, deltaZ, rho000, rho100, rho010, rho001, rho110, rho011, rho101, rho111, val;
            TYPE Pas3 = Pas*Pas*Pas;

            if (!flag)
            {
                Image = 0.0;
            }

            for (int k=0; k<Nx; k++)
            {
                c0X = (int) ((xcont[3*k]   - Origine(0))/Pas);
                c0Y = (int) ((xcont[1+3*k] - Origine(1))/Pas);
                c0Z = (int) ((xcont[2+3*k] - Origine(2))/Pas);

                if (( ( (c0X+1)*N+(c0Y+1)*nz+c0Z+1 ) >= (N*nx) )||( c0X*N+c0Y*nz+c0Z<0 ))
                {
                    cerr << "proj3DLin: ERREUR PROJECTION flag = " << flag << endl;
                    exit(1);
                }

                deltaX = xcont[3*k]   - (Origine(0) + c0X*Pas);
                deltaY = xcont[1+3*k] - (Origine(1) + c0Y*Pas);
                deltaZ = xcont[2+3*k] - (Origine(2) + c0Z*Pas);

                rho000 = (Pas-deltaX) * (Pas-deltaY)  * (Pas-deltaZ);
                rho100 = deltaX       * (Pas-deltaY)  * (Pas-deltaZ);
                rho010 = (Pas-deltaX) * deltaY        * (Pas-deltaZ);
                rho001 = (Pas-deltaX) * (Pas-deltaY)  * deltaZ      ;
                rho110 = deltaX       * deltaY        * (Pas-deltaZ);
                rho011 = (Pas-deltaX) * deltaY        * deltaZ      ;
                rho101 = deltaX       * (Pas-deltaY)  * deltaZ      ;
                rho111 = deltaX       * deltaY        * deltaZ      ;

                if (flag)
                {
                    alpha[k] = ( ImagePt[c0X*N     + c0Y*nz     + c0Z    ]*rho000+
                                 ImagePt[c0X*N     + c0Y*nz     + c0Z+1  ]*rho001+
                                 ImagePt[c0X*N     + (c0Y+1)*nz + c0Z    ]*rho010+
                                 ImagePt[c0X*N     + (c0Y+1)*nz + c0Z+1  ]*rho011+
                                 ImagePt[(c0X+1)*N + c0Y*nz     + c0Z    ]*rho100+
                                 ImagePt[(c0X+1)*N + c0Y*nz     + c0Z+1  ]*rho101+
                                 ImagePt[(c0X+1)*N + (c0Y+1)*nz + c0Z    ]*rho110+
                                 ImagePt[(c0X+1)*N + (c0Y+1)*nz + c0Z+1  ]*rho111 )/Pas3;

                }
                else
                {
                    val = alpha[k]/Pas3;
                    ImagePt[c0X*N     + c0Y*nz     + c0Z    ] += rho000*val;
                    ImagePt[c0X*N     + c0Y*nz     + c0Z+1  ] += rho001*val;
                    ImagePt[c0X*N     + (c0Y+1)*nz + c0Z    ] += rho010*val;
                    ImagePt[c0X*N     + (c0Y+1)*nz + c0Z+1  ] += rho011*val;
                    ImagePt[(c0X+1)*N + c0Y*nz     + c0Z    ] += rho100*val;
                    ImagePt[(c0X+1)*N + c0Y*nz     + c0Z+1  ] += rho101*val;
                    ImagePt[(c0X+1)*N + (c0Y+1)*nz + c0Z    ] += rho110*val;
                    ImagePt[(c0X+1)*N + (c0Y+1)*nz + c0Z+1  ] += rho111*val;
                }
            }
            return;
        }



        void ConvImage(RegGrid &fftKer)
        {
            TinyVector<int,DIMPOINT> LongKer = fftKer.getLength();
            if ( ( LongKer(0)!=(Long(0)/2+1) )||( LongKer(1)!=(Long(1)/2+1) )||( LongKer(2)!=(Long(2)/2) ) )
            {
                cerr << "ConvImage: dimensions mismatch: " << "Long ker = " << LongKer << endl << "Long image = " << Long << endl;
                cerr << "Warning: Padded input image needed for in-place FFT" << endl;
                exit(1);
            }

            int LongX = Long(0), LongY = Long(1), LongZ = Long(2) - 2;
            int LongX2 = LongKer(0), LongY2 = LongKer(1), LongZ2 = LongKer(2);
            int N = LongX*LongY*LongZ;

            TYPE *imagePt = Image.data(), *fftKerPt = fftKer.getImage();
            TYPE val;
            int ii,jj,kk;
            fftw_plan forward, backward;

            // if (fftw_wisdom != NULL) {fftw_import_wisdom_from_string(fftw_wisdom);}
            // fftw_plan_with_nthreads(2);
            forward = fftw_plan_dft_r2c_3d(LongX,LongY,LongZ, imagePt, (fftw_complex *) imagePt, FFTW_ESTIMATE);
            fftw_execute(forward);


            // fftw_wisdom = fftw_export_wisdom_to_string();
            // fftw_import_wisdom_from_string(fftw_wisdom);


            for (int i=0; i<LongX; i++)
            {
                ii = (i<LongX2)?i:(LongX-i);
                for (int j=0; j<LongY; j++)
                {
                    jj = (j<LongY2)?j:(LongY-j);
                    for (int k=0; k<LongZ2; k++)
                    {
                        val = fftKerPt[ii*LongZ2*LongY2 + jj*LongZ2 + k]/N;
                        imagePt[i*LongY*(LongZ+2) + j*(LongZ+2) + 2*k] *= val; // real part
                        imagePt[i*LongY*(LongZ+2) + j*(LongZ+2) + 2*k+1] *= val; // imaginary part
                    }
                }
            }

            // fftw_plan_with_nthreads(2);
            backward = fftw_plan_dft_c2r_3d(LongX,LongY,LongZ, (fftw_complex *) imagePt, imagePt, FFTW_ESTIMATE);
            fftw_execute(backward);

            // fftw_wisdom = fftw_export_wisdom_to_string();

            // fftw_cleanup_threads();
            fftw_destroy_plan(forward);
            fftw_destroy_plan(backward);

            return;
        }


        void projConvProj(ArrVectPoint &xIn, Array<TYPE,1> &CoeffIn, ArrVectPoint &xOut, Array<TYPE,1> &CoeffOut, RegGrid &fftKer)
        {
            if ((xIn.extent(firstDim)!=CoeffIn.extent(firstDim))||(xOut.extent(firstDim)!=CoeffOut.extent(firstDim)))
            {
                cerr << "dimensions mismatched in projConvProj : " << xIn.extent(firstDim) << " " << CoeffIn.extent(firstDim) << " " << xOut.extent(firstDim) << " " << CoeffOut.extent(firstDim) << endl;
                exit(1);
            }
            proj3DLin(xIn, (TYPE *) CoeffIn.data(), 0);
            ConvImage(fftKer);
            proj3DLin(xOut, (TYPE *) CoeffOut.data(), 1);
        }




};

#endif //REGGRID