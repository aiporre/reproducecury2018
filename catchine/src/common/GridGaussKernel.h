#ifndef GRIDGAUSSKERNEL
#define GRIDGAUSSKERNEL

#include "Kernel.h"
#include "Utils.h"
#include "RegGrid.h"
#include <omp.h>

template < typename TYPE, int DIMPOINT, int DIMVECT >
class GridGaussKernel : public Kernel<TYPE,DIMPOINT,DIMVECT>
{
        typedef TinyVector<TYPE,DIMVECT> Vect;
        typedef Array<Vect,1> ArrVect;
        typedef TinyVector<TYPE,DIMPOINT> VectPoint;
        typedef Array<VectPoint,1> ArrVectPoint;

        TYPE Sigma, ooSigma2, Ratio, Margin, Step;
        RegGrid<TYPE,DIMPOINT> fftKer;
        RegGrid<TYPE,DIMPOINT> Image;


        ArrVectPoint boundingBox(ArrVectPoint &x, ArrVectPoint &y)
        {
            ArrVectPoint  minmax(2);
            ArrVectPoint minmaxX = minMaxPerDim<TYPE,DIMPOINT>(x);
            ArrVectPoint minmaxY = minMaxPerDim<TYPE,DIMPOINT>(y);
            for (int d=0; d<DIMPOINT; d++)
            {
                minmax(firstDim)(d) = min( minmaxX(firstDim)(d), minmaxY(firstDim)(d) );
                minmax(secondDim)(d) = max( minmaxX(secondDim)(d), minmaxY(secondDim)(d) );
            }
            return minmax;
        }

        void adjustGridSize(ArrVectPoint &minmax)
        {
            TinyVector<int,DIMPOINT> leng;
            VectPoint length = (minmax(1) - minmax(0))/Step + Margin;
            for (int d=0; d<DIMPOINT; d++)
            {
                leng(d) = (length(d)<=16.0)*16 + (length(d)>16.0)*((length(d)<=32.0)*32 + (length(d)>32.0)*((length(d)<=64.0)*64 + (length(d)>64.0)*4*( (int) (length(d)/4) + 1)));
            }
            if ( !Image.isSameSize(leng(0),leng(1),leng(2)+2) )
            {
                cout << "mise a jour du noyau : " << leng(0) << " " << leng(1) << " " << leng(2) << endl;
                fftKer.setGridFFTFct(leng(0)/2 +1, leng(1)/2+1, leng(2)/2+1, Sigma);
                Image.resizeGrid(leng(0),leng(1),leng(2)+2); // padded grid for further fft
            }
            fftKer.changeOrigin( minmax(0) - ( leng*Step - minmax(1) + minmax(0) )/2 );
            Image.changeOrigin( minmax(0) - ( leng*Step - minmax(1) + minmax(0) )/2 );
        }


    public:

        void Write(ofstream &f)
        {
            f << "GridGauss,sigma,ratio=" << endl << Sigma << " " << Ratio <<  endl;
        }

        GridGaussKernel(TYPE sigma, TYPE ratio) : Sigma(sigma), ooSigma2(1.0/pow2(sigma)), Ratio(ratio), fftKer(ratio * sigma), Image(ratio * sigma), Step(ratio*sigma), Margin(10.0/Ratio)
        {
            if (DIMPOINT != 3)
            {
                cerr << "GridGaussKernel not implemented for dimension different from 3!" << endl;
                exit(1);
            }
        }

        ArrVect ConvKer(ArrVectPoint &x, ArrVectPoint &y, ArrVect &alpha)
        {
            int nx = x.extent(firstDim);
            int ny = y.extent(firstDim);

            ArrVectPoint minmax = boundingBox(x,y);
            adjustGridSize(minmax);

            ArrVect gamma(ny);
            Array<TYPE,1> CoeffX(nx), CoeffY(ny);

            for (int d=0; d<DIMVECT; d++)
            {
                for (int i=0; i<nx; i++)
                    CoeffX(i) = alpha(alpha.base(0) + i)(d);
                Image.projConvProj(x, CoeffX, y, CoeffY, fftKer);
                for (int j=0; j<ny; j++)
                    gamma(gamma.base(0) + j)(d) = CoeffY(j);
            }
            // cout << "GridGaussKernel: gamma = " << gamma << endl;
            return gamma;
        }



        ArrVectPoint ConvGradKer(ArrVect &gamma, ArrVectPoint &x, ArrVectPoint &y, ArrVect &alpha)
        {
            int nx = x.extent(firstDim);
            int ny = y.extent(firstDim);
            ArrVectPoint beta(Range(1,nx));

            ArrVectPoint minmax = boundingBox(x,y);
            adjustGridSize(minmax);

            int dimc = (DIMVECT*(1+DIMPOINT));
            Array<TYPE,1> Coeffs(ny*dimc);
            for(int i=0; i<ny; i++)
            {
                for(int d=0; d<DIMVECT; d++)
                {
                    Coeffs(i*dimc+d) = gamma(i+1)(d);
                    for(int dp=0; dp<DIMPOINT; dp++)
                        Coeffs(i*dimc+DIMVECT+DIMPOINT*d+dp) = gamma(i+1)(d)*y(i+1)(dp);
                }
            }

            Array<TYPE,1> betacont(nx*dimc);
            #pragma omp parallel for
            for (int k=0; k<dimc; k++)
            {
                Array<TYPE,1> CoeffX(nx), CoeffY(ny);

                // Array<TYPE,1> CoeffY = Coeffs(Range(k,k+(ny-1)*dimc,dimc));
                // Array<TYPE,1> CoeffX = betacont(Range(k,k+(nx-1)*dimc,dimc));
                for (int i=0; i<ny; i++)
                    CoeffY(i) = Coeffs(k + i*dimc);
                Image.projConvProj(y, CoeffY, x, CoeffX, fftKer);
                for (int i=0; i<nx; i++)
                    betacont(k+i*dimc) = CoeffX(i);
            }

            int base1, base2;
            for(int i=1; i<nx+1; i++)
            {
                base1 = (i-1)*dimc;
                base2 = base1+DIMVECT;
                TYPE scalprod = 0;
                for(int dv=0; dv<DIMVECT; dv++)
                    scalprod += alpha(i)(dv)*betacont(base1+dv);
                for(int dp=0; dp<DIMPOINT; dp++)
                {
                    beta(i)(dp) = - scalprod * x(i)(dp);
                    for(int dv=0; dv<DIMVECT; dv++)
                        beta(i)(dp) += betacont(base2+dv*DIMPOINT+dp)*alpha(i)(dv);
                }
                beta(i) *= 2*ooSigma2;
            }
            return beta;
        }


        ArrVectPoint ConvSpecDiffTKer(ArrVectPoint &x, ArrVect &alpha, ArrVect &eta)
        {
            int nx = x.extent(firstDim);
            ArrVectPoint beta(Range(1,nx));

            ArrVectPoint minmax = minMaxPerDim(x);
            adjustGridSize(minmax);

            int dimc = 2*(DIMVECT*(1+DIMPOINT));
            Array<TYPE,1> Coeffs(nx*dimc);
            for(int i=0; i<nx; i++)
            {
                for(int d=0; d<DIMVECT; d++)
                {
                    int base = i*dimc;
                    Coeffs(base+d) = eta(i+1)(d);
                    base += DIMVECT;
                    for(int dp=0; dp<DIMPOINT; dp++)
                        Coeffs(base+DIMPOINT*d+dp) = eta(i+1)(d)*x(i+1)(dp);
                    base = i*dimc + (DIMPOINT+1)*DIMVECT;
                    Coeffs(base+d) = alpha(i+1)(d);
                    base += DIMVECT;
                    for(int dp=0; dp<DIMPOINT; dp++)
                        Coeffs(base+DIMPOINT*d+dp) = alpha(i+1)(d)*x(i+1)(dp);
                }
            }

            Array<TYPE,1> betacont(nx*dimc);
            #pragma omp parallel for
            for (int k=0; k<dimc; k++)
            {
                Array<TYPE,1> CoeffX(nx);
                for (int i=0; i<nx; i++)
                    CoeffX(i) = Coeffs(k + i*dimc);
                Image.projConvProj(x, CoeffX, x, CoeffX, fftKer);
                for (int i=0; i<nx; i++)
                    betacont(k + i*dimc) = CoeffX(i);
            }

            for(int i=1; i<nx+1; i++)
            {
                int base = (i-1)*dimc;
                TYPE scalprod = 0;
                for(int dv=0; dv<DIMVECT; dv++)
                    scalprod += alpha(i)(dv)*betacont(base+dv);
                base += DIMVECT;
                for(int dp=0; dp<DIMPOINT; dp++)
                {
                    beta(i)(dp) = 0;
                    for(int dv=0; dv<DIMVECT; dv++)
                        beta(i)(dp) += betacont(base+dv*DIMPOINT+dp)*alpha(i)(dv);
                }
                base = (i-1)*dimc + (DIMPOINT+1)*DIMVECT;
                for(int dv=0; dv<DIMVECT; dv++)
                    scalprod += eta(i)(dv)*betacont(base+dv);
                base += DIMVECT;
                for(int dp=0; dp<DIMPOINT; dp++)
                {
                    beta(i)(dp) -= scalprod * x(i)(dp);
                    for(int dv=0; dv<DIMVECT; dv++)
                        beta(i)(dp) += betacont(base+dv*DIMPOINT+dp)*eta(i)(dv);
                }

                beta(i) *= 2*ooSigma2;
            }
            return beta;
        }



        // ArrVectPoint ConvGradKer(ArrVect &gamma, ArrVectPoint &x, ArrVectPoint &y, ArrVect &alpha)
        // {
        //   int nx = x.extent(firstDim);
        //   int ny = y.extent(firstDim);
        //   ArrVectPoint beta(nx);
        //
        //   TYPE ooSigma2 = 1.0/pow2(Sigma);
        //
        //   ArrVectPoint minmax = boundingBox(x,y);
        //   adjustGridSize(minmax);
        //
        //   Array<TYPE,1> CoeffX(nx), CoeffY(ny), Aux1(nx), Aux2(nx);
        //
        //   Aux1 = 0.0;
        //   for (int p=0; p<DIMVECT; p++)
        //   {
        //     for (int j=0; j<ny; j++)
        //       CoeffY(j) = gamma(gamma.base(0) + j)(p);
        //     Image.projConvProj(y, CoeffY, x, CoeffX, fftKer);
        //     for (int i=0; i<nx; i++)
        //       Aux1(i) += CoeffX(i)*alpha(alpha.base(0) + i)(p);
        //   }
        //
        //   for (int dp=0; dp<DIMPOINT; dp++)
        //   {
        //     for (int i=0; i<nx; i++) {beta(beta.base(0)+i)(dp) = -Aux1(i)*x(x.base(0) + i)(dp);}
        //     Aux2 = 0.0;
        //     for (int dv=0; dv<DIMVECT; dv++)
        //     {
        //       for (int j=0; j<ny; j++)
        //         CoeffY(j) = gamma(gamma.base(0) + j)(dv)*y(y.base(0) + j)(dp);
        //       Image.projConvProj(y, CoeffY, x, CoeffX, fftKer);
        //       for (int i=0; i<nx; i++)
        //         Aux2(i) += CoeffX(i)*alpha(alpha.base(0) + i)(dv);
        //     }
        //     for (int i=0; i<nx; i++) {beta(beta.base(0)+i)(dp) += Aux2(i);}
        //   }
        //   beta *= 2*ooSigma2;
        //
        //   return beta;
        //
        // }


        //  // Ce code est vraisemblablement faux...
        // ArrVectPoint ConvSpecDiffTKer(ArrVectPoint &x, ArrVect &alpha, ArrVect &eta)
        // {
        //   int nx = x.extent(firstDim);
        //   ArrVectPoint beta(nx);
        //
        //   TYPE ooSigma2 = 1.0/pow2(Sigma);
        //
        //   ArrVectPoint minmax = minMaxPerDim(x);
        //   adjustGridSize(minmax);
        //
        //   Array<TYPE,1> CoeffY(nx), Aux1(nx), Aux2(nx);
        //
        //   Aux1 = 0.0;
        //   for (int p=0; p<DIMVECT; p++)
        //   {
        //     for (int j=0; j<nx; j++)
        //       CoeffY(j) = alpha(alpha.base(0) + j)(p);
        //     Image.projConvProj(x, CoeffY, x, CoeffY, fftKer);
        //     for (int i=0; i<nx; i++)
        //       Aux1(i) += CoeffY(i)*eta(eta.base(0) + i)(p);
        //   }
        //
        //   for (int dp=0; dp<DIMPOINT; dp++)
        //   {
        //     for (int i=0; i<nx; i++) {beta(beta.base(0)+i)(dp) = -Aux1(i)*x(x.base(0) + i)(dp);}
        //     Aux2 = 0.0;
        //     for (int dv=0; dv<DIMVECT; dv++)
        //     {
        //       for (int j=0; j<nx; j++)
        //         CoeffY(j) = alpha(alpha.base(0) + j)(dv)*x(x.base(0) + j)(dp);
        //       Image.projConvProj(x, CoeffY, x, CoeffY, fftKer);
        //       for (int i=0; i<nx; i++)
        //         Aux2(i) += CoeffY(i)*eta(eta.base(0) + i)(dv);
        //     }
        //     for (int i=0; i<nx; i++) {beta(beta.base(0)+i)(dp) += Aux2(i);}
        //   }
        //   beta *= 2*ooSigma2;
        //
        //   // cout << "ConvSpecDiffTKer = " << beta << endl;
        //   return beta;
        // }
};

#endif // GRIDGAUSSKERNEL
