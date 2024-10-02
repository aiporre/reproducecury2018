
#define SHOW(x) cout<<#x "="<<x<<endl;

#include <ctime>
#include <fstream>
#include <random/uniform.h>
using namespace ranlib;

#define __DIMPOINT__ 3
#define __DIMVECT__ 3
#define __TYPE__ float

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
using namespace blitz;

#include "SqDistScalarKernel.h"
#include "CommonFunctions.h"
#include "FastGaussKernel.h"

typedef TinyVector<__TYPE__,__DIMVECT__> Vect;
typedef Array<Vect,1> ArrVect;
typedef TinyVector<__TYPE__,__DIMPOINT__> Point;
typedef Array<Point,1> ArrPoint;
typedef TinyVector<__TYPE__,__DIMPOINT__> VectPoint;
typedef Array<Point,1> ArrVectPoint;

int main (int argc, char *argv[])
{

    (void)argc;

    cout << endl << "Starting program " << argv[0] << endl;
    cout << endl << "Parameters: Sigma Nx Ny epsilon" << endl << endl;

    double tic, toc;

    int Nx=1000, Ny=1000;
    double epsilon=1e-5;
    __TYPE__ Sigma = .25;

    if(argc==5)
    {
        int cnt = 1;
        Sigma = atof (argv[cnt++]);
        Nx = atoi (argv[cnt++]);
        Ny = atoi (argv[cnt++]);
        epsilon = atof (argv[cnt++]);
    }
    else
        cout << "Using default values:" << endl;
    cout << "  Sigma=" << Sigma;
    cout << ", Nx=" << Nx;
    cout << ", Ny=" << Ny;
    cout << ", epsilon=" << epsilon;
    cout << endl << endl;

    // create a Gauss kernel (direct computation)
    GaussFunction<__TYPE__> *Fct = new GaussFunction<__TYPE__>(Sigma);
    SqDistScalarKernel<__TYPE__,__DIMPOINT__,TinyVector<__TYPE__,__DIMVECT__> > DirKer(Fct);

    // create a Gauss kernel (FGT)
    FastGaussKernel<__TYPE__,__DIMPOINT__,__DIMVECT__> FgtKer(Sigma,epsilon);

    // set a random generator
    Uniform<__TYPE__> unif;
    unif.seed((unsigned int)time(0));
    for(int i=1; i<10000; i++)
        unif.random();

    // create x,y,alpha vectors with random values
    ArrPoint x(Range(1,Nx)), y(Range(1,Ny));
    ArrVect alpha(Range(1,Nx)), gamma(Range(1,Ny));
    ArrVect gammadir(Range(1,Ny)), gammafgt(Range(1,Ny));
    ArrVectPoint betadir(Range(1,Nx)), betafgt(Range(1,Nx));

    for(int i=1; i<Nx+1; i++)
        for(int d=0; d<__DIMPOINT__; d++)
            x(i)(d) = unif.random();
    for(int i=1; i<Ny+1; i++)
        for(int d=0; d<__DIMPOINT__; d++)
            y(i)(d) = unif.random();
    for(int i=1; i<Nx+1; i++)
        for(int d=0; d<__DIMVECT__; d++)
            alpha(i)(d) = unif.random();
    for(int i=1; i<Ny+1; i++)
        for(int d=0; d<__DIMVECT__; d++)
            gamma(i)(d) = unif.random();


    cout << endl << "EvalConv:" << endl;
    // compute direct convolution

    tic = (double) clock () / (double) CLOCKS_PER_SEC;
    gammadir = DirKer.EvalConv(y,x,alpha);
    toc = (double) clock () / (double) CLOCKS_PER_SEC;
    cout << toc - tic << " seconds for direct computation" << endl;

    // compute FGT convolution
    tic = (double) clock () / (double) CLOCKS_PER_SEC;
    gammafgt = FgtKer.EvalConv(y,x,alpha);
    toc = (double) clock () / (double) CLOCKS_PER_SEC;
    cout << toc - tic << " seconds for FGT" << endl;

    // compare results
    double err = 0, nrm = 0;
    for(int i = 1; i<Ny+1; i++)
    {
        Vect v = gammadir(i)-gammafgt(i);
        err += sqrt(sum(v*v));
        v = gammadir(i);
        nrm += sqrt(sum(v*v));
    }
    err /= Ny;
    nrm /= Ny;
    cout << "Erreur moyenne (relative): " << err << " (" << err/nrm << ")" << endl;



    cout << endl << "Grad1Conv:" << endl;
    // compute direct convolution
    tic = (double) clock () / (double) CLOCKS_PER_SEC;
    betadir = DirKer.Grad1Conv(alpha,x,y,gamma);
    toc = (double) clock () / (double) CLOCKS_PER_SEC;
    cout << toc - tic << " seconds for direct computation" << endl;

    // compute fgt convolution
    tic = (double) clock () / (double) CLOCKS_PER_SEC;
    betafgt = FgtKer.Grad1Conv(alpha,x,y,gamma);
    toc = (double) clock () / (double) CLOCKS_PER_SEC;
    cout << toc - tic << " seconds for FGT" << endl;

    // compare results
    err = 0;
    nrm = 0;
    for(int i = 1; i<Nx+1; i++)
    {
        VectPoint v = (betadir(i)-betafgt(i));
        err += sqrt(sum(v*v));
        v = betadir(i);;
        nrm += sqrt(sum(v*v));
    }
    err /= Nx;
    nrm /= Nx;
    cout << "Erreur moyenne (relative): " << err << " (" << err/nrm << ")" << endl;

    if(Nx==Ny)
    {
        cout << endl << "GradConv:" << endl;
        // compute direct convolution
        tic = (double) clock () / (double) CLOCKS_PER_SEC;
        betadir = DirKer.GradConv(alpha,x,gamma);
        toc = (double) clock () / (double) CLOCKS_PER_SEC;
        cout << toc - tic << " seconds for direct computation" << endl;

        // compute FGT convolution
        tic = (double) clock () / (double) CLOCKS_PER_SEC;
        betafgt = FgtKer.GradConv(alpha,x,gamma);
        toc = (double) clock () / (double) CLOCKS_PER_SEC;
        cout << toc - tic << " seconds for FGT" << endl;

        // compare results
        err = 0;
        nrm = 0;
        for(int i = 1; i<Nx+1; i++)
        {
            VectPoint v = (betadir(i)-betafgt(i));
            err += sqrt(sum(v*v));
            v = betadir(i);
            nrm += sqrt(sum(v*v));
        }
        err /= Nx;
        nrm /= Nx;
        cout << "Erreur moyenne (relative): " << err << " (" << err/nrm << ")" << endl;
    }
    cout << endl;
}
