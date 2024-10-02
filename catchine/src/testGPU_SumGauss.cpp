
// Test program for kernel convolutions on the GPU

#define SHOW(x) cout<<#x "="<<x<<endl;

#include <ctime>
#include <fstream>
#include <random/uniform.h>
using namespace ranlib;

// Dimension of data in all code (usually 2 or 3)
#ifndef __Dim__
#define __Dim__ 3
#endif
#define __DIMPOINT__ __Dim__
#define __DIMVECT__ __Dim__

#define __TYPE__ double

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
using namespace blitz;

#include "SqDistScalarKernel.h"
#include "CommonFunctions.h"
#include "SumGaussGpuKernel.h"

typedef TinyVector<__TYPE__,__DIMVECT__> Vect;
typedef Array<Vect,1> ArrVect;
typedef TinyVector<__TYPE__,__DIMPOINT__> Point;
typedef Array<Point,1> ArrPoint;
typedef TinyVector<__TYPE__,__DIMPOINT__> VectPoint;
typedef Array<Point,1> ArrVectPoint;

int main (int argc, char *argv[])
{
	
    (void) argc;
	
    cout << endl << "Starting program " << argv[0] << endl;
    cout << endl << "Parameters: SigmaMin SigmaMax NSigmas Nx Ny" << endl << endl;
	
    double tic, toc;
	
    int Nx=10000, Ny=10000;
    int NSigmas = 3;
    __TYPE__ SigmaMin = .05;
    __TYPE__ SigmaMax = .5;
	
    if(argc==6)
    {
        int cnt = 1;
        SigmaMin = atof (argv[cnt++]);
        SigmaMax = atof (argv[cnt++]);
        NSigmas = atoi (argv[cnt++]);
        Nx = atoi (argv[cnt++]);
        Ny = atoi (argv[cnt++]);
    }
    else
    cout << "Using default values:" << endl;
    cout << "  SigmaMin=" << SigmaMin;
    cout << "  SigmaMax=" << SigmaMax;
    cout << "  NSigmas=" << NSigmas;
    cout << ", Nx=" << Nx;
    cout << ", Ny=" << Ny;
    cout << endl;
	
    if(sizeof(__TYPE__)==sizeof(float))
	cout << endl << "Using single-precision floating-point numbers" << endl << endl;
    else	
	cout << endl << "Using double-precision floating-point numbers" << endl << endl;

    // create a sum of Gauss kernels (direct computation)
    Array<Function<__TYPE__>*,1> FctArr(Range(1,NSigmas));
    Array<__TYPE__,1> Weights(Range(1,NSigmas));
    Array<__TYPE__,1> Sigmas(Range(1,NSigmas));
    Weights = 1;
    Sigmas(1) = SigmaMin;
    for(int i=2; i<NSigmas+1; i++)
	Sigmas(i) = SigmaMin+((__TYPE__)(i-1)/(NSigmas-1))*(SigmaMax-SigmaMin);
    for(int i=1; i<NSigmas+1; i++)
	FctArr(i) = new GaussFunction<__TYPE__>(Sigmas(i));
    SumFunction<__TYPE__>* Fct = new SumFunction<__TYPE__>(FctArr,Weights);
    SqDistScalarKernel<__TYPE__,__DIMPOINT__,TinyVector<__TYPE__,__DIMVECT__> > DirKer(Fct);

    // create a sum of Gauss kernels (GPU)
    SumGaussGpuKernel<__TYPE__,__DIMPOINT__,__DIMVECT__> GpuKer(Weights,Sigmas);

    // set a random generator
    Uniform<__TYPE__> unif;
    unif.seed((unsigned int)time(0));
    for(int i=1; i<10000; i++)
    	unif.random();
    // create x,y,alpha,eta vectors with random values
    ArrPoint x(Range(1,Nx)), y(Range(1,Ny));
    ArrVect alpha(Range(1,Nx)), gamma(Range(1,Ny));
    ArrVect gammadir(Range(1,Ny)), gammagpu(Range(1,Ny));
    ArrVectPoint eta(Range(1,Nx)), betadir(Range(1,Nx)), betagpu(Range(1,Nx));
	
    for(int i=1; i<Nx+1; i++)
    for(int d=0; d<__DIMPOINT__; d++)
    x(i)(d) = unif.random();
    for(int i=1; i<Ny+1; i++)
    for(int d=0; d<__DIMPOINT__; d++)
    y(i)(d) = unif.random();
    for(int i=1; i<Nx+1; i++)
    for(int d=0; d<__DIMPOINT__; d++)
    eta(i)(d) = unif.random();
    for(int i=1; i<Nx+1; i++)
    for(int d=0; d<__DIMVECT__; d++)
    alpha(i)(d) = unif.random();
    for(int i=1; i<Ny+1; i++)
    for(int d=0; d<__DIMVECT__; d++)
    gamma(i)(d) = unif.random();
    
    
	/// GPU Initialization
	
	tic = (double) clock () / (double) CLOCKS_PER_SEC;
    // do a first call to the GPU for initializing (quick turnaround; there should be something better to do !!!)

	ArrPoint dum1(Range(1,100)); ArrVect dum2(Range(1,100)); dum1 = 0.0; dum2 = 0.0; dum2 = GpuKer.EvalConv(dum1,dum1,dum2);
	toc = (double) clock () / (double) CLOCKS_PER_SEC;
	double tmp = toc - tic;
	cout << tmp << " seconds for GPU init" << endl;
    

    
    
	/// Test of EvalConv
	
	cout << endl << "EvalConv:" << endl;
    // compute gpu convolution
    tic = (double) clock () / (double) CLOCKS_PER_SEC;
    gammagpu = GpuKer.EvalConv(y,x,alpha);
    toc = (double) clock () / (double) CLOCKS_PER_SEC;
	double timeGpu = toc - tic;
    cout << timeGpu << " seconds for GPU" << endl;
    
    // compute direct convolution
    tic = (double) clock () / (double) CLOCKS_PER_SEC;
    gammadir = DirKer.EvalConv(y,x,alpha);
    toc = (double) clock () / (double) CLOCKS_PER_SEC;
	double timeDir = toc - tic;
	cout << timeDir << " seconds for direct computation" << endl;
	cout << "speedup : " << timeDir/timeGpu << " x" << endl;
	
    // compare results
    double err = 0, errmax = 0, nrm = 0;
    for(int i = 1; i<Ny+1; i++)
    {
        Vect v = gammadir(i)-gammagpu(i);
        tmp = sqrt(sum(v*v));
        err += tmp;
        errmax = max(tmp,errmax);
        v = gammadir(i);
        nrm += sqrt(sum(v*v));
    }
    err /= Ny;
    nrm /= Ny;
    cout << "Erreur moyenne (relative): " << err << " (" << err/nrm << ")" << endl;
    cout << "Erreur max (relative): " << errmax << " (" << errmax/nrm << ")" << endl;
	

	/// Test of Grad1Conv
	
    cout << endl << "Grad1Conv:" << endl;
    // compute direct convolution
    tic = (double) clock () / (double) CLOCKS_PER_SEC;
    betadir = DirKer.Grad1Conv(alpha,x,y,gamma);
    toc = (double) clock () / (double) CLOCKS_PER_SEC;
	timeDir = toc - tic;
	cout << timeDir << " seconds for direct computation" << endl;
	
    // compute gpu convolution
    tic = (double) clock () / (double) CLOCKS_PER_SEC;
    betagpu = GpuKer.Grad1Conv(alpha,x,y,gamma);
    toc = (double) clock () / (double) CLOCKS_PER_SEC;
	timeGpu = toc - tic;
    cout << timeGpu << " seconds for GPU" << endl;
	cout << "speedup : " << timeDir/timeGpu << " x" << endl;
	
    // compare results
    err = 0;
    nrm = 0;
    for(int i = 1; i<Nx+1; i++)
    {
        VectPoint v = (betadir(i)-betagpu(i));
        err += sqrt(sum(v*v));
        v = betadir(i);;
        nrm += sqrt(sum(v*v));
    }
    err /= Nx;
    nrm /= Nx;
    cout << "Erreur moyenne (relative): " << err << " (" << err/nrm << ")" << endl;
    
    
    
	/// Test of GradConv
	
    if(Nx==Ny)
    {
        cout << endl << "GradConv:" << endl;
        // compute direct convolution
        tic = (double) clock () / (double) CLOCKS_PER_SEC;
        betadir = DirKer.GradConv(alpha,x,gamma);
        toc = (double) clock () / (double) CLOCKS_PER_SEC;
		timeDir = toc - tic;
		cout << timeDir << " seconds for direct computation" << endl;
        // compute GPU convolution
        tic = (double) clock () / (double) CLOCKS_PER_SEC;
        betagpu = GpuKer.GradConv(alpha,x,gamma);
        toc = (double) clock () / (double) CLOCKS_PER_SEC;
		timeGpu = toc - tic;
		cout << timeGpu << " seconds for GPU" << endl;
		
		cout << "speedup : " << timeDir/timeGpu << " x" << endl;
		
        // compare results
        err = 0;
        nrm = 0;
        for(int i = 1; i<Nx+1; i++)
        {
            VectPoint v = (betadir(i)-betagpu(i));
            err += sqrt(sum(v*v));
            v = betadir(i);
            nrm += sqrt(sum(v*v));
        }
        err /= Nx;
        nrm /= Nx;
        cout << "Erreur moyenne (relative): " << err << " (" << err/nrm << ")" << endl;
    }
    
	
	/// Test of DiffConv(ArrPoint &x, ArrVect &beta, ArrVectPoint &eta)
	
    if(Nx==Ny)
    {
        cout << endl << "DiffConv:" << endl;
        // compute direct convolution
        tic = (double) clock () / (double) CLOCKS_PER_SEC;
        gammadir = DirKer.DiffConv(x,alpha,eta);
        toc = (double) clock () / (double) CLOCKS_PER_SEC;
		timeDir = toc - tic;
		cout << timeDir << " seconds for direct computation" << endl;
        // compute GPU convolution
        tic = (double) clock () / (double) CLOCKS_PER_SEC;
        gammagpu = GpuKer.DiffConv(x,alpha,eta);
        toc = (double) clock () / (double) CLOCKS_PER_SEC;
		timeGpu = toc - tic;
		cout << timeGpu << " seconds for GPU" << endl;
		
		cout << "speedup : " << timeDir/timeGpu << " x" << endl;
		
        // compare results
        err = 0;
        nrm = 0;
        for(int i = 1; i<Nx+1; i++)
        {
            Vect v = (gammadir(i)-gammagpu(i));
            err += sqrt(sum(v*v));
            v = betadir(i);
            nrm += sqrt(sum(v*v));
        }
        err /= Nx;
        nrm /= Nx;
        cout << "Erreur moyenne (relative): " << err << " (" << err/nrm << ")" << endl;
		
	}
    
	/// Test of GradDiffConv(ArrPoint &x, ArrVect &beta, ArrVectPoint &eta)
	
    if(Nx==Ny)
    {
        cout << endl << "GradDiffConv:" << endl;
        // compute direct convolution
        tic = (double) clock () / (double) CLOCKS_PER_SEC;
        betadir = DirKer.GradDiffConv(x,alpha,eta);
        toc = (double) clock () / (double) CLOCKS_PER_SEC;
		timeDir = toc - tic;
		cout << timeDir << " seconds for direct computation" << endl;
        // compute GPU convolution
        tic = (double) clock () / (double) CLOCKS_PER_SEC;
        betagpu = GpuKer.GradDiffConv(x,alpha,eta);
        toc = (double) clock () / (double) CLOCKS_PER_SEC;
		timeGpu = toc - tic;
		cout << timeGpu << " seconds for GPU" << endl;
		
		cout << "speedup : " << timeDir/timeGpu << " x" << endl;
		
        // compare results
        err = 0;
        nrm = 0;
        for(int i = 1; i<Nx+1; i++)
        {
            VectPoint v = (betadir(i)-betagpu(i));
            err += sqrt(sum(v*v));
            v = betadir(i);
            nrm += sqrt(sum(v*v));
        }
        err /= Nx;
        nrm /= Nx;
        cout << "Erreur moyenne (relative): " << err << " (" << err/nrm << ")" << endl;
		
	}
	
    cout << endl;
    delete(Fct);
}
