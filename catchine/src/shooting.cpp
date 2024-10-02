#define SHOW(x) cout<<#x "="<<x<<endl;

#include <ctime>
#include <fstream>
#include <string.h>
#include <random/uniform.h>
using namespace ranlib;

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
using namespace blitz;

// Dimension of data in all code (usually 2 or 3)
#ifndef __Dim__
#define __Dim__ 3
#endif

// Use of Cuda
#ifndef CudaFlag
#define CudaFlag 0
#endif

// Use of FGT
#ifndef FgtFlag
#define FgtFlag 1
#endif

// types of functions for kernels
#include "CommonFunctions.h" // Gaussian, Cauchy, Cauchy2
//#include "SobolevFunctions.h" // kernel functions for H1 and H2 Sobolev norms
#include "ReadFunction.h" // leave this after includes of all function types

// types of kernels
#include "SqDistScalarKernel.h"
#include "TriKernel.h"
/*
#include "CurlFreeKernel.h"
#include "DivFreeKernel.h"
*/
#if FgtFlag==1
#include "FastGaussKernel.h"
#endif
#if CudaFlag==1
#include "CauchyGpuKernel.h"
#include "GaussGpuKernel.h"
#include "SumGaussGpuKernel.h"
#include "SumCauchyGpuKernel.h"
#endif
#include "ReadKernel.h" // leave this after includes of all kernel types


// types of deformations
#include "LargeDef_InitParam.h"
#include "ReadEvol.h" // leave this after includes of all deformation types

#include "shooting.h"

int main (int argc, char *argv[])
{

    cout << "Starting program " << argv[0] << endl;

    if (argc < 4)   // Check the value of argc. If not enough parameters have been passed, inform user and exit.
    {
        cout << "Usage is TypeFloat DataFile ResFile" << endl; // Inform the user of how to use the program
        exit(0);
    }

    int cnt = 1;
    string TypeFloat = argv[cnt++];
    cout << "TypeFloat=" << TypeFloat;
    char* DataFile = argv[cnt++];
    cout << ", DataFile=" << DataFile;
    char* ResFile = argv[cnt++];
    cout << ", ResFile=" << ResFile;


    cout << endl;

    if(!TypeFloat.compare("float"))
        return shooting<float>(DataFile, ResFile);
    else if(!TypeFloat.compare("double"))
        return shooting<double>(DataFile, ResFile);

}
