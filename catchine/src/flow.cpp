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
#define FgtFlag 0
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
#include "LargeDef.h"
#include "SmallDef.h"
/*
#include "LargeDefSpec.h"
#include "FreeEvol.h"
*/
#include "ReadEvol.h" // leave this after includes of all deformation types

#include "flow.h"

int main (int argc, char *argv[])
{
    cout << "Starting program " << argv[0] << endl;

    string TypeFloat("double");
    string DataFile("points.mch");
    string ResFile("out.mch");
    string PointsFileOut("pout.mch");
	int tdeb = -1;
	int tfin = -1;

    for(int i = 1; i < argc; ++i)
    {
        // Type of float or double
        if (strcmp(argv[i], "-f") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid type float " << argv[i] << std::endl;
                return -1;
            }
            TypeFloat = argv[++i];
        }
        // Input data file
        else if (strcmp(argv[i], "-d") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid input data file " << argv[i] << std::endl;
                return -1;
            }
            DataFile = argv[++i];
        }
        // Deformation parameter file
        else if (strcmp(argv[i], "-m") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid deformation parameter file " << argv[i] << std::endl;
                return -1;
            }
            ResFile = argv[++i];
        }
        // Output data file
        else if (strcmp(argv[i], "-o") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid output data file " << argv[i] << std::endl;
                return -1;
            }
            PointsFileOut = argv[++i];
        }
        // initial time step
        else if (strcmp(argv[i], "-b") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid tdeb " << argv[i] << std::endl;
                return -1;
            }
            tdeb = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-e") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid tfin " << argv[i] << std::endl;
                return -1;
            }
            tfin = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            std::cout << "flow options:" << std::endl;
            std::cout << "\t-f ARG\t:\ttype of floating point accuracy used in all code: float or double -- default=double" << std::endl;
            std::cout << "\t-d ARG\t:\tname of  input data file -- default=tmp.mch" << std::endl;
            std::cout << "\t-o ARG\t:\tname of output data file -- default=pout.mch" << std::endl;
            std::cout << "\t-m ARG\t:\tname of deformation parameter file -- default=out.mch" << std::endl;
            std::cout << "\t-b ARG\t:\tinitial time step for flow integration -- default=-1" << std::endl;
            std::cout << "\t-e ARG\t:\t final  time step for flow integration -- default=-1" << std::endl;
            std::cout << "\t-h\t:\tprint this output" << std::endl;
            return 0;
        }
    }

    cout << endl;

    if(!TypeFloat.compare("float"))
        return flow<float>(ResFile.c_str(), DataFile.c_str(), PointsFileOut.c_str(),tdeb,tfin);
    else if(!TypeFloat.compare("double"))
        return flow<double>(ResFile.c_str(), DataFile.c_str(), PointsFileOut.c_str(),tdeb,tfin);

}
