
#define DLIB_CATCHINE 1

#ifdef DLIB_CATCHINE 
#include <dlib/optimization.h>
using namespace dlib;
#endif


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
#endif

#include "ReadKernel.h" // leave this after includes of all kernel types

// types of deformations
#include "LargeDef_InitParam.h"
#include "SmallDef.h"
#include "LargeDef.h"
#include "FreeEvol.h"
/*
#include "LargeDefSpec.h"
*/
#include "ReadEvol.h" // leave this after includes of all deformation types

// types of matchings
#include "Landmarks.h" // landmark matching
#include "Measure.h" // measure matching
#include "CurveCurr.h" // curve matching
#if __Dim__==3
#include "L2Image.h" // L2 image matching
#include "SurfCurr.h" // surface matching
#endif
#include "ReadTarget.h" // leave this after includes of all target types

// general matching class
#include "Match.h"

#include "matching.h"

int main (int argc, char *argv[])
{

    // All default argument values
    string TypeFloat("double");
    string DataFile("tmp.mch");
    string ResFile("out.mch");
    double RegWeight = 0.;
    int useoptim = 1;
    int niters = 500;
    double StepSize = 1.;
    double breakratio = 1.e-10;
    int loopbreak = 10;
    int verbosemode = 0;
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
        // Output data file
        else if (strcmp(argv[i], "-o") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid output data file " << argv[i] << std::endl;
                return -1;
            }
            ResFile = argv[++i];
        }
        // Weight of regularization term
        else if (strcmp(argv[i], "-w") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid weight of regularization term " << argv[i] << std::endl;
                return -1;
            }
            RegWeight = atof(argv[++i]);
        }
        // type of gradient descent
        else if (strcmp(argv[i], "-u") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid type of gradient descent " << argv[i] << std::endl;
                return -1;
            }
            useoptim = atoi(argv[++i]);
        }
        // Maximum number of iterations
        else if (strcmp(argv[i], "-i") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid number of iterations " << argv[i] << std::endl;
                return -1;
            }
            niters = atoi(argv[++i]);
        }
        // Initial gradient descent step
        else if (strcmp(argv[i], "-s") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid initial gradient descent step " << argv[i] << std::endl;
                return -1;
            }
            StepSize = atof(argv[++i]);
        }
        // Parameter 1 for stopping rule in gradient descent
        else if (strcmp(argv[i], "-b") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid parameter 1 for stopping rule in gradient descent " << argv[i] << std::endl;
                return -1;
            }
            breakratio = atof(argv[++i]);
        }
        // Parameter 2 for stopping rule in gradient descent
        else if (strcmp(argv[i], "-l") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid parameter 2 for stopping rule in gradient descent " << argv[i] << std::endl;
                return -1;
            }
            loopbreak = atoi(argv[++i]);
        }
        // Verbose Mode : display info at each iteration of gradient descent (0,1,2)
        else if (strcmp(argv[i], "-v") == 0)
        {
            if (i+1 == argc)
            {
                std::cout << "Invalid verbose mode " << argv[i] << std::endl;
                return -1;
            }
            verbosemode = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            std::cout << "match options:" << std::endl;
            std::cout << "\t-f ARG\t:\ttype of floating point accuracy used in all code: float or double -- default=double" << std::endl;
            std::cout << "\t-d ARG\t:\tname of  input data file -- default=tmp.mch" << std::endl;
            std::cout << "\t-o ARG\t:\tname of output data file -- default=out.mch" << std::endl;
            std::cout << "\t-w ARG\t:\tweight of regularization term -- default=0" << std::endl;
            std::cout << "\t-u ARG\t:\ttype of gradient descent used: 0 (fixed gradient descent) or 1 (adaptive stepsize) or 2 (LBFGS) -- default=1" << std::endl;
            std::cout << "\t-v ARG\t:\tverbose mode: 0 (no output except on error), 1 (report status of process), or 2 (report status and detail optimization steps)   -- default=1" << std::endl;
            std::cout << "\t-i ARG\t:\tmaximum number of iterations in gradient descent -- default=500" << std::endl;
            std::cout << "\t-s ARG\t:\tinitial gradient descent step -- default=1" << std::endl;
            std::cout << "\t-b ARG\t:\tparameter 1 for stopping rule in gradient descent -- default=1e-10" << std::endl;
            std::cout << "\t-l ARG\t:\tparameter 2 for stopping rule in gradient descent -- default=10" << std::endl;
            std::cout << "\t-h\t:\tprint this output" << std::endl;
            return 0;
        }
    }

    if(verbosemode) cout << "Starting program : match" << endl;

// Rest of the code is templatized ; hence put in a separate .h file

    if(!TypeFloat.compare("float"))
        return matching<float>(DataFile.c_str(), ResFile.c_str(), niters, StepSize, RegWeight, useoptim, breakratio, loopbreak, verbosemode);
    else if(!TypeFloat.compare("double"))
        return matching<double>(DataFile.c_str(), ResFile.c_str(), niters, StepSize, RegWeight, useoptim, breakratio, loopbreak, verbosemode);
}

