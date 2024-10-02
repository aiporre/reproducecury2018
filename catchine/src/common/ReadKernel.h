#ifdef KERNEL

template < typename TYPE, int DIMPOINT, int DIMVECT >
Kernel< TinyVector<TYPE,DIMPOINT> , TinyVector<TYPE,DIMVECT> , TinyVector<TYPE,DIMPOINT> >* ReadKernel(ifstream &f)
{
    typedef TinyVector<TYPE,DIMPOINT> Point;
    typedef Array<Point,1> ArrPoint;
    typedef TinyVector<TYPE,DIMVECT> Vect;
    typedef Array<Vect,1> ArrVect;
    typedef TinyVector<TYPE,DIMPOINT> VectPoint;
    typedef Array<VectPoint,1> ArrVectPoint;
    typedef Array<TYPE,1> ArrType;

    Kernel<Point,Vect,VectPoint> *Ker = NULL;
    string Tag;
    f >> Tag;
#ifdef FASTGAUSSKERNEL
    if(!Tag.compare("FastGauss,sigma,epsilon="))
    {
        TYPE sigma, epsilon;
        f >> sigma >> epsilon;
        Ker = new FastGaussKernel<TYPE,DIMPOINT,DIMVECT>(sigma,epsilon);
    }
    else
#endif // FASTGAUSSKERNEL
#ifdef GRIDGAUSSKERNEL
    if(!Tag.compare("GridGauss,sigma,ratio="))
    {
        double sigma, ratio;
       	f >> sigma >> ratio;
       	Ker = new GridGaussKernel<TYPE,DIMPOINT,DIMVECT>(sigma,ratio);
    }
    else
#endif // GRIDGAUSSKERNEL
#ifdef CAUCHYGPUKERNEL
        if(!Tag.compare("CauchyGpu,sigma="))
        {
            TYPE sigma;
            f >> sigma;
            Ker = new CauchyGpuKernel<TYPE,DIMPOINT,DIMVECT>(sigma);
        }
        else
#endif // CAUCHYGPUKERNEL
#ifdef GAUSSGPUKERNEL
            if(!Tag.compare("GaussGpu,sigma="))
            {
                TYPE sigma;
                f >> sigma;
                Ker = new GaussGpuKernel<TYPE,DIMPOINT,DIMVECT>(sigma);
            }
            else
#endif // GAUSSGPUKERNEL
#ifdef SUMGAUSSGPUKERNEL
            if(!Tag.compare("SumGaussGpu,Weights="))
            {
                ArrType weights, sigmas;
                ReadArr(weights,f);
                f >> Tag;
                ReadArr(sigmas,f);
                Ker = new SumGaussGpuKernel<TYPE,DIMPOINT,DIMVECT>(weights,sigmas);
            }
            else
#endif // SUMGAUSSGPUKERNEL
#ifdef SUMCAUCHYGPUKERNEL
            if(!Tag.compare("SumCauchyGpu,Weights="))
            {
                ArrType weights, sigmas;
                ReadArr(weights,f);
                f >> Tag;
                ReadArr(sigmas,f);
                Ker = new SumCauchyGpuKernel<TYPE,DIMPOINT,DIMVECT>(weights,sigmas);
            }
            else
#endif // SUMCAUCHYGPUKERNEL
#ifdef SQDISTSCALARKERNEL
                if(!Tag.compare("SqDistScalar,function="))
					Ker = new SqDistScalarKernel<TYPE,DIMPOINT,Vect>(f);
                else
#endif // SQDISTSCALARKERNEL 
#ifdef TRIKERNEL
                    if(!Tag.compare("Tri,functionTilde="))
                        Ker = new TriKernel<TYPE,DIMPOINT,DIMVECT>(f);
                    else
#endif // TRIKERNEL                
#ifdef CURLFREEKERNEL
                    if(!Tag.compare("CurlFree,function="))
                    {
                        Function<TYPE> *fct = ReadFunction<TYPE>(f);
                        Ker = new CurlFreeKernel<TYPE,DIMPOINT,DIMVECT>(fct);
                    }
                    else
#endif // CURLFREEKERNEL                
#ifdef DIVFREEKERNEL
                        if(!Tag.compare("DivFree,function="))
                        {
                            Function<TYPE> *fct = ReadFunction<TYPE>(f);
                            Ker = new DivFreeKernel<TYPE,DIMPOINT,DIMVECT>(fct);
                        }
                        else
#endif // DIVFREEKERNEL        
						{
                            cout << "Error constructing Kernel from file: unknown type '" << Tag << "'." << endl;
							throw -1;
						}

    return Ker;
}


#endif // KERNEL

