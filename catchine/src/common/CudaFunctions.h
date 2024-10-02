




template < typename TYPE >
class RadialFunction
{
    
public:
    
    virtual __device__ __forceinline__ TYPE Eval(TYPE u) = 0;
    virtual __device__ __forceinline__ TYPE Diff(TYPE u) = 0;
    virtual __device__ __forceinline__ TYPE Diff2(TYPE u) = 0;
    virtual __device__ __forceinline__ void DiffDiff2(TYPE u, TYPE* d1, TYPE* d2)
    {
        *d1 = Diff(u);
        *d2 = Diff2(u);
    }
};




template < typename TYPE >
class CauchyFunction : public RadialFunction<TYPE>
{
    TYPE Sigma, ooSigma2, ooSigma4;
    
    public:
    
    CauchyFunction() { }
    
    CauchyFunction(TYPE sigma)
    {
        Sigma = sigma;
        ooSigma2 = 1.0/pow(sigma,2);
        ooSigma4 = pow(ooSigma2,2);
    }
    
    __device__ __forceinline__ TYPE Eval(TYPE r2)
    {
        return 1.0/(1.0+r2*ooSigma2);
    }
    
    __device__ __forceinline__ TYPE Diff(TYPE r2)
    {
        TYPE u = 1.0+r2*ooSigma2;
        return - ooSigma2 / (u*u);
    }
    
    __device__ __forceinline__ TYPE Diff2(TYPE r2)
    {
        TYPE u = 1.0+r2*ooSigma2;
        return 2.0 * ooSigma4 / (u*u*u);
    }
    
    __device__ __forceinline__ void DiffDiff2(TYPE r2, TYPE* d1, TYPE* d2)
    {
        TYPE u = 1.0/(1.0+r2*ooSigma2);
        *d1 = - ooSigma2 * u * u;
        *d2 = - 2.0 * ooSigma2 * *d1 * u;
    }
};

template < typename TYPE >
class GaussFunction : public RadialFunction<TYPE>
{
    TYPE Sigma, ooSigma2, ooSigma4;
    
    public:
    
    GaussFunction() { }
    
    GaussFunction(TYPE sigma)
    {
        Sigma = sigma;
        ooSigma2 = 1.0/pow(sigma,2);
        ooSigma4 = pow(ooSigma2,2);
    }
    
    __device__ __forceinline__ TYPE Eval(TYPE r2)
    {
        return exp(-r2*ooSigma2);
    }
    
    __device__ __forceinline__ TYPE Diff(TYPE r2)
    {
        return - ooSigma2 * exp(-r2*ooSigma2);
    }
    
    __device__ __forceinline__ TYPE Diff2(TYPE r2)
    {
        return ooSigma4 * exp(-r2*ooSigma2);
    }
    
    __device__ __forceinline__ void DiffDiff2(TYPE r2, TYPE* d1, TYPE* d2)
    {
        *d1 = - ooSigma2 * exp(-r2*ooSigma2);
        *d2 = - ooSigma2 * *d1;
    }
};

template < typename TYPE >
class SumGaussFunction : public RadialFunction<TYPE>
{
	int Nfuns;
        TYPE *Sigmas, *ooSigma2s, *ooSigma4s;
        TYPE *Weights;
    
    public:
    
    SumGaussFunction() { }

    SumGaussFunction(int nfuns, TYPE* weights, TYPE* sigmas)
    {
        Nfuns = nfuns;
        cudaMalloc((void**)&Weights, sizeof(TYPE)*Nfuns);
        cudaMalloc((void**)&Sigmas, sizeof(TYPE)*Nfuns);
        cudaMalloc((void**)&ooSigma2s, sizeof(TYPE)*Nfuns);
        cudaMalloc((void**)&ooSigma4s, sizeof(TYPE)*Nfuns);
    	cudaMemcpy(Weights, weights, sizeof(TYPE)*Nfuns, cudaMemcpyHostToDevice);
    	cudaMemcpy(Sigmas, sigmas, sizeof(TYPE)*Nfuns, cudaMemcpyHostToDevice);
	TYPE *oosigma2s = (TYPE*)malloc(sizeof(TYPE)*Nfuns);
	TYPE *oosigma4s = (TYPE*)malloc(sizeof(TYPE)*Nfuns);
        for(int i=0; i<Nfuns; i++)
        {
            oosigma2s[i] = 1.0/pow(sigmas[i],2);
	    oosigma4s[i] = pow(oosigma2s[i],2);
        }
    	cudaMemcpy(ooSigma2s, oosigma2s, sizeof(TYPE)*Nfuns, cudaMemcpyHostToDevice);
    	cudaMemcpy(ooSigma4s, oosigma4s, sizeof(TYPE)*Nfuns, cudaMemcpyHostToDevice);
	free(oosigma2s);
	free(oosigma4s);
    }
      
    ~SumGaussFunction()
    {
    	cudaFree(Weights);
    	cudaFree(Sigmas);
    	cudaFree(ooSigma2s);
    	cudaFree(ooSigma4s);
    }
    
    __device__ __forceinline__ TYPE Eval(TYPE r2)
    {
        TYPE res = 0.0;
        for(int i=0; i<Nfuns; i++)
            res += Weights[i] * exp(-r2*ooSigma2s[i]);
        return res;
    }
    
    __device__ __forceinline__ TYPE Diff(TYPE r2)
    {
        TYPE res = 0.0;
        for(int i=0; i<Nfuns; i++)
            res += Weights[i] * (- ooSigma2s[i] * exp(-r2*ooSigma2s[i]));
        return res;
    }

    __device__ __forceinline__ TYPE Diff2(TYPE r2)
    {
        TYPE res = 0.0;
        for(int i=0; i<Nfuns; i++)
            res += Weights[i] * (ooSigma4s[i] * exp(-r2*ooSigma2s[i]));
        return res;
    }

    __device__ __forceinline__ void DiffDiff2(TYPE r2, TYPE* d1, TYPE* d2)
    {
        TYPE tmp;
        *d1 = 0.0;
        *d2 = 0.0;
        for(int i=0; i<Nfuns; i++)
        {
            tmp = - ooSigma2s[i] * exp(-r2*ooSigma2s[i]);
	    *d1 += Weights[i] * tmp;
            *d2 += Weights[i] * (- ooSigma2s[i] * tmp);
        }
    }

};

template < typename TYPE >
class SumCauchyFunction : public RadialFunction<TYPE>
{
	int Nfuns;
        TYPE *Sigmas, *ooSigma2s, *ooSigma4s;
        TYPE *Weights;
    
    public:
    
    SumCauchyFunction() { }

    SumCauchyFunction(int nfuns, TYPE* weights, TYPE* sigmas)
    {
        Nfuns = nfuns;
        cudaMalloc((void**)&Weights, sizeof(TYPE)*Nfuns);
        cudaMalloc((void**)&Sigmas, sizeof(TYPE)*Nfuns);
        cudaMalloc((void**)&ooSigma2s, sizeof(TYPE)*Nfuns);
        cudaMalloc((void**)&ooSigma4s, sizeof(TYPE)*Nfuns);
    	cudaMemcpy(Weights, weights, sizeof(TYPE)*Nfuns, cudaMemcpyHostToDevice);
    	cudaMemcpy(Sigmas, sigmas, sizeof(TYPE)*Nfuns, cudaMemcpyHostToDevice);
	TYPE *oosigma2s = (TYPE*)malloc(sizeof(TYPE)*Nfuns);
	TYPE *oosigma4s = (TYPE*)malloc(sizeof(TYPE)*Nfuns);
        for(int i=0; i<Nfuns; i++)
        {
            oosigma2s[i] = 1.0/pow(sigmas[i],2);
	    oosigma4s[i] = pow(oosigma2s[i],2);
        }
    	cudaMemcpy(ooSigma2s, oosigma2s, sizeof(TYPE)*Nfuns, cudaMemcpyHostToDevice);
    	cudaMemcpy(ooSigma4s, oosigma4s, sizeof(TYPE)*Nfuns, cudaMemcpyHostToDevice);
	free(oosigma2s);
	free(oosigma4s);
    }
      
    ~SumCauchyFunction()
    {
    	cudaFree(Weights);
    	cudaFree(Sigmas);
    	cudaFree(ooSigma2s);
    	cudaFree(ooSigma4s);
    }
    
    __device__ __forceinline__ TYPE Eval(TYPE r2)
    {
        TYPE res = 0.0;
        for(int i=0; i<Nfuns; i++)
            res += Weights[i] / (1.0+r2*ooSigma2s[i]);
        return res;
    }
           
    __device__ __forceinline__ TYPE Diff(TYPE r2)
    {
        TYPE res = 0.0, u;
        for(int i=0; i<Nfuns; i++)
	{
		u = 1.0+r2*ooSigma2s[i];
	        res += Weights[i] * (- ooSigma2s[i] / (u*u));
	}
        return res;
    }

    __device__ __forceinline__ TYPE Diff2(TYPE r2)
    {
        TYPE res = 0.0, u;
        for(int i=0; i<Nfuns; i++)
	{
		u = 1.0+r2*ooSigma2s[i];
        	res += (Weights[i] * 2.0 * ooSigma4s[i]) / (u*u*u);
	}
        return res;
    }

    __device__ __forceinline__ void DiffDiff2(TYPE r2, TYPE* d1, TYPE* d2)
    {
        TYPE u, tmp;
        *d1 = 0.0;
        *d2 = 0.0;
        for(int i=0; i<Nfuns; i++)
        { 
		u = 1.0/(1.0+r2*ooSigma2s[i]);
        	tmp = - ooSigma2s[i] * u * u;
		*d1 += Weights[i] * tmp;
        	*d2 += Weights[i] * (-2.0 * ooSigma2s[i] * tmp * u);
        }
    }

 

};



