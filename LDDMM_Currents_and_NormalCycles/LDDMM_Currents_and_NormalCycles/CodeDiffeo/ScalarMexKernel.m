function ker = ScalarMexKernel(name,sigma,weight,tag_gpu,dimpoint,dimvect,deviceID)

% scalar kernel using Mex files

if nargin < 7
    deviceID = 0;
end

if nargin < 3
    weight = 1;
end

if nargin < 4
    tag_gpu = '';
end

pathstr = fileparts(mfilename('fullpath'));
folder = [pathstr,'/Mex',tag_gpu,'Kernels'];
addpath(folder)
MexKernelEval = [];
MexKernelDiff = [];
MexKernelGrad = [];
MexKernelHess = [];

% compile mex functions if not done already
Compile(1)

    function Compile(docheck,item)
        if nargin<1
            docheck = 0;
        end
        if nargin<2
            for item = {'Eval','Diff','Grad','Hess'}
                Compile(docheck,item{1});
            end
        else
            switch tag_gpu
                case '' 
                    funame = [name,'Kernel',item];
                    if ~docheck || exist(funame)~=3
                        disp(['Compiling ',funame,'.c'])
                        mex(['-DFUN=s',name],'-outdir',folder,'-output',funame,[folder,'/Kernel',item,'.c'])
                    end
                    eval(['MexKernel',item,' = @',name,'Kernel',item,';']);
                case 'Gpu' 
                    funame = [name,'Gpu_',num2str(dimpoint),'_',num2str(dimvect),'_',item];
                    sourcefile = [folder,'/KernelGpuConv'];
                    if ~docheck || exist(funame)~=3
                        disp(['Compiling ',funame])
                        mexinc = '/usr/local/MATLAB/R2017b/extern/include';
                        if system(['nvcc', ...
                            ' -DKERNEL=SCALARRADIAL', ...
                            ' -DEVAL=s',item, ...
                            ' -D__TYPE__=double', ...
                            ' -DRADIALFUN=',name,'Function', ...
                            ' -D__DIMPOINT__=',num2str(dimpoint), ...
                            ' -D__DIMVECT__=',num2str(dimvect), ...
                            ' -c "',sourcefile,'.cu" ','-o "',sourcefile,'.o" -arch sm_20 -std=c++11 -Wno-deprecated-gpu-targets -Xcompiler -fPIC -I ',mexinc])
                            error('nvcc command not available; aborting compilation. Maybe your system does not support GPU computations with Cuda.')
                        end
                        mex('-lcudart','-lcufft','-L/usr/local/cuda/lib64','-cxx','-outdir',folder,'-output',funame,[sourcefile,'.o'])
                        eval(['!rm -f "',sourcefile,'.o"'])
                    end
                    eval(['MexKernel',item,' = @',funame,';']);
            end
        end
    end
ker.Compile = @Compile;

oosigma = 1/sigma;

switch tag_gpu
    case ''
        Eval = @EvalCpu;
        Grad = @GradCpu;
        Diff = @DiffCpu;
        Hess = @HessCpu;
    case 'Gpu'
        Eval = @EvalGpu;
        Grad = @GradGpu;
        Diff = @DiffGpu;
        Hess = @HessGpu;
        % first call to initialize GPU
        Eval(rand(100,dimpoint),rand(100,dimpoint),rand(100,dimvect));
end
        
    function gamma = EvalCpu(y,x,alpha) 
        % computes K(y,x) alpha
        gamma = weight * MexKernelEval(oosigma*y,oosigma*x,alpha);
    end

    function gamma = EvalGpu(y,x,alpha) 
        % computes K(y,x) alpha
        gamma = weight * MexKernelEval(oosigma*y',oosigma*x',alpha',deviceID)';
    end

ker.Eval = Eval;

    function G = GradCpu(x,alpha,beta) 
        % gradient of alpha^T*K(x,x)*beta wrt x
        G = (weight*oosigma)*MexKernelGrad(oosigma*x,alpha,beta);
    end

    function G = GradGpu(x,alpha,beta) 
        % gradient of alpha^T*K(x,x)*beta wrt x
        G = (weight*oosigma)*MexKernelGrad(alpha',oosigma*x',beta',alpha',oosigma*x',beta',deviceID)';
    end

ker.Grad = Grad;

    function D = DiffCpu(x,alpha,eta) 
        % differential of K(x,x)*alpha wrt x, applied to eta
        D = weight * MexKernelDiff(oosigma*x,alpha,oosigma*eta);
    end

    function D = DiffGpu(x,alpha,eta) 
        % differential of K(x,x)*alpha wrt x, applied to eta
        D = weight * MexKernelDiff(oosigma*x',oosigma*eta',oosigma*x',oosigma*eta',alpha',deviceID)';
    end

ker.Diff = Diff;

    function H = HessCpu(x,alpha,eta) 
        % Hessian matrix of alpha^T K(x,x) alpha, multiplied by eta
        H = (weight*oosigma) * MexKernelHess(oosigma*x,alpha,oosigma*eta);
    end

    function H = HessGpu(x,alpha,eta) 
        % Hessian matrix of alpha^T K(x,x) alpha, multiplied by eta
        H = (weight*oosigma) * MexKernelHess(oosigma*x',alpha',oosigma*eta',oosigma*x',alpha',oosigma*eta',deviceID)';
    end

ker.Hess = Hess;





% the following is the precomputation part. It does exactly the same in this code, just left for 
% compatibility with ScalarKernel

x = [];

    function Precomp(x_)
        % Given x, precomputes quantities which are shared by all
        % convolutions (does nothing here, just for compatibility)
        x = x_;
    end

ker.Precomp = @Precomp;

    function gamma = Eval_Precomp(alpha)
        gamma = Eval(x,x,alpha);
    end

ker.Eval_Precomp = @Eval_Precomp;


    function G = Grad_Precomp(alpha,beta)
        G = Grad(x,alpha,beta);
    end

ker.Grad_Precomp = @Grad_Precomp;

    function D = Diff_Precomp(alpha,eta)        
        D = Diff(x,alpha,eta);
    end

ker.Diff_Precomp = @Diff_Precomp;

    function H = Hess_Precomp(alpha,eta)
        H = Hess(x,alpha,eta);
    end

ker.Hess_Precomp = @Hess_Precomp;




end
