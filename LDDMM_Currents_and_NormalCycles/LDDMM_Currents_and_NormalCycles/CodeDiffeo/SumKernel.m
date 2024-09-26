function ker = SumKernel(kers)

% sum of kernels
nk = length(kers);

nx = [];
dx = [];

dxo = [];
onex = [];
U = [];

    function gamma = Eval(y,x,alpha)
        gamma = 0;
        for k = 1:nk
            gamma = gamma + kers{k}.Eval(y,x,alpha);
        end
    end
ker.Eval = @Eval;

    function K = Matrix(y,x)
        K = 0;
        for k = 1:nk
            K = K + kers{k}.Matrix(y,x);
        end
    end
ker.Matrix = @Matrix;

    function G = Grad(x,alpha,beta) 
        G = 0;
        for k = 1:nk
            G = G + kers{k}.Grad(x,alpha,beta);
        end
    end
ker.Grad = @Grad;

    function D = Diff(x,alpha,eta) 
        D = 0;
        for k = 1:nk
            D = D + kers{k}.Diff(x,alpha,eta);
        end
    end
ker.Diff = @Diff;


    function H = Hess(x,alpha,eta) 
        H = 0;
        for k = 1:nk
            H = H + kers{k}.Hess(x,alpha,eta);
        end
    end
ker.Hess = @Hess;



% the following functions make use of the precomputation part of the kernel
% function to save time when computing the kernel and its derivatives on
% the same data

    function Precomp(x)
        for k = 1:nk
            kers{k}.Precomp(x);
        end
    end
ker.Precomp = @Precomp;

    function gamma = Eval_Precomp(alpha)
        gamma = 0;
        if nargin==0
            for k=1:nk
                gamma = gamma + kers{k}.Eval_Precomp();
            end
        else
            for k=1:nk
                gamma = gamma + kers{k}.Eval_Precomp(alpha);
            end
        end
    end
ker.Eval_Precomp = @Eval_Precomp;


    function G = Grad_Precomp(alpha,beta)
        G = 0;
        for k=1:nk
            G = G + kers{k}.Grad_Precomp(alpha,beta);
        end
    end
ker.Grad_Precomp = @Grad_Precomp;


    function D = Diff_Precomp(alpha,eta)
        D = 0;
        for k=1:nk
            D = D + kers{k}.Diff_Precomp(alpha,eta);
        end
    end
ker.Diff_Precomp = @Diff_Precomp;


    function H = Hess_Precomp(alpha,eta)
        H = 0;
        for k=1:nk
            H = H + kers{k}.Hess_Precomp(alpha,eta);
        end
    end
ker.Hess_Precomp = @Hess_Precomp;


end
