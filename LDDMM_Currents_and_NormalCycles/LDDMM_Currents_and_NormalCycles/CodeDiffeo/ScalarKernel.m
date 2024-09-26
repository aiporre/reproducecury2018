function ker = ScalarKernel(fun,sigma,weight)

% scalar kernel : this implements convolutions with a kernel of the form
% k(x,y)= weight * h(-|x-y|^2/sigma^2)
% where x,y are in R^d and h = fun.Eval

if nargin < 3
    weight = 1;
end

oosigma = 1/sigma;

nx = [];
dx = [];

dxo = [];
onex = [];
U = [];

    function gamma = Eval(y,x,alpha)
        % computes gamma_i = sum_j K(y_i,x_j) alpha_j, for two sets
        % of points y_i, x_j in R^d, and vectors alpha_j in R^k
        % This convolution writes gamma = K(y,x) alpha in matrix form.        
        onex = ones(1,size(x,1));
        oney = ones(size(y,1),1);
        Uxy = 0;
        for k = 1:size(x,2)
            Uxy = Uxy + (oney*(oosigma*x(:,k)')-(oosigma*y(:,k))*onex).^2;
        end
        gamma = weight * (fun.Eval(Uxy)*alpha);
    end

ker.Eval = @Eval;

    function K = Matrix(y,x)
        % computes the matrix K(y,x) of the K(y_i,x_j) for two sets
        % of points y_i, x_j in R^d
        onex = ones(1,size(x,1));
        oney = ones(size(y,1),1);
        Uxy = 0;
        for k = 1:size(x,2)
            Uxy = Uxy + (oney*(oosigma*x(:,k)')-(oosigma*y(:,k))*onex).^2;
        end
        K = weight * fun.Eval(Uxy);
    end

ker.Matrix = @Matrix;

    function G = Grad(x,alpha,beta) 
        % gradient of alpha^T*K(x,x)*beta wrt x
        Precomp(x);
        G = Grad_Precomp(alpha,beta);
    end

ker.Grad = @Grad;


    function D = Diff(x,alpha,eta) 
        % differential of K(x,x)*alpha wrt x, applied to eta
        Precomp(x);
        D = Diff_Precomp(alpha,eta);
    end

ker.Diff = @Diff;


    function H = Hess(x,alpha,eta) 
        % Hessian matrix of alpha^T K(x,x) alpha, multiplied by eta
        Precomp(x);
        H = Hess_Precomp(alpha,eta);
    end

ker.Hess = @Hess;



% the following functions make use of the precomputation part of the kernel
% function to save time when computing the kernel and its derivatives on
% the same data

    function Precomp(x)
        % Given x, precomputes quantities which are shared by all
        % convolutions
        [nx,dx] = size(x);
        onex = ones(1,nx);
        U = 0;
        % the follwing loop computes the matrix U of the |x_i-x_j|^2/sigma^2
        for k = 1:dx
            xk = (oosigma*x(:,k))*ones(1,nx);
            dxo{k} = xk - xk';
            U = U + dxo{k}.^2;
        end
        % computes the h(|x_i-x_j|^2/sigma^2)
        fun.Precomp(U);
    end


ker.Precomp = @Precomp;

    function gamma = Eval_Precomp(alpha)
        % computes K(x,x) alpha where x was previously
        % loaded via function Precomp
        if nargin==0
            gamma = weight * fun.Eval_Precomp();
        else
            gamma = weight * (fun.Eval_Precomp()*alpha);
        end
    end

ker.Eval_Precomp = @Eval_Precomp;


    function G = Grad_Precomp(alpha,beta)
        G = zeros(nx,dx);
        AB = 0;
        for k = 1:size(alpha,2)
            AB = AB + alpha(:,k)*beta(:,k)';
        end
        T = (2*weight*oosigma) * fun.Diff_Precomp().*(AB+AB');
        for k = 1:dx
            G(:,k) = sum(T.*dxo{k},2);
        end
    end

ker.Grad_Precomp = @Grad_Precomp;


    function D = Diff_Precomp(alpha,eta)
        XE = 0;
        for k = 1:dx
            etaok = (oosigma*eta(:,k))*onex;
            detaok = etaok - etaok';
            XE = XE + detaok.*dxo{k};
        end
        D = (2*weight) * ((fun.Diff_Precomp().*XE)*alpha);
    end

ker.Diff_Precomp = @Diff_Precomp;


    function H = Hess_Precomp(alpha,eta)
        da = size(alpha,2);
        dU = fun.Diff_Precomp();
        A = 0;
        for k=1:da
            A = A + alpha(:,k)*alpha(:,k)';
        end
        XE = 0;
        for k=1:dx
            etaok = (oosigma*eta(:,k))*onex;
            detao{k} = etaok - etaok';
            XE = XE + detao{k}.*dxo{k};
        end
        T1 = 2 * (fun.Diff2_Precomp() .* XE);
        H = zeros(nx,dx);
        for k = 1:dx
            H(:,k) = (4*oosigma) * sum((T1 .* dxo{k} + dU .* detao{k}) .* A , 2);
        end
    end

ker.Hess_Precomp = @Hess_Precomp;


end
