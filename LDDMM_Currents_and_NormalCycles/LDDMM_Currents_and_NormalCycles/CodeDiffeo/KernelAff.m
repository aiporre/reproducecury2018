function ker = KernelAff(ker0)

ker.id = 'KernelAff';

    function gamma = Eval(y,x,alpha)
        % computes gamma = K(y,x) alpha.        
        gamma = Matrix(y,x) * alpha;
    end

ker.Eval = @Eval;

    function K = Matrix(y,x)
        % computes K(y,x).        
        d = size(x,2);
        nx = size(x,1)-d-1;
        ny = size(y,1)-d-1;
        x = x(1:nx,:);
        y = y(1:ny,:);
        onex = ones(1,nx);
        oney = ones(ny,1);
        K = [ker0.Matrix(y,x),y,oney;[x';onex],zeros(d+1)];
    end

ker.Matrix = @Matrix;

    function G = Grad(x,alpha,beta)         
        % gradient of alpha^T*K(x,x)*beta wrt x
        d = size(x,2);
        nx = size(x,1)-d-1;
        alphax = alpha(1:nx,:);
        alphad = alpha(nx+1:end-1,1:d);
        betax = beta(1:nx,:);
        betad = beta(nx+1:end-1,:);
        x = x(1:nx,:);
        G = [ker0.Grad(x,alphax,betax) + alphax*betad' + betax*alphad';...
            zeros(d+1,d)];        
    end

ker.Grad = @Grad;

% dummy functions Precomp...
x_ = [];
    function Precomp(x)
        x_ = x;
    end
ker.Precomp = @Precomp;

    function gamma = Eval_Precomp(alpha)
        gamma = Eval(x_,x_,alpha);
    end
ker.Eval_Precomp = @Eval_Precomp;

    function G = Grad_Precomp(alpha,beta)
        G = Grad(x_,alpha,beta);
    end
ker.Grad_Precomp = @Grad_Precomp;

end