function Fun = TPSFunction

% This implements the function f(u)=.5*u*log(u) and its first and second order derivatives
%
% example of use :
%
% fun = TPSFunction;
% u = rand(1000);   % 1000*1000 matrix of random numbers
% fu = fun.Eval(u);    % computes f(u_ij) for all entries of u
% fpu = fun.Diff(u);    % computes f'(u_ij) for all entries of u
% fppu = fun.Diff2(u);   % computes f''(u_ij) for all entries of u
%
% To speed up computations when evaluating the function and its derivatives
% on the same data, you can use the "Precomp" functions. For example, the
% following lines does the same as previously but could be potentially
% faster for large matrices
%
% fun = TPSFunction;
% u = rand(1000);   % 1000*1000 matrix of random numbers
% fun.Precomp(u);
% fu = fun.Eval_Precomp();
% fpu = fun.Diff_Precomp();
% fppu = fun.Diff2_Precomp();

    function f = Eval(u)
        f = .5*u.*log(u);
        f(u==0) = 0;
    end

Fun.Eval = @Eval;

    function d = Diff(u)
        d = .5*(log(u)+1);
        d(u==0) = 0;
    end

Fun.Diff = @Diff;

    function d2 = Diff2(u)
        d2 = .5./u;
        d2(u==0) = 0;
    end

Fun.Diff2 = @Diff2;

eval_p = [];

    function Precomp(u)
        eval_p = u;
    end

Fun.Precomp = @Precomp;

    function f = Eval_Precomp
        f = Eval(eval_p);
    end

Fun.Eval_Precomp = @Eval_Precomp;

    function f = Diff_Precomp
        f = Diff(eval_p);
    end

Fun.Diff_Precomp = @Diff_Precomp;

    function f = Diff2_Precomp
        f = Diff2(eval_p);
    end

Fun.Diff2_Precomp = @Diff2_Precomp;


end
