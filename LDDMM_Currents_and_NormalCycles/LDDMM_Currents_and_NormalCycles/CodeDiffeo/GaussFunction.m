function Fun = GaussFunction

% This implements the function f(u)=exp(-u) and its first and second order derivatives
%
% example of use :
%
% fun = GaussFunction;
% u = rand(1000);   % 1000*1000 matrix of random numbers
% fu = fun.Eval(u);    % computes f(u_ij) = exp(-u_ij) for all entries of u
% fpu = fun.Diff(u);    % computes f'(u_ij) = -exp(-u_ij) for all entries of u
% fppu = fun.Diff2(u);   % computes f''(u_ij) = exp(u_ij) for all entries of u
%
% To speed up computations when evaluating the function and its derivatives
% on the same data, you can use the "Precomp" functions. For example, the
% following lines does the same as previously but could be potentially
% faster for large matrices
%
% fun = GaussFunction;
% u = rand(1000);   % 1000*1000 matrix of random numbers
% fun.Precomp(u);
% fu = fun.Eval_Precomp();
% fpu = fun.Diff_Precomp();
% fppu = fun.Diff2_Precomp();

    function f = Eval(u)
        f = exp(-u);
    end

Fun.Eval = @Eval;

    function d = Diff(u)
        d = -exp(-u);
    end

Fun.Diff = @Diff;

    function d2 = Diff2(u)
        d2 = exp(-u);
    end

Fun.Diff2 = @Diff2;

eval_p = [];

    function Precomp(u)
        eval_p = Eval(u);
    end

Fun.Precomp = @Precomp;

    function f = Eval_Precomp
        f = eval_p;
    end

Fun.Eval_Precomp = @Eval_Precomp;

    function f = Diff_Precomp
        f = -eval_p;
    end

Fun.Diff_Precomp = @Diff_Precomp;

    function f = Diff2_Precomp
        f = eval_p;
    end

Fun.Diff2_Precomp = @Diff2_Precomp;


end
