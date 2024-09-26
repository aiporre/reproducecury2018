function Fun = ScaledFunction(fct,lambda)

% implements g(u)=f(lambda*u)
% fct is the input function f

lambda2 = lambda^2;

    function f = Eval(u)
        f = fct.Eval(lambda*u);
    end

Fun.Eval = @Eval;

    function d = Diff(u)
        d = lambda*fct.Diff(lambda*u);
    end

Fun.Diff = @Diff;

    function d2 = Diff2(u)
        d2 = lambda2*fct.Diff2(lambda*u);
    end

Fun.Diff2 = @Diff2;

    function Precomp(u)
        fct.Precomp(lambda*u);
    end

Fun.Precomp = @Precomp;

    function f = Eval_Precomp
        f = fct.Eval_Precomp();
    end

Fun.Eval_Precomp = @Eval_Precomp;

    function d = Diff_Precomp
        d = lambda * fct.Diff_Precomp();
    end

Fun.Diff_Precomp = @Diff_Precomp;

    function d2 = Diff2_Precomp
        d2 = lambda2 * fct.Diff2_Precomp();
    end

Fun.Diff2_Precomp = @Diff2_Precomp;


end
