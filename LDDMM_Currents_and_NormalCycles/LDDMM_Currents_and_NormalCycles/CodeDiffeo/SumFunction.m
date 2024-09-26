function Fun = SumFunction(fct)

n = length(fct);

    function f = Eval(u)
        f = 0;
        for k = 1:n
            f = f + fct{k}.Eval(u);
        end
    end

Fun.Eval = @Eval;

    function d = Diff(u)
        d = 0;
        for k = 1:n
            d = d + fct{k}.Diff(u);
        end
    end

Fun.Diff = @Diff;

    function d2 = Diff2(u)
        d2 = 0;
        for k = 1:n
            d2 = d2 + fct{k}.Diff2(u);
        end
    end

Fun.Diff2 = @Diff2;

    function Precomp(u)
        for k = 1:n
        	fct{k}.Precomp(u);
        end
    end

Fun.Precomp = @Precomp;

    function f = Eval_Precomp
        f = 0;
        for k = 1:n
            f = f + fct{k}.Eval_Precomp();
        end
    end

Fun.Eval_Precomp = @Eval_Precomp;

    function d = Diff_Precomp
        d = 0;
        for k = 1:n
            d = d + fct{k}.Diff_Precomp();
        end
    end

Fun.Diff_Precomp = @Diff_Precomp;

    function d2 = Diff2_Precomp
        d2 = 0;
        for k = 1:n
            d2 = d2 + fct{k}.Diff2_Precomp();
        end
    end

Fun.Diff2_Precomp = @Diff2_Precomp;


end
