function WriteLbfgs_liblbfgs_Quentin(s,f)
    
fprintf(f,'Lbfgs_liblbfgs_Quentin\n');

if ~isfield(s,'optim_M')
        s.optim_M = 25;
end

fprintf(f,'M=\n%d\n',s.optim_M);

