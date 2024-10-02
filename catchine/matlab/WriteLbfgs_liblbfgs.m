function WriteLbfgs_liblbfgs(s,f)
    
fprintf(f,'Lbfgs_liblbfgs\n');

if ~isfield(s,'optim_M')
        s.optim_M = 10;
end

if ~isfield(s,'optim_linesearch')
        s.optim_linesearch = 'Default';
end


fprintf(f,'M=\n%d\n',s.optim_M);
fprintf(f,'LineSearch=\n%s\n',s.optim_linesearch);

