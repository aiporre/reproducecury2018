function WriteLbfgs_Dlib(s,f)
    
fprintf(f,'Optim_Dlib\n');
fprintf(f,'Method=\n');
fprintf(f,'LBFGS\n');

if ~isfield(s,'optim_M')
        s.optim_M = 10;
end

if ~isfield(s,'optim_breakratio')
        s.optim_breakratio = '1e-7';
end


fprintf(f,'M=\n%d\n',s.optim_M);
fprintf(f,'BreakRatio=\n%s\n',s.optim_breakratio);

