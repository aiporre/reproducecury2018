function WriteBfgs_Dlib(s,f)
    
fprintf(f,'Optim_Dlib\n');
fprintf(f,'Method=\n');
fprintf(f,'BFGS\n');

if ~isfield(s,'optim_breakratio')
        s.optim_breakratio = '1e-7';
end


fprintf(f,'BreakRatio=\n%s\n',s.optim_breakratio);

