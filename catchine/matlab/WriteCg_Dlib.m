function WriteCg_Dlib(s,f)
    
fprintf(f,'Optim_Dlib\n');
fprintf(f,'Method=\n');
fprintf(f,'CG\n');

if ~isfield(s,'optim_breakratio')
        s.optim_breakratio = '1e-7';
end


fprintf(f,'BreakRatio=\n%s\n',s.optim_breakratio);

