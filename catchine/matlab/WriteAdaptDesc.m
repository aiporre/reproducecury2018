function WriteAdaptDesc(s,f)
    
fprintf(f,'AdaptGradDescent\n');

if ~isfield(s,'optim_maxiter')
        s.optim_maxiter = 500;
end

if ~isfield(s,'optim_stepsize')
        s.optim_stepsize = 0.1;
        disp('no auto mode for gradient step size in C++ code. Setting it to 1')
end

if ~isfield(s,'optim_breakratio')
        s.optim_breakratio = 1e-6;
    end
    if ~isfield(s,'optim_loopbreak')
        s.optim_loopbreak = 40;
    end

fprintf(f,'Niters=\n%d\n',s.optim_maxiter);
fprintf(f,'StepSize=\n%f\n',s.optim_stepsize);
fprintf(f,'BreakRatio=\n%d\n',s.optim_breakratio);
fprintf(f,'LoopBreak=\n%f\n',s.optim_loopbreak);

