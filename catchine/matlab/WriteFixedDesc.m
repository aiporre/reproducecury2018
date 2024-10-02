function WriteFixedDesc(s,f)
    
fprintf(f,'FixedGradDescent\n');

if ~isfield(s,'optim_maxiter')
        s.optim_maxiter = 500;
end

if ~isfield(s,'optim_stepsize')
        s.optim_stepsize = 1;
        disp('no auto mode for gradient step size in C++ code. Setting it to 1')
end
   
fprintf(f,'Niters=\n%d\n',s.optim_maxiter);
fprintf(f,'StepSize=\n%f\n',s.optim_stepsize);

