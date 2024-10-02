function writeLargeDef_InitParam(s,f)

fprintf(f,'LargeDef_InitParam\n');
fprintf(f,['Dimension=\n',num2str(size(s.x,1)),'\n']);
if ~isfield(s,'T')
    s.T = 30;
end
fprintf(f,'NumTimeSteps=\n%d\n',s.T);
fprintf(f,'NumPoints=\n%d\n',size(s.x,2));
if isfield(s,'sigmaV')
    if length(s.sigmaV) > 1
        error('varying sigmaV is not allowed with LargeDef_InitParam')
    else
        s.CppKer.Sigma = s.sigmaV;
    end
end
fprintf(f,'Solver=\n');
writeSolver(f,s);
fprintf(f,'Kernel=\n');
writeKernel(f,s.CppKer);
fprintf(f,'InitPos,InitMom=\n');
writeArr(s.x,f,0);
if ~isfield(s,'mom')
    s.mom = zeros(size(s.x));
end
writeArr(s.mom(:,:,1),f,0);
