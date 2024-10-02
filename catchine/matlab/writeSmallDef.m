function writeSmallDef(s,f)

fprintf(f,'SmallDef\n');
fprintf(f,['Dimension=\n',num2str(size(s.x,1)),'\n']);
fprintf(f,'NumPoints=\n%d\n',size(s.x,2));
fprintf(f,'Kernel=\n');
if isfield(s,'sigmaV')
    if length(s.sigmaV)>1
        disp('s.sigmaV should be of length 1 for small deformations; keeping only first value')
    end
    s.CppKer.Sigma = s.sigmaV(1);
end
writeKernel(f,s.CppKer);
if isfield(s,'mom')
    if isfield(s,'phix')
        fprintf(f,'Position,Momentum,Phi=\n');
        writeArr(s.x,f);
        writeArr(s.mom,f);
        writeArr(s.phix,f);
    else
        fprintf(f,'Position,Momentum=\n');
        writeArr(s.x,f);
        writeArr(s.mom,f);
    end
else
    fprintf(f,'Position=\n');
    writeArr(s.x,f);
end
