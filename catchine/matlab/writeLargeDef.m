function writeLargeDef(s,f)

fprintf(f,'LargeDef\n'); 
fprintf(f,['Dimension=\n',num2str(size(s.x,1)),'\n']);
if ~isfield(s,'T')
    s.T = 30;
end
fprintf(f,'NumTimeSteps=\n%d\n',s.T); 
fprintf(f,'NumPoints=\n%d\n',size(s.x,2)); 
if isfield(s,'sigmaV')
    if length(s.sigmaV) > 1
        if length(s.CppKer)==1
            tmp = s.CppKer;
            s = rmfield(s,'CppKer');
            for k=1:length(s.sigmaV)
                s.CppKer{k} = tmp;
            end
        end
        for k=1:length(s.sigmaV)
            s.CppKer{k}.Sigma = s.sigmaV(k);
        end
    else
        s.CppKer.Sigma = s.sigmaV;
    end
end
if length(s.CppKer)==1
        fprintf(f,'Kernel=same\n');
        writeKernel(f,s.CppKer);
else
        fprintf(f,'Kernel=\n');
        for t=1:s.T
            writeKernel(f,s.CppKer{t});
        end
end
if isfield(s,'X') && isfield(s,'mom') && size(s.X,3)==s.T && size(s.mom,3)==s.T
    fprintf(f,'Position,Momentum=\n');
    writeArr(s.X,f);
    writeArr(s.mom,f);
elseif isfield(s,'x') && isfield(s,'mom') && (size(s.mom,3)==1 || sum(sum(sum(abs(s.mom(:,:,2:end))))))
    fprintf(f,'Position(1),Momentum(1)=\n');
    writeArr(s.x,f);
    writeArr(s.mom,f);
else
    fprintf(f,'Position(1)=\n');
    writeArr(s.x,f);
end
