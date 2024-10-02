function writeFunction(f,fun,sigma)
if iscell(fun)
    fprintf(f,'Sum\n%d\n',length(fun));
    for i=1:length(fun)
        fprintf(f,'%f\n',fun{i}{1});
        writeFunction(f,fun{i}{2},sigma{i}) ;
    end
else
    fprintf(f,'%s',fun);
    if ~strcmp(fun,'Cubic')
        fprintf(f,',sigma=\n%f',sigma);
        fprintf(f,'\n');
    end
end
    