function [Function,coef] = readFunction(f)

coef = [];
Tag = fscanf(f,'%s',1);
if strcmp(Tag,'Sum')
    Nb = fscanf(f,'%d',1);
    Function = cell(Nb,1);
    coef = cell(Nb,1);
    for i=1:Nb
        Function{i}{1} = fscanf(f,'%f',1);
        [Function{i}{2},coef{i}] = readFunction(f);
    end
else
    if strcmp(Tag,'WeightedGaussian,coef=')
        Function = 'WeightedGaussian';
        coef = fscanf(f,'%f',1);
        fscanf(f,'%s',1);
    elseif strcmp(Tag,'SpecGaussian,coef0=')
        Function = 'SpecGaussian';
        coef(1) = fscanf(f,'%f',1);
        fscanf(f,'%s',1);
        coef(2) = fscanf(f,'%f',1);
        fscanf(f,'%s',1);
    else
        Function = strrep(Tag,',sigma=','');
        coef = fscanf(f,'%f',1);
    end
end
