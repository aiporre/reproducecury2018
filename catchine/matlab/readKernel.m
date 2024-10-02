function CppKer = readKernel(f)

CppKer.Type = strrep(fscanf(f,'%s',1),',function=','');

if strcmp(CppKer.Type,'FastGauss,sigma,order,K,e=')
    CppKer.Type = 'FastGauss';
    CppKer.Sigma = fscanf(f,'%f',1);
    CppKer.order = fscanf(f,'%d',1);
    CppKer.K = fscanf(f,'%d',1);
    CppKer.e = fscanf(f,'%f',1);
elseif strcmp(CppKer.Type,'FastGauss,sigma,epsilon=')
    CppKer.Type = 'FastGauss';
    CppKer.Sigma = fscanf(f,'%f',1);
    CppKer.epsilon = fscanf(f,'%f',1);
elseif strcmp(CppKer.Type,'CauchyGpu,sigma=')
    CppKer.Type = 'CauchyGpu';
   CppKer.Sigma = fscanf(f,'%f',1);
elseif strcmp(CppKer.Type,'GaussGpu,sigma=')
    CppKer.Type = 'GaussGpu';
    CppKer.Sigma = fscanf(f,'%f',1);
elseif strcmp(CppKer.Type,'SumGaussGpu,Weights=')
    CppKer.Type = 'SumGaussGpu';
    CppKer.Weights = readArr(f);
    fscanf(f,'%s',1);
    CppKer.Sigmas = readArr(f);
elseif strcmp(CppKer.Type,'SumCauchyGpu,Weights=')
    CppKer.Type = 'SumCauchyGpu';
    CppKer.Weights = readArr(f);
    fscanf(f,'%s',1);
    CppKer.Sigmas = readArr(f);
elseif strcmp(CppKer.Type,'Tri,functionTilde=')
    CppKer.Type = 'Tri';
    [CppKer.FunctionTilde,CppKer.CoefTilde] = readFunction(f);
    CppKer.Sigma = fscanf(f,'%f',1);
    fscanf(f,'%s',1);
    [CppKer.FunctionOrtho,CppKer.CoefOrtho] = readFunction(f);
    CppKer.Sigma = fscanf(f,'%f',1);
elseif strcmp(CppKer.Type,'SqDistScalar')
    [CppKer.Function,CppKer.Sigma] = readFunction(f);
end