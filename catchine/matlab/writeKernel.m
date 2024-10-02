function writeKernel(f,ker)

if strcmp(ker.Type,'FastGauss')
    if ~isfield(ker,'epsilon')
        warning('no required precision given for Fast Gauss Kernel; using default value 1e-3')
        ker.epsilon = 1e-3;
    end
    fprintf(f,'FastGauss,sigma,epsilon=\n%f %f\n',ker.Sigma,ker.epsilon);
elseif strcmp(ker.Type,'GridGauss')
    fprintf(f,'GridGauss,sigma,ratio=\n%f %f\n',ker.Sigma,ker.ratio);
elseif strcmp(ker.Type,'CauchyGpu')
    fprintf(f,'CauchyGpu,sigma=\n%f\n',ker.Sigma);
elseif strcmp(ker.Type,'GaussGpu')
    fprintf(f,'GaussGpu,sigma=\n%f\n',ker.Sigma);
elseif strcmp(ker.Type,'SumGaussGpu')
    fprintf(f,'SumGaussGpu,Weights=\n');
    writeArr(ker.Weights,f);
    fprintf(f,'Sigmas=\n');
    writeArr(ker.Sigmas,f);
elseif strcmp(ker.Type,'SumCauchyGpu')
    fprintf(f,'SumCauchyGpu,Weights=\n');
    writeArr(ker.Weights,f);
    fprintf(f,'Sigmas=\n');
    writeArr(ker.Sigmas,f);
elseif strcmp(ker.Type,'Tri')
    fprintf(f,'Tri,functionTilde=\n');
    if strcmp(ker.FunctionTilde,'WeightedGaussian')
        fprintf(f,'WeightedGaussian,coef=\n%f\nsigma=\n%f\n',ker.CoefTilde,ker.Sigma);
    else
        fprintf(f,'%s,sigma=\n%f\n',ker.FunctionTilde,ker.Sigma);
    end
    fprintf(f,'functionOrtho=\n');
    if strcmp(ker.FunctionOrtho,'WeightedGaussian')
        fprintf(f,'WeightedGaussian,coef=\n%f\nsigma=\n%f\n',ker.CoefOrtho,ker.Sigma);
    elseif strcmp(ker.FunctionOrtho,'SpecGaussian')
        fprintf(f,'SpecGaussian,coef0=\n%f\ncoef1=\n%f\nsigma=\n%f\n',ker.CoefOrtho(1),ker.CoefOrtho(2),ker.Sigma);
    else
        fprintf(f,'%s,sigma=\n%f\n',ker.FunctionOrtho,ker.Sigma);
    end    
else
    fprintf(f,'%s,function=\n',ker.Type);
    writeFunction(f,ker.Function,ker.Sigma);
end

