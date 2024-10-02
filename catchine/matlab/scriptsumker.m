% matching random landmarks in 2D with large deformations
% with sum of gauss kernels

clear s target

% given
nx = 15;
ech = linspace(0,1,nx);
x = [ech;zeros(1,nx)];
t = pi/3;
y = [cos(t),-sin(t);sin(t),-cos(t)]*[ech;.25+(ech-.5).^2];

s.x = x;
target{1}.y = y;

if isfield(s,'mom')
    s.mom = s.mom(:,:,1);
end

s.gammaR = 0;

s.optim_maxiter = 10000;
s.optim_stepsize = .1;
s.optim_verbosemode = 2;
s.optim_breakratio = 1e-10;
s.optim_loopbreak = 20;
s.typefloat = 'double';
s.rigidmatching = 0;
%s.useoptim = 'adaptdesc';
s.useoptim = 'lbfgs_Quentin';
s.T = 20;

target{1}.method = 'landmarks';
target{1}.vx = 1:nx;

s.targetweights = 1;

s.useDef = 'LargeDef_InitParam';

% div-free kernel
% s.CppKer.Sigma = 5;
% a = 1;
% c = 1/s.CppKer.Sigma^2;
% b = 1/(2*c);
% s.CppKer.Type = 'Tri';
% s.CppKer.FunctionTilde = 'WeightedGaussian';
% s.CppKer.CoefTilde = a;
% s.CppKer.FunctionOrtho = 'SpecGaussian';
% s.CppKer.CoefOrtho(1) = b;
% s.CppKer.CoefOrtho(2) = -a;


% sum of Gaussian kernels (CPU)
s.CppKer.Type = 'SqDistScalar';
s.CppKer.Function = {{1,'Gaussian'},{1,'Gaussian'},{1,'Gaussian'},{1,'Gaussian'},{1,'Gaussian'}};
s.CppKer.Sigma = {2 1 .5 .2 .1};

% sum of Gaussian kernels (GPU)
% s.CppKer.Type = 'SumGaussGpu';
% s.CppKer.Weights = [1 1];
% s.CppKer.Sigmas = [.5 .05];

s = matchCpp(s,target);


clf
if strcmp(s.useDef,'SmallDef')
    s.T = 2;
    s.X(:,:,1) = s.x;
    s.X(:,:,2) = s.phix;
end
s.showtraj = 1;
s.transmatrix = eye(3);
s.transvector = ones(3,1);
s.showgrid = 1;
s.gridsize = 50;
s.usefgt = 0;
s.tau = 1/(s.T-1);
s.normcoefV = ones(1,s.T);
affiche(s);
axis equal




