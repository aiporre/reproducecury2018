% circle to ellipse via measure matching 

clear

Nx = 17;
x = zeros(3,Nx);
x(1,:) = cos(2*pi*(1:Nx)/(Nx));
x(2,:) = sin(2*pi*(1:Nx)/(Nx));

s.x = x;


Ny = 15;
y = zeros(3,Ny);
y(1,:) = 2*cos(2*pi*(1:Ny)/(Ny));
y(2,:) = sin(2*pi*(1:Ny)/(Ny));
target{1}.y = y;


s.gammaR = 0;
s.sigmaV = 1;
s.optim_maxiter = 500;
s.optim_stepsize = .01;
s.optim_verbosemode = 2;
s.optim_breakratio = 1e-10;
s.optim_loopbreak = 10;
s.rigidmatching = 0;
s.useoptim = 'adaptdesc';
%s.useoptim = 'fixedesc';
s.T = 30;


target{1}.method = 'measures';
%target{1}.method = 'landmarks';
target{1}.vx = 1:Nx;
target{1}.CppKer.Sigma = 2;
target{1}.CppKer.Type = 'SqDistScalar';
target{1}.CppKer.Function = 'Gaussian';

s.targetweights = [1];

%s.useDef = 'LargeDef';
%s.useDef = 'SmallDef';
s.CppKer.Type = 'SqDistScalar';
s.CppKer.Function = 'Gaussian';

s.useDef = 'LargeDef_InitParam';
s.useoptim = 'lbfgs_Quentin';

%s.mom = .5*s.x / 3.21;

s = matchCpp(s,target);



clf
s.showtraj = 1;
s.optim_verbosemode = 1;
s.transmatrix = eye(3);
s.transvector = ones(3,1);
s.showgrid = 1;
s.gridsize = 30;
s.usefgt = 0;
s.tau = 1/(s.T-1);
s.normcoefV = ones(1,s.T);
affiche(s);
axis equal

% s.show = {0};
% s.target{1}.vx = 1:size(s.x,2);
% s.target{1}.vy = 1:size(s.target{1}.y,2);
% makewrl('ess.wrl',s);





