
% matching two curves via scalar measure model and large deformations

clear

Nx = 50;
x = zeros(2,Nx);
x(1,:) = (1:Nx)/(Nx);
x(2,:) = 0*sin(2*pi*(1:Nx)/(Nx));
s.x = x;
s.T = 10;
% s.X = repmat(s.x,[1 1 s.T]);
% s.mom = zeros(size(s.X));

Ny = 50;
y = zeros(2,Ny);
y(1,:) = (1:Ny)/(Ny);
y(2,:) = .5*sin(2*pi*(1:Ny)/(Ny));
target{1}.y = y;

s.gammaR = 0;
s.sigmaV = .25;
s.optim_maxiter = 100;
s.optim_stepsize = 1;
s.optim_verbosemode = 1;
s.optim_breakratio = 1e-10;
s.optim_loopbreak = 10;
s.rigidmatching = 0;
s.useoptim = 'adaptdesc';


target{1}.method = 'measures';
target{1}.vx = 1:Nx;

s.targetweights = [1];

s.useDef = 'LargeDef';
s.CppKer.Type = 'SqDistScalar';
s.CppKer.Function = 'Gaussian';
target{1}.CppKer.Type = 'SqDistScalar';
target{1}.CppKer.Function = 'Gaussian';
target{1}.CppKer.Sigma = .25;
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
s.xmarker = ' ';
s.ymarker = ' ';
s.phimarker = ' ';
affiche(s);
axis equal

s.target{1}.vx = [1:size(s.x,2)-1;2:size(s.x,2)];
s.target{1}.vy = [1:size(s.target{1}.y,2)-1;2:size(s.target{1}.y,2)];

