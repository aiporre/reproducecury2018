
clear s target

s.T = 10;
%s.x = [1,0,0;0,1,0;0,0,1];%rand(3,3);
%[s.x,target{1}.vx] = readbyu('001_lpt_r.byu');
%[target{1}.y,target{1}.vy] = readbyu('000_lpt_r.byu');

Nx = 200;
%s.x = [rand(3,Nx),[0,1,0,1;0,0,1,1;0,0,0,0]];
s.x = 2*(rand(3,Nx)-.5);
s.x = s.x(:,sum(s.x(1:2,:).^2)<1);
s.x(3,:) = 0*s.x(1,:);
target{1}.vx = delaunay(s.x(1,:),s.x(2,:))';
%target{1}.vx = orient(target{1}.vx);

Ny = 200;
%target{1}.y = [rand(3,Ny),[0,1,0,1;0,0,1,1;0,0,0,0]];
target{1}.y = 2*(rand(3,Ny)-.5);
target{1}.y = target{1}.y(:,sum(target{1}.y(1:2,:).^2)<1)
target{1}.y(3,:) = .5 + sum((target{1}.y(1:2,:)).^2);
target{1}.y(1,:) = 2*target{1}.y(1,:);
target{1}.y(2,:) = .5*target{1}.y(2,:);
target{1}.vy = delaunay(target{1}.y(1,:),target{1}.y(2,:))';
%s.vy = orient(s.vy);

nfx = size(target{1}.vx,2);
nfy = size(target{1}.vy,2);

s.useDef = 'LargeDef';
s.targetweights = 1;
s.gammaR = 0;
s.optim_maxiter = 5000;
s.optim_stepsize = 1;
s.optim_verbosemode = 2;
s.optim_breakratio = 1e-6;
s.optim_loopbreak = 10;
s.rigidmatching = 0;
s.useoptim = 'adaptdesc';


target{1}.method = 'surfcurr';

s.sigmaV = .5;

% s.CppKer.Type = 'SqDistScalar';
% s.CppKer.Function = 'Cauchy';

% % curl free kernel, ancien style
% target{1}.CppKer.Type = 'CurlFree';
% target{1}.CppKer.Function = 'Cauchy2';

% curl free kernel, TRI style
% TRI kernels with gaussian functions
% 2 cas envisagés :
% 1) example 1 du papier :
%   k_tilde(r) = a*exp(-c*r^2)  et  k_ortho(r) = (b-a*r^2)*exp(-c*r^2)
%   avec a >= 0 et b >= (d-1)*a/(2*c) avec d=2 ici
%   a=0 donne un noyau scalaire et b=(d-1)*a/(2*c) donne un noyau div-free
% 1) example 2 du papier, mais le coef a ici est l'opposé de celui du papier :
%   k_tilde(r) = a*exp(-c*r^2)  et  k_ortho(r) = b*exp(-c*r^2)
%   avec a <= 0 et b >= -a/(2*c)
%   a=0 donne un noyau scalaire et b=-a/(2*c) donne un noyau curl-free
c = 1/s.sigmaV^2;
a = -1;
b = -a/(2*c);
s.CppKer.Type = 'Tri';
s.CppKer.FunctionTilde = 'WeightedGaussian';
if a<0
    s.CppKer.CoefTilde = a;
    s.CppKer.FunctionOrtho = 'WeightedGaussian';
    s.CppKer.CoefOrtho = b;
else
    s.CppKer.CoefTilde = a;
    s.CppKer.FunctionOrtho = 'SpecGaussian';
    s.CppKer.CoefOrtho(1) = b;
    s.CppKer.CoefOrtho(2) = -a;
end



target{1}.sigmaW = .5;
%target{1}.y = s.x+.2;%rand(3,3);
%target{1}.vx = [1;2;3];
%target{1}.vy = [1;2;3];
target{1}.wx = ones(1,nfx);
target{1}.wy = ones(1,nfy);

s = matchCpp(s,target);

clf
s.showpoints = 0;
affiche(s);

s.transmatrix = eye(3);
s.transvector = zeros(3,1);

s.showtraj = 1;

%s.showgrid = 1;
makewrl('ess.wrl',s);
