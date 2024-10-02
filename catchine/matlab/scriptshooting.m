
% shooting with 500 random points and random momentum

clear

Nx = 500;
s.x = rand(3,Nx);
s.mom = .5*(rand(3,Nx)-.5);

s.sigmaV = 1;
s.T = 30;

s.typefloat = 'float';


%s.CppKer.Type = 'SqDistScalar';
%s.CppKer.Function = 'Cauchy';

s.CppKer.Type = 'CauchyGpu';

%s.CppKer.Type = 'FastGauss';
%s.CppKer.e = 16;
%s.CppKer.K = 20;
%s.CppKer.order = 9;

s = shootingCpp(s);

s.showtraj = 1;
s.optim_verbosemode = 1;
s.transmatrix = eye(3);
s.transvector = ones(3,1);
s.typefloat = 'float';
s.showgrid = 1;
s.gridsize = 1;
s.usefgt = 0;
s.tau = 1/(s.T-1);
s.normcoefV = ones(1,s.T);
s.target = {};
s.show = {0};
%affiche(s);
plot3(squeeze(s.X(1,:,:))',squeeze(s.X(2,:,:))',squeeze(s.X(3,:,:))')
axis equal

s.show = {0};
s.target{1}.vx = 1:size(s.x,2);
%s.target{1}.vy = 1:size(s.target{1}.y,2);
makewrl('ess.wrl',s);



