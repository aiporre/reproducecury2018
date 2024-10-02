
% simple landmark matching

clear s

% % random
% nx = 10;
% Dim = 3;
% x = rand(Dim,nx);
% y = x+randn(Dim,nx)/10;

% % random
% nx = 2;
% Dim = 3;
% x = rand(2,nx);
% y = rand(2,nx);
% x(3,:) = 0;
% y(3,:) = 0;

% % given
% nx = 2;
% Dim = 3;
% x = [[0;0],[1;.1]];
% y = [[1;0],[0;.1]];
% x(3,:) = 0;
% y(3,:) = 0;


% % given
% nx = 3;
% Dim = 3;
% x = [0.106761861607241   0.779051723231275   0.890922504330789
%    0.653757348668560   0.715037078400694   0.334163052737496
%                    0                   0                   0];         
% y = [0.197809826685929   0.500022435590201   0.609866648422558
%    0.030540946304637   0.479922141146060   0.617666389588455
%                    0                   0                   0];

% given
nx = 5;
ech = linspace(0,1,nx);
x = [ech;zeros(2,nx)];
y = [ech;.25+(ech-.5).^2;zeros(1,nx)];

% % given (measures/currents)
% nx = 13;
% ech = linspace(0,1,nx);
% x = [ech;zeros(2,nx)];
% ny = 39;
% ech = linspace(0,1,ny);
% y = [ech;.25*(1+(cos(3*ech)).^2);zeros(1,ny)];

s.x = x;
target{1}.y = y;

if isfield(s,'mom')
    s.mom = s.mom(:,:,1);
end

s.gammaR = 0;
s.sigmaV = .5;

s.optim_maxiter = 10000;
s.optim_stepsize = .1;
s.optim_verbosemode = 1;
s.optim_breakratio = 1e-10;
s.optim_loopbreak = 20;
s.typefloat = 'float';
s.rigidmatching = 0;
%s.useoptim = 'adaptdesc';
s.useoptim = 'lbfgs_Quentin';
s.T = 20;

target{1}.method = 'landmarks';
target{1}.vx = 1:nx;

% target{1}.method = 'measures';
% target{1}.sigmaI = .5;
% target{1}.vx = 1:nx;

% target{1}.method = 'curvecurr';
% target{1}.sigmaW = .5;
% target{1}.vx = [1:nx-1;2:nx];
% target{1}.vy = [1:ny-1;2:ny];

s.targetweights = 1;

s.useDef = 'LargeDef_InitParam';
%s.useDef = 'SmallDef';
s.CppKer.Type = 'SqDistScalar';
%s.CppKer.Type = 'CauchyGpu';
s.CppKer.Function = 'Cauchy';
s = matchCpp(s,target);


clf
% hold on
% y = [target{1}.y,target{2}.y];
% plot3(y(1,:),y(2,:),y(3,:),'o')
% plot3(s.x(1,:),s.x(2,:),s.x(3,:),'*')
% plot3(squeeze(s.X(1,:,:))',squeeze(s.X(2,:,:))',squeeze(s.X(3,:,:))','LineWidth',3)
if strcmp(s.useDef,'SmallDef')
    s.T = 2;
    s.X(:,:,1) = s.x;
    s.X(:,:,2) = s.phix;
end
s.showtraj = 1;
s.optim_verbosemode = 2;
s.transmatrix = eye(3);
s.transvector = ones(3,1);
s.showgrid = 1;
s.gridsize = 50;
s.usefgt = 0;
s.tau = 1/(s.T-1);
s.normcoefV = ones(1,s.T);
affiche(s);
axis equal




