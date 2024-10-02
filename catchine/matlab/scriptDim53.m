
% matching random landmarks in 3D with large deformations

clear

% random
N = 3;
Dim = 53;
x = rand(Dim,N);
y = rand(Dim,N);
%y = x+randn(size(x))/5;
%x(3,:) = 0;
%y(3,:) = 0;

% % given
% N = 3;
% Dim = 3;
% x = [0.106761861607241   0.779051723231275   0.890922504330789
%    0.653757348668560   0.715037078400694   0.334163052737496
%                    0                   0                   0];         
% y = [0.197809826685929   0.500022435590201   0.609866648422558
%    0.030540946304637   0.479922141146060   0.617666389588455
%                    0                   0                   0];

s.x = x;
target{1}.y = y;

s.gammaR = 0;%.001;
s.sigmaV = 5;

s.optim_maxiter = 10000;
s.optim_stepsize = 1;
s.optim_verbosemode = 0;
s.optim_breakratio = 1e-10;
s.optim_loopbreak = 20;
s.rigidmatching = 0;
s.useoptim = 'adaptdesc';
%s.useoptim = 'fixedesc';
s.T = 50;


target{1}.method = 'landmarks';
target{1}.vx = 1:N;

s.targetweights = 1;

s.useDef = 'LargeDef';
s.CppKer.Type = 'SqDistScalar';%'CauchyGpu';
s.CppKer.Function = 'Cauchy';



sout{1} = matchCpp(s,target);


s.useDef = 'LargeDef_InitParam';
s.useoptim = 'lbfgs_Quentin';

sout{2} = matchCpp(s,target);

s.mom = sout{1}.mom(:,:,1);

sout{3} = matchCpp(s,target);

s.mom = sout{2}.mom;
s.X = sout{2}.X;

s.useDef = 'LargeDef';
s.useoptim = 'adaptdesc';

sout{4} = matchCpp(s,target);

s.mom = sout{4}.mom(:,:,1);
s.useDef = 'LargeDef_InitParam';
s.useoptim = 'lbfgs_Quentin';

sout{5} = matchCpp(s,target);


clf
hold on
for k=1:length(sout)
    s=sout{k};
% hold on
% y = [target{1}.y,target{2}.y];
% plot3(y(1,:),y(2,:),y(3,:),'o')
% plot3(s.x(1,:),s.x(2,:),s.x(3,:),'*')
% plot3(squeeze(s.X(1,:,:))',squeeze(s.X(2,:,:))',squeeze(s.X(3,:,:))','LineWidth',3)
if s.useDef==2
    s.T = 2;
    s.X(:,:,1) = s.x;
    s.X(:,:,2) = s.phix;
end
s.showtraj = 1;
switch k
    case 1
        s.trajtag = 'b';
    case 2
        s.trajtag = 'r';
    case 3
        s.trajtag = 'g';
    case 4
        s.trajtag = 'k';
    case 5
        s.trajtag = 'c';
end
h(k)=plot3(s.x(1,1),s.x(2,1),s.x(3,1),s.trajtag);
s.optim_verbosemode = 1;
s.transmatrix = eye(3);
s.transvector = ones(3,1);
s.showgrid = 0;
s.gridsize = 10;
s.usefgt = 0;
s.tau = 1/(s.T-1);
s.normcoefV = ones(1,s.T);
affiche(s);
axis equal

end

legend(h,{'LargeDef','LargeDefInitParam','LargeDef puis LargeDefInitParam','LargeDefInitParam puis LargeDef','LargeDefInitParam puis LargeDef puis LargeDefInitParam'})

