
% testing accuracy vs time for several time steps for shooting
% shooting with 100 random points and random momentum


clear

% s.x = [[0;0;0],[1;.1;0]];
% s.mom = [[1;0;0],[-1;0;0]];

Nx = 100;
s.x = rand(3,Nx);
s.mom = .1*randn(3,Nx);

s.sigmaV = 1;

s.CppKer.Type = 'SqDistScalar';
s.CppKer.Function = 'Gaussian';
Solvers = {'Euler','EulerTrapezoidal','RungeKutta4'};

echt = [5:30,40:10:100,200,500,1000];
for i=1:length(echt)
    s.T = echt(i);
    for k=1:length(Solvers)
        s.Solver = Solvers{k};
        tic
        stmp = shootingCpp(s);
        timeshoot{k}(i) = toc;
        phix{k}(:,:,i) = stmp.phix;
    end
end
for i=1:length(echt)
    for k=1:length(Solvers)
        diffnormphi{k}(i) = sqrt(sum(sum((phix{k}(:,:,i)-phix{length(Solvers)}(:,:,end)).^2)));
    end
end
clf
hold on
clr = 'rbgk';
for k=1:length(Solvers)
    h(k)=plot(timeshoot{k},diffnormphi{k},clr(k));
end
for i=1:length(echt)
    for k=1:length(Solvers)
        text(timeshoot{k}(i),diffnormphi{k}(i),num2str(echt(i)),'color',clr(k))
    end
end
legend(h,Solvers)
title('computing time (in s) vs error (as distance to most accurate) for shooting, for different methods and time steps');

% s.showtraj = 1;
% s.typefloat = 'double';
% s.optim_verbosemode = 1;
% s.transmatrix = eye(3);
% s.transvector = ones(3,1);
% s.showgrid = 1;
% s.gridsize = 1;
% s.usefgt = 0;
% s.sigmaV2 = s.sigmaV.^2;
% s.tau = 1/(s.T-1);
% s.normcoefV = ones(1,s.T);
% s.target = {};
% s.show = {0};



