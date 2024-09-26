function s = MatchLandmarks(x,y,sigma)

% landmark matching via LDDMM
% x and y must be n*d matrices giving coordinates of n points in R^d
% sigma is width of kernel. The kernel is gaussian by default, this can be
% changed in the body of this function.
% output is a structure s which encodes the optimal deformation found.

[n,d] = size(x);

% choose kernel ; here we choose a gaussian kernel with width sigma
ker = ScalarKernel(GaussFunction,sigma);    % direct implementation in Matlab
% ker = ScalarMexKernel('Gauss',sigma);     % using Mex files to speed up
% ker = ScalarMexGpuKernel('Gauss',sigma);  % using Mex files and Gpu

% initialize positions and momenta
q0 = x;
p0 = zeros(n,d);

% gammaR parameter (0 for exact matching)
gammaR = 0;

% initialize matching
H = HamiltSys(ker,n,d);
DM = DiffeoMatch(H,gammaR);

OD = OptimDiffeo(DM);

OD.SetOptimOption('Display','iter');
OD.SetOptimOption('HessUpdate','lbfgs');

OD.SetOptimOption('TolFun',1e-4);
OD.SetOptimOption('TolX',1e-8);

method = 1;
switch method
    case 1
        % optimize on initial momenta only
        DM.AddTarget(TargetLandmarks(x,y));
        p0 = OD.Optimize_onlymom(p0,q0);
    case 2
        % optimize on initial positions and momenta (symetric formulation)
        DM.AddTarget(TargetLandmarks(x,x),0);
        DM.AddTarget(TargetLandmarks(x,y),1);
        [p0,q0] = OD.Optimize(p0,q0);
end

s.H = H;
s.p0 = p0;
s.q0 = q0;
s.dist = sqrt(H.Energy(p0,q0,1));






