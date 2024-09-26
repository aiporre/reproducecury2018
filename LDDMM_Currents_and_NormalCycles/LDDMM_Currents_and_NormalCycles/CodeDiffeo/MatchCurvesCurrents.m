function [phi,phiS] = MatchCurves(S,T,sigmaV,sigmaW)

% curve matching in R^d via LDDMM and currents
% S and T are source and target curves, given either as :
% - nvS*d and nvT*d arrays of coordinates in R^d of the two sequences of 
%   points forming the curves (curve S has nvS points and curve T has nvT points) 
% - two structures, each with fields Vertices (nv*d array) and Faces (nf*2 array)
%   giving vertices coordinates and face connectivity of segments in R^d
% sigmaV is width of deformation kernel. 
% sigmaW is width of metching kernel. 
%   These two kernels are gaussian by default, this can be
%   changed in the body of this function.
% output is a structure s which encodes the optimal deformation found.

if ~isstruct(S)
    S = curve2mesh(S);
    T = curve2mesh(T);
end

[n,d] = size(S.Vertices);

% initialize positions and momenta
q0 = S.Vertices;
p0 = zeros(n,d);

for k=1:length(sigmaW)
    
    % choose deformation kernel ; here we choose a gaussian kernel with width sigma
    kerV = ScalarKernel(GaussFunction,sigmaV);    % direct implementation in Matlab
    % kerV = ScalarMexKernel('Gauss',sigmaV);     % using Mex files to speed up
    % kerV = ScalarMexGpuKernel('Gauss',sigmaV);  % using Mex files and Gpu
    
    % choose matching kernel ; here we choose a gaussian kernel with width sigma
    kerW = ScalarKernel(GaussFunction,sigmaW(k));    % direct implementation in Matlab
    % kerW = ScalarMexKernel('Gauss',sigmaW(k));     % using Mex files to speed up
    % kerW = ScalarMexGpuKernel('Gauss',sigmaW(k));  % using Mex files and Gpu
    
    % gammaR parameter
    gammaR = 0.01;
    
    % initialize matching
    H = HamiltSys(kerV,n,d);
    DM = DiffeoMatch(H,gammaR);
    
    OD = OptimDiffeo(DM);
    
    OD.SetOptimOption('Display','iter');
    OD.SetOptimOption('HessUpdate','lbfgs');
    
    OD.SetOptimOption('TolFun',1e-4);
    OD.SetOptimOption('TolX',1e-8);
    
    DM.AddTarget(TargetCurveCurr(S,T,kerW));
    p0 = OD.Optimize_onlymom(p0,q0);
    
    phi.H = H;
    phi.p0 = p0;
    phi.q0 = q0;
    phi.dist = sqrt(H.Energy(p0,q0,1));
    
    if nargout==2
        phiS = S;
        phiS.Vertices = Flow(phi,S.Vertices);
    end
    
end


function S = curve2mesh(x)
n = size(x,1);
S.Vertices = x;
S.Faces = [1:n-1;2:n]';



