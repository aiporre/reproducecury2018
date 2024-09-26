function [phi,phiS] = MatchCurves_NC(S,T,sigmaV,sigmaW,p0)

% curve matching in R^d via LDDMM and currents
% S and T are source and target curves, given either as :
% - nvS*d and nvT*d arrays of coordinates in R^d of the two sequences of 
%   points forming the curves (curve S has nvS points and curve T has nvT points) 
% - two structures, each with fields Vertices (nv*d array) and Faces (nf*2 array)
%   giving vertices coordinates and face connectivity of segments in R^d
% sigmaV is width of deformation kernel. 
% sigmaW is width of matching kernel. 
%   These two kernels are gaussian by default, this can be
%   changed in the body of this function.
% p0 is optional nv*d array which gives initialization of momentum vectors (default is zero values)
% output is a structure s which encodes the optimal deformation found.

if ~isstruct(S)
    S = curve2mesh(S);
    T = curve2mesh(T);
end

[n,d] = size(S.Vertices);

if d==2
    S.Vertices(:,3) = 0;
    T.Vertices(:,3) = 0;
    d = 3;
end

% initialize positions and momenta
q0 = S.Vertices;
if nargin < 5
   p0 = zeros(n,d);
end

for k=1:length(sigmaW)
    
    % choose deformation kernel ; here we choose a sum of 4 Cauchy kernels
    % with width sigmaV
    kerV = ScalarMexKernel('Sum4Cauchy',sigmaV); 
    
    % gammaR parameter
    gammaR = 0.01;
    
    % initialize matching
    H = HamiltSys(kerV,n,d);
    DM = DiffeoMatch(H,gammaR);
    
    OD = OptimDiffeo(DM);
    
    %OD.SetOptimOption('Display','iter');
    %OD.SetOptimOption('HessUpdate','lbfgs');
    
    OD.SetOptimOption('TolFun',1e-4);
    OD.SetOptimOption('TolX',1e-4);
    
    DM.AddTarget(TargetNormalCycleCurves_Constant_linear(S,T,sigmaW(k),100));
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



