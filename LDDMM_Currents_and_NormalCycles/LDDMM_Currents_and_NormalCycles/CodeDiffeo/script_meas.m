
% Measure matching example : match two sets of points as Dirac measures

% first set of points : 15 points forming a circle
nx = 15; 
th = linspace(0,2*pi,nx+1)';
th = th(1:end-1);
x = .7*[cos(th),sin(th)];

% second set of points : 25 points forming an kind of deformed ellipse
ny = 25;
th = linspace(0,2*pi,ny+1)';
th = th(1:end-1);
mod = 1+sin(4*th)/8;
y = [2*cos(th).*mod,sin(th).*mod];

% set parameters
d = 2; % dimension
gammaR = 0.00001; % weight of the energy term

% kernel for LDDMM
sigmaV = .5;
%kerV = ScalarMexKernel('Gauss',sigmaV);
kerV = ScalarKernel(GaussFunction,sigmaV);

% kernel for measures
sigmaI = .5;
%kerI = ScalarMexKernel('Gauss',sigmaI);
kerI = ScalarKernel(GaussFunction,sigmaI);

% initialize positions and momenta
q0 = x;
p0 = zeros(nx,d);

% initialize LDDMM functions
H = HamiltSys(kerV,nx,d);
DM = DiffeoMatch(H,gammaR);
OD = OptimDiffeo(DM);

% initialize measure target
DM.AddTarget(TargetMeasures(x,y,kerI));

% perform the optimization
p0 = OD.Optimize_onlymom(p0,q0);

% plot results
figure(1)
clf
H.PlotGrid(p0,q0);
DM.Plot(p0,q0);
axis equal


