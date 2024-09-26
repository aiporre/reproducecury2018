
% This matlab scripts performs curve matching via normal cycles on two 
% fish contours

addpath CodeDiffeo

% load data
load('TestData/FishContours');

% main parameters (kernel widths)
sigmaV = .2; % deformation kernel
sigmaW = .5; % spatial kernel for Normal Cycles

% call the matching algorithm
[phi,phiC1] = MatchCurvesNormalCycles(C1,C2,sigmaV,sigmaW);

% display registration result
clf
hold on
plotcurve(C1,'b');
plotcurve(C2,'r');
plotcurve(phiC1,'g');
axis equal