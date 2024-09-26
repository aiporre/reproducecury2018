
% This matlab scripts performs surface matching via currents on two hippocampal surfaces
% using GPU accelaration for kernel convolutions.

addpath CodeDiffeo

%% load data
clear S T

S = readbyu('TestData/s1002_hippoL.byu',1);
T = readbyu('TestData/s1007_hippoL.byu',1);

%% parameters
sV = 10;
sW = [20,5,1];

%% matching
[~,phiS] = MatchSurf(S,T,sV,sW);

%% display
load
figure
hold on
plotsurf(S,'b');
plotsurf(T,'r');
plotsurf(phiS,'g');
camlight
axis equal
alpha .5

