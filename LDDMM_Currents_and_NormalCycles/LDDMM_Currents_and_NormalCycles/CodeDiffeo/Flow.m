function [phiz,Z] = Flow(phi,z,ti,tf)

% This function computes the image of any given set of points zi by the LDDMM 
% deformation encoded in structure phi.
% z must be a m*d matrix giving coordinates of m points zi in R^d
% phi must be the ouput of a call to MatchLandmarks or MatchCurves function
% output phiz contains coordinates of the phi(zi)
% optional output Z contains the full trajectories of the points zi through
% the flow
% optional ti and tf scalar inputs give initial and final time of
% integration. Default to ti=0 and tf=1.

if nargin<3
    ti = 0;
end
if nargin<4
    tf = 1;
end

[phiz,Z] = phi.H.Flow(z,phi.p0,phi.q0,ti,tf);

        

