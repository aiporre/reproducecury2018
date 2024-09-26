function [ptf,qtf] = Shoot(phi,tf)

% This function computes the image of any given set of points zi by the LDDMM 
% deformation encoded in structure phi.
% phi must be the ouput of a call to MatchLandmarks or MatchCurves function
% output phiz contains coordinates of the phi(zi)
% optional tf scalar input gives final time of
% integration. Defaults to tf=1.
% ptf and qtf outputs give momentum vectors and positions at time tf

if nargin<2
    tf = 1;
end

[ptf,qtf] = phi.H.Shoot(phi.p0,phi.q0,0,tf);

        

