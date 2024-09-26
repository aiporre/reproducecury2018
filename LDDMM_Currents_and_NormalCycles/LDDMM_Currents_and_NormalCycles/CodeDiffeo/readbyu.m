function [V,F] = readbyu(byufile,reducefactor)

if nargin<2
    reducefactor = 1;
end

% load vertices of byu surface

fbyu = fopen(byufile,'r');

% read header
ncomponents = fscanf(fbyu,'%d',1);	% number of components
npoints = fscanf(fbyu,'%d',1);		% number of vertices
nfaces = fscanf(fbyu,'%d',1);		% number of faces
nedges = fscanf(fbyu,'%d',1);		% number of edges
fscanf(fbyu,'%d',2*ncomponents);	% components (ignored)

% read data
V = fscanf(fbyu,'%f',[3,npoints]);		% vertices
F = fscanf(fbyu,'%d',[3,nfaces]);		% faces

fclose(fbyu);

ind = [find(F<0);nedges+1];
F = abs(F);

if reducefactor~=1
    [F,V] = reducepatch(F',V',reducefactor);
    F = F';
    V = V';
end

if nargout==1
    S.Vertices = V';
    S.Faces = F';
    V = S;
end
