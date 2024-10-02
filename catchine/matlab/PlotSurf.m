function h=PlotSurf(V,F,clr)

if isstruct(V) && nargin ==2
    tmp = V;
    clr = F;
    clear V;
    V = tmp.Vertices;
    F = tmp.Faces;
end

if iscell(V)
    [V F] = ConcatenateSurf(V,F);
end

S.Vertices = V';
S.Faces = F';
h = patch(S,'FaceColor',clr,'EdgeColor','none');
view(-98,12),axis equal, camlight, alpha .2,