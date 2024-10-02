function h = surfplot(S,clr,plotnormals)

if nargin<3
    plotnormals = 0;
end

if nargin<2
    clr = 'r'
end

if isfield(S,'Weights')
    S = rmfield(S,'Weights');
end

if size(S.Vertices,2)~=3
    S.Vertices = S.Vertices';
end

if size(S.Faces,2)~=3
    S.Faces = S.Faces';
end

h = patch(S,'FaceColor',clr,'EdgeColor','none');
view(140,0)
axis off
axis equal

if plotnormals
    hold on
    N = .5*cross(S.Vertices(S.Faces(:,2),:)-S.Vertices(S.Faces(:,1),:),S.Vertices(S.Faces(:,3),:)-S.Vertices(S.Faces(:,1),:));
    C = (S.Vertices(S.Faces(:,1),:)+S.Vertices(S.Faces(:,2),:)+S.Vertices(S.Faces(:,3),:))/3;
    h = [h,quiver3(C(:,1),C(:,2),C(:,3),N(:,1),N(:,2),N(:,3))];
end
