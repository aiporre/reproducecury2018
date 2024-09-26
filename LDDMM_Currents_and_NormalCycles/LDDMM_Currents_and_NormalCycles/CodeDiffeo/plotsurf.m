function h = plotsurf(S,clr,plotNormals)

if nargin<3
    plotNormals = 0;
end

if nargin<2
    clr = 'r';
end


if iscell(S)
    
    h = cell(1,length(S));
    for i=1:length(S)
        h{i} = plotsurf(S{i},clr,plotNormals);
    end
    
else
    
    if isfield(S,'Weights')
        S = rmfield(S,'Weights');
    end
    
    if isfield(S,'Label')
        S = rmfield(S,'Label');
    end
    
    N = [];
    if isfield(S,'Normals')
        N = S.Normals;
        S = rmfield(S,'Normals');
    end
    
    Sfields = fieldnames(S);
    S = rmfield(S,Sfields(~strcmp(Sfields,'Vertices')&~strcmp(Sfields,'Faces')));
    if size(S.Faces,2)==2
        h = patch(S,'EdgeColor',clr,'LineWidth',3);
    else
        h = patch(S,'FaceColor',clr,'EdgeColor','none');
    end
    axis off
    axis equal
    
    if plotNormals
        hold on
        if isempty(N)
            C = (S.Vertices(S.Faces(:,1),:)+S.Vertices(S.Faces(:,2),:)+S.Vertices(S.Faces(:,3),:))/3;
            N = .5*cross(S.Vertices(S.Faces(:,2),:)-S.Vertices(S.Faces(:,1),:),S.Vertices(S.Faces(:,3),:)-S.Vertices(S.Faces(:,1),:));
            %N = N./repmat(sqrt(sum(N.^2,2)),1,3);
            h = [h,quiver3(C(:,1),C(:,2),C(:,3),N(:,1),N(:,2),N(:,3),0,'Color',clr)];
        else
            h = [h,quiver3(S.Vertices(:,1),S.Vertices(:,2),S.Vertices(:,3),N(:,1),N(:,2),N(:,3),0,'Color',clr)];
        end
    end
    
end