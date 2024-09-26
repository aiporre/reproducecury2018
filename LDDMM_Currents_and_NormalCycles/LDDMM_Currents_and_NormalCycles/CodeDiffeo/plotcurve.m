function h = plotcurve(C,clr)

if nargin<2
    clr = 'b';
end

if iscell(C)
    for k = 1:length(C)
        h{k} = plotcurve(C{k},clr);
    end
else
    d = size(C.Vertices,2);
    
    if d == 2
         h = plot([C.Vertices(C.Faces(:,1),1)';C.Vertices(C.Faces(:,2),1)'],[C.Vertices(C.Faces(:,1),2)';C.Vertices(C.Faces(:,2),2)'],clr,'LineWidth',3);
    else 
        
    h = plot3([C.Vertices(C.Faces(:,1),1)';C.Vertices(C.Faces(:,2),1)'],[C.Vertices(C.Faces(:,1),2)';C.Vertices(C.Faces(:,2),2)'],[C.Vertices(C.Faces(:,1),3)';C.Vertices(C.Faces(:,2),3)'],clr,'LineWidth',3);
    end
end
