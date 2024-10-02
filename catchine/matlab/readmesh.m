function [V,F,T] = readmesh(file)

%intwarning('on')

f = fopen(file,'r');
while ~strcmp(fgetl(f),'Dimension')
end
dim = fscanf(f,'%d',1);
while ~strcmp(fgetl(f),'Vertices')
end
nV = fscanf(f,'%d',1);
V = fscanf(f,'%f',[4,nV]);
T = logical(V(4,:));
V = V(1:3,:);
while ~strcmp(fgetl(f),'Triangles')
end
nF = fscanf(f,'%d',1);
F = int32(fscanf(f,'%d',[4,nF]));
F = F(1:3,:);
% while ~strcmp(fgetl(f),'Corners')
% end
% nC = fscanf(f,'%d',1);
% C = fscanf(f,'%f',[1,nF]);
% while ~strcmp(fgetl(f),'Edges')
% end
% nE = fscanf(f,'%d',1);
% E = fscanf(f,'%f',[3,nE]);
% while ~strcmp(fgetl(f),'Ridges')
% end
% nR = fscanf(f,'%d',1);
% R = fscanf(f,'%f',[1,nR]);
% while ~strcmp(fgetl(f),'Normals')
% end
% nN = fscanf(f,'%d',1);
% N = fscanf(f,'%f',[3,nN]);
% while ~strcmp(fgetl(f),'NormalAtVertices')
% end
% nNV = fscanf(f,'%d',1);
% NV = fscanf(f,'%f',[2,nNV]);
% while ~strcmp(fgetl(f),'NormalAtTriangleVertices')
% end
% nNFV = fscanf(f,'%d',1);
% NFV = fscanf(f,'%f',[3,nNFV]);
% while ~strcmp(fgetl(f),'Tangents')
% end
% nT = fscanf(f,'%d',1);
% T = fscanf(f,'%f',[3,nT]);
% 
% fgetl(f)
% fgetl(f)
% fgetl(f)
% fgetl(f)
% fgetl(f)
% fgetl(f)
% fgetl(f)
% fgetl(f)
% fgetl(f)
% fgetl(f)

fclose(f);
