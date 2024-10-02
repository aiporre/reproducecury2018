function writeArr(A,f,realsize)

if nargin < 3
    realsize = 1;
end
sz = size(A);
if realsize
    if length(sz)==2
        if sz(2)==1
            sz = sz(1);
        elseif sz(1)==1
            sz = sz(2);
        end
    end
end
fprintf(f,'%d \n',length(sz));
fprintf(f,'%d ',sz);
fprintf(f,'\n');
fprintf(f,'%f ',A);
fprintf(f,'\n');

% function writeArr(A,f,groupes)
% 
% if nargin==2
%     groupes = ones(1,ndims(A));
% end
% 
% index = 1;
% for i = 1:length(groupes)
%     A = permute(A,index+groupe(i)-1:-1:index);
%     index = index + groupe(i);
% end
% 
% fprintf(f,'%d \n',length(size(A)));
% fprintf(f,'%d ',size(A));
% fprintf(f,'\n');
% fprintf(f,'%f ',A);
% fprintf(f,'\n');