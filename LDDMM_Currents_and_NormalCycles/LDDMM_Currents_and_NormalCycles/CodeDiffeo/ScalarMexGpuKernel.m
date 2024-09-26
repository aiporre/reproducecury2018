function ker = ScalarMexGpuKernel(name,sigma,dimpoint,dimvect,weight,deviceID)

% scalar kernel using GPU Mex files

if nargin < 6
    deviceID = 0;
end

if nargin < 5
    weight = 1;
end

ker = ScalarMexKernel(name,sigma,weight,'Gpu',dimpoint,dimvect,deviceID);

