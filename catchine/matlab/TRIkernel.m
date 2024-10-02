function CppKer = TRIkernel(tag,sigma)

% TRI kernels with gaussian functions
% 2 cas envisagés :
% 1) example 1 du papier :
%   k_tilde(r) = a*exp(-c*r^2)  et  k_ortho(r) = (b-a*r^2)*exp(-c*r^2)
%   avec a >= 0 et b >= (d-1)*a/(2*c) avec d=dimension of ambient space
%   a=0 donne un noyau scalaire et b=(d-1)*a/(2*c) donne un noyau div-free
% 1) example 2 du papier, mais le coef a ici est l'opposé de celui du papier :
%   k_tilde(r) = a*exp(-c*r^2)  et  k_ortho(r) = b*exp(-c*r^2)
%   avec a <= 0 et b >= -a/(2*c)
%   a=0 donne un noyau scalaire et b=-a/(2*c) donne un noyau curl-free

c = 1/sigma^2;

switch tag
    case 'Scalar'
        a = 0;
        b = 1;
    case 'CurlFree'
        a = -1;
        b = -a/(2*c);
    case 'DivFree'
        d = 3;
        a = 1;
        b = (d-1)*a/(2*c);
end
CppKer.Type = 'Tri';
CppKer.FunctionTilde = 'WeightedGaussian';
if a<0
    CppKer.CoefTilde = a;
    CppKer.FunctionOrtho = 'WeightedGaussian';
    CppKer.CoefOrtho = b;
else
    CppKer.CoefTilde = a;
    CppKer.FunctionOrtho = 'SpecGaussian';
    CppKer.CoefOrtho(1) = b;
    CppKer.CoefOrtho(2) = -a;
end

