function [tabcoeff,mmax,lmax] = SobolevCoeff(s,mmax,lmax)

tabcoeff = [];

for k = 0:mmax + lmax

    tabcoeff = [tabcoeff, (1 +k*(k+1))^s];

end