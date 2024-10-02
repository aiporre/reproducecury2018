function Kxy = KernelMatrix(x,y,h,sigma)

[d,m] = size(x);
[d,n] = size(y);

Kxy = zeros(m,n);
for i=1:m
    for j=1:n
        Kxy(i,j) = h(sum((x(:,i)-y(:,j)).^2)/sigma^2);
    end
end

