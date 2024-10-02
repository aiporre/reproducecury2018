
% erreurs sur les convolutions par réarrangement aléatoire

m = 1;
ns = 10.^(2:7);
d = 3;
k = 3;
sigma = 0.25;
types = {@single,@double};
tagtypes = {'single','double'};

for t = 1:length(types)
    type = types{t};
    for k = 1:length(ns)
        n = ns(k);
        
        
        x = type(rand(m,d));
        y =  type(rand(n,d));
        beta = type(rand(n,k));
        
        Kxy = type(0);
        for i = 1:d
            Kxy = Kxy + (repmat(x(:,i),1,n)-repmat(y(:,i)',m,1)).^2;
        end
        Kxy = 1./(1+Kxy/sigma^2);
        
        gamma1 = Kxy*beta;
        ind = randperm(n);
        gamma2 = Kxy(:,ind)*beta(ind,:);
        err = sum(sqrt(sum((gamma1-gamma2).^2,2))) / sum(sqrt(sum((gamma1).^2,2)));
        
        disp(['m=',num2str(m),', n=',num2str(n),', type=',tagtypes{t},', relative error : ',num2str(err)])
    end
end


