function [mu_best, sigma_best, a_best, CHI]=vonMisesInterp(X,Y,T,mu_int,sigma_int,flag_offset)

if size(X,1)==1
    X=X';
    Y=Y';
end

Y=Y./sum(Y,1);
N=size(Y,1);
eps_fact=1/1000;
effective_eps=1/N*eps_fact;

if nargin<5
    mu=[-179:180]/360*2*pi;
    sigma=[1:360]/360*2*pi;
else
    mu=mu_int/T*2*pi;
    sigma=sigma_int/T*2*pi;
end

sigma_min=min(sigma);
mu_best=nan;
sigma_best=nan;
a_best=nan;
CHI=inf;
X=X/T*2*pi;
for i=1:numel(mu)
    for j=1:numel(sigma)
        if flag_offset
            y=[ones(size(X)) exp(cos(X-mu(i))/(sigma(j)^2))];
            a=y\Y;
        else
            y=vonMisesDistribution(X,mu(i),sigma(j));
            a=1;
        end

        var=y*a+effective_eps;

        d=sum(((Y-var).^2)./var);
        if d<CHI
            CHI=d;
            mu_best=mu(i);
            sigma_best=sigma(j);
            a_best=a;
        end
    end
end

mu_best=mu_best/(2*pi)*T;
sigma_best=sigma_best/(2*pi)*T;
CHI=CHI/sum(Y);










