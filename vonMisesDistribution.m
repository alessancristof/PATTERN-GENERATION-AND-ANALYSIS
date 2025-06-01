function Y=vonMisesDistribution(X,mu,sigma)


var2=cos(X-mu)/(sigma^2);
var2=var2-max(var2,[],"all");
Y=exp(var2);
Y=Y/sum(Y,"all");
