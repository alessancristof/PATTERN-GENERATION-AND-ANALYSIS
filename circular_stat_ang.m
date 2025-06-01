function [me,sd,med,sdmed]=circular_stat_ang(v,k,w)

% v vector with values to compute statistics
% modulus of circular statistics
% w vector with weights for statistivs

v=v(:);

if nargin<2
   k=[]; 
end
if nargin<3
   w=[]; 
end


if not(isempty(k))
   v=(v/k)*2*pi; 
end

x=cos(v);
y=sin(v);

if isempty(w)
    x=mean(x);
    y=mean(y);
else
    var=sum(w);
    x=sum(x.*w)/var;
    y=sum(y.*w)/var;
end

me=atan2(y,x);
R=sqrt(x^2+y^2);
sd=sqrt(max(-2*log(R),0));

if nargout>2
    xmed=median(x);
    ymed=median(y);
    med=atan2(ymed,xmed);
    R=sqrt(xmed^2+ymed^2);
    sdmed=sqrt(-2*log(R));
end


if not(isempty(k))
    me=me/(2*pi)*k;
    sd=sd/(2*pi)*k;
    if nargout>2
        med=med/(2*pi)*k;
        sdmed=sdmed/(2*pi)*k;
    end
end






