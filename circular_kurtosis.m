function [exkurt,skew]=circular_kurtosis(v,k,w)

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
x2=cos(2*v);
y2=sin(2*v);

if isempty(w)
    xm=mean(x);
    ym=mean(y);
    xm2=mean(x2);
    ym2=mean(y2);
else
    var=sum(w);
    xm=sum(x.*w)/var;
    ym=sum(y.*w)/var;
    xm2=sum(x2.*w)/var;
    ym2=sum(y2.*w)/var;
end

me=atan2(ym,xm);
R=sqrt(xm^2+ym^2);
me2=atan2(ym2,xm2);
R2=sqrt(xm2^2+ym2^2);

skew=(R2*sin(me2-2*me))/(1-R)^(3/2);
exkurt=(R2*cos(me2-2*me)-R^4)/(1-R)^2;



