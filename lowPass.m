function [I,ker]=lowPass(I,p,n,beta,padding)

q=1;
maxpq = max(p,q);
fc = 1/maxpq;
order = 2*round(n*maxpq);
b = fir1(order,fc,kaiser(order+1,beta));
% b=[0 0 0 0 0 0 0 0 0 0 0 0 0.5 1 0.5 0 0 0 0 0 0 0 0 0 0 0 0];
% b=[0 0 0 0 0 0 0 0 0 0 0 -0.1 0.25 0.7 0.25 -0.1 0 0 0 0 0 0 0 0 0 0 0];
b=b/sum(b);
ker=b'*b;

if isempty(I)
    figure
    var=1:numel(b);
    var=var-mean(var);
    plot(var,b,'Marker','.')
    var=round((51-numel(b))/2);
    if var>0
        ker=padarray(ker,[var var],0);
    end
    OTF = fftshift(psf2otf(ker));
    figure
    surf(abs(OTF));
    title('Corresponding |OTF|');
    axis square;
    axis tight
else

    I=imfilter(I,b,padding);
    I=imfilter(I,b',padding);
end



