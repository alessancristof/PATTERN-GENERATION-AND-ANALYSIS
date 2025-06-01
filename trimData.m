function [xtrim,ytrim,xdis,ydis,x_me,y_me,x_sd,y_sd,x_se,y_se]=trimData(x,y,n,p,trimside)
% Progressive trimming of data (x, y) by binning. The data is ordered according to x,
%  binned in overlapped bins of size n, and each bin is trimmed based on a percentile p of y.
% Statistics of kept data in each bin are optionally computed.

% x (input): vector of data.
% y (input): second vector of data dependent on x (same size as x).
% n (input): number of values in each bin.
% p (input): [0 - 100] percentile of data to be removed in each bin
% trimside (input): type of trimming ('both' = p/2 percent of both lowest and highest
% values, 'lower' = p percent of lowest values, 'upper' = p percent of highest values).

% xtrim and ytrim (output): the trimmed data, i.e., that was kept at least in one bin
% xdis and ydis (output): the discarded data, i.e., that was never kept in any bins
% x_me y_me (optional output): nx1 vectors, the mean of x and y ever the kept data in each
% of the n bins.
% x_sd y_sd (optional output): nx1 vectors, the standard deviation of x and y over the kept data in each
% of the n bins.
% x_se y_se (optional output): nx1 vectors, the standard error of x and y over the kept data in each
% of the n bins.


if nargin<5
    trimside='both';
end

[x,ind]=sort(x);
y=y(ind);

n2=round(n/2);
mask=false(size(x));
for i=1:floor(numel(x)/n2)
    a=(i-1)*n2+1;
    b=a+n-1;
    if b>numel(y)
        break
    end
    ind=a:b;
    y2=y(ind);
    switch trimside
        case 'both'
            ylow=prctile(y2,p/2);
            yup=prctile(y2,100-p/2);
            mask2=and(y2>=ylow,y2<=yup);
        case 'lower'
            ylow=prctile(y2,p);
            mask2=y2>=ylow;
        case 'upper'
            yup=prctile(y2,100-p);
            mask2=y2<=yup;
    end
    mask(ind)=or(mask(ind),mask2);
    if nargout>4
        x2=x(ind);
        x2=x2(mask2);
        y2=y2(mask2);
        if numel(x2)<2
            x_me(i,1)=nan;
            x_sd(i,1)=nan;
            x_se(i,1)=nan;
            y_me(i,1)=nan;
            y_sd(i,1)=nan;
            y_se(i,1)=nan;
        else
            x_me(i,1)=mean(x2);
            x_sd(i,1)=std(x2);
            x_se(i,1)=x_sd(i)/sqrt(numel(x2));
            y_me(i,1)=mean(y2);
            y_sd(i,1)=std(y2);
            y_se(i,1)=y_sd(i)/sqrt(numel(y2));
        end
    end
end
xtrim=x(mask);
ytrim=y(mask);
xdis=x(not(mask));
ydis=y(not(mask));







