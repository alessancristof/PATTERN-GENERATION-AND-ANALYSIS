function [Lmax,Lmin,Vmax1,Vmax2,Vmin1,Vmin2]=matrixDiagonalize2D(M11,M22,M12)
% diagonalize 2x2  real simmetric matrix

% INPUT
% M11, M22, and M12 are the element of the matrix
% M11, M22, and M12 can be matrices of any size (n x m x s)

% OUTPUT
% Lmax and Lmin: maximaland minimal eigenvalues, of the same size as M11 (n x m x s).
% Vmax1 and Vmax2 (n x m x s): first and second component of the normalized eigenvector associated
% to Lmax.
% Vmin1 and Vmin2 (n x m x s, optional): first and second component of the normalized eigenvector associated
% to Lmin.

T=M11+M22;
if nargout>2
    var=M12.*M12;
    mask=var<=eps;
    D=M11.*M22-var;
    clear var
else
    D=M11.*M22-M12.*M12;
end

var1=T/2;
var2=max(var1.*var1-D,0);
var2=sqrt(var2);

Lmax=var1+var2;
Lmin=var1-var2;

if nargout>2
    Vmax1=Lmax-M22;
    Vmax2=M12;
    Vmax=sqrt(Vmax1.^2+Vmax2.^2)+eps;
    Vmax1=Vmax1./Vmax;
    Vmax2=Vmax2./Vmax;
    Vmax1(mask)=1;
    Vmax2(mask)=0;
end

if nargout>4
    Vmin1=Lmin-M22;
    Vmin2=M12;
    Vmin=sqrt(Vmin1.^2+Vmin2.^2)+eps;
    Vmin1=Vmin1./Vmin;
    Vmin2=Vmin2./Vmin;
    Vmin1(mask)=0;
    Vmin2(mask)=1;
end