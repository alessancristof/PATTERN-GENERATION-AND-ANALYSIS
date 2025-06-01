function [G1,G2,G1var,G2var,gmax]=firGradient2D(V,DepthDimension,GradientDiscretization,Variance)

V=double(V);

switch GradientDiscretization
    case 'Roberts'
        p1=1;
        h1=[-1 0; 0 1]/sqrt(2);
    case 'Prewitt'
        p1=[1; 1]/2;
        h1=[-1; 1];
    case 'Sobel Tap 3'
        p1=[1; 2; 1]/4;
        h1=[-1; 0; 1]/2;
    case 'Sobel Tap 5'
        p1=1;
        h1=[-5 -8 -10 -8 -5; -4 -10 -20 -10 -4; 0 0 0 0 0; -4 -10 -20 -10 -4; -5 -8 -10 -8 -5]/240;
    case 'Gaussian Tap 5'
        sigma=0.96;
        x=[-2; -1; 0; 1; 2];
        p1=exp(-(x.^2)/(2*sigma^2));
        p1=p1/sum(p1);
        h1=(x/sigma^2).*p1;
    case 'Gaussian Tap 7'
        sigma=1.12;
        x=[-3; -2; -1; 0; 1; 2; 3];
        p1=exp(-(x.^2)/(2*sigma^2));
        p1=p1/sum(p1);
        h1=(x/sigma^2).*p1;
    case 'Sharr Tap 3'
        p1=[3; 10; 3]/16;
        h1=[-1; 0; 1]/2;
    case 'Sharr Tap 5'
        p1=[0.0231; 0.2413; 0.4713; 0.2413; 0.0231];
        h1=[-0.0831; -0.3339; 0; 0.3339; 0.0831];
    case 'Farid Tap 3' % av error < 2 degrees, err max 5 degrees
        p1=[0.229879; 0.540242; 0.229879];
        h1=[-0.425287; 0; 0.425287];
    case 'Farid Tap 4'
        p1=[0.092645; 0.407355; 0.407355; 0.092645];
        h1=[-0.236506; -0.267576; 0.267576; 0.236506];
    case 'Farid Tap 5' % av error < 0.25 degrees, err max 0.4 degrees
        p1=[0.037659; 0.249153; 0.426375; 0.249153; 0.037659];
        h1=[-0.109604; -0.276691; 0; 0.276691; 0.109604];
end



p2=p1';
h2=h1';
h1f=h1*p2;
h2f=p1*h2;

switch DepthDimension
    case 1
        p1=permute(p1,[3 1 2]);
        p2=permute(p2,[3 1 2]);
        h1=permute(h1,[3 1 2]);
        h2=permute(h2,[3 1 2]);
        h1f=permute(h1f,[3 1 2]);
        h2f=permute(h2f,[3 1 2]);
    case 2
        p1=permute(p1,[1 3 2]);
        p2=permute(p2,[1 3 2]);
        h1=permute(h1,[1 3 2]);
        h2=permute(h2,[1 3 2]);
        h1f=permute(h1f,[1 3 2]);
        h2f=permute(h2f,[1 3 2]);
    case 3
       
end

G1=imfilter(V,h1,'replicate');
G1=imfilter(G1,p2,'replicate');

G2=imfilter(V,h2,'replicate');
G2=imfilter(G2,p1,'replicate');


if isempty(Variance)
    Gvar=[];
    gmax=[];
else
    Gsq=G1.^2+G2.^2;
    gmax=sqrt(max(Gsq,[],"all"));
    G1var=imfilter(Variance,h1.^2,'replicate');
    G1var=imfilter(G1var,p2.^2,'replicate');
    % G1var=imfilter(Variance,h1f.^2,'replicate');
    G2var=imfilter(Variance,h2.^2,'replicate');
    G2var=imfilter(G2var,p1.^2,'replicate');
    % G2var=imfilter(Variance,h2f.^2,'replicate');
    

    % c=(Gsq-Gvar)./(Gsq+eps);
    % c=sqrt(max(c,0));
    % G1=G1.*c;
    % G2=G2.*c;
end

