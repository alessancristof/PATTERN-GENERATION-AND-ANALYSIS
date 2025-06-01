function [theta_mean,theta_std,ShE,ShEspread,spread_mean,theta_exkurt,Histo]=circularHistogramStatistics(V1,V2,mask,flag_norm)

% circular statistics of non-directional vectors weighted by vector norm
% INPUT
% V1 and V2: first and second component of vectors defining direction,
% their norm define the weight used in statistical descriptors computation,
% they can be matrices of any size.
% mask: mask to select the vectors for analysis
% flag_norm: normalize the vectors discarding their norm

% OUTPUT
% theta_mean: angle circular mean (deg), vector n x 1 with n the number of
% element of V1.
% theta_std (n x 1): angle circular SD (deg)
% ShE: Shannon Entropy
% ShEspread (deg): exponential of Shannon entropy, it represent a conventional
% spread. 
% spread_mean (deg): measure of spread based on thresholding the histogram.
% theta_exkurt: circular kurtosis
% Histo: structure containing data and parameters relative to the histogram of direction angles, necessary for plotting.

if not(isempty(mask))
    V1=V1(mask);
    V2=V2(mask);
end
V1=V1(:);
V2=V2(:);

if isempty(V1)
    theta_exkurt=nan;
    theta_mean=nan;
    theta_std=nan;
    rho_m=nan;
    rho_sd=nan;
    ShE=nan;
    ShEspread=nan;
    spread_mean=nan;
    exclusion_percentage=nan;
    Histo=[];
    return
end

rho=sqrt(V1.^2+V2.^2);
rho_m=mean(rho);
rho_sd=std(rho);

% remove zero norm vectors
mask=rho>10*eps;
exclusion_percentage=100-sum(mask)/numel(mask)*100;
Histo.ExclusionPercentage=exclusion_percentage;

rho=rho(mask);
V1=V1(mask);
V2=V2(mask);


if flag_norm==1 % normalized
    V1=V1./rho;
    V2=V2./rho;
    rho=ones(size(rho));
end

V1(abs(V1)<=eps)=eps;
theta=atan(V2./V1);
theta=theta/pi*180;
[theta_mean,theta_std]=circular_stat_ang(theta,180,rho);
[theta_exkurt,theta_skew]=circular_kurtosis(theta,180,rho);

bins=[-89:1:90]-0.5;
Histo.Bins=bins;
h=double(w_hist(bins,theta,rho));
h=h/sum(h);
Histo.AngularDistribution=h;

% vonMises interpolation
mu_int=round(theta_mean)+(-10:0.1:10);
mu_int(mu_int<-90)=mu_int(mu_int<-90)+180;
mu_int(mu_int>90)=mu_int(mu_int>90)-180;
sigma_int=(max(round(theta_std*0.5),1):0.1:min(round(theta_std*1.5),180));
[mu_best, sigma_best, c, CHI]=vonMisesInterp(bins,h,180,mu_int,sigma_int,false);
Histo.VMmu=mu_best;
Histo.VMsigma=sigma_best;
Histo.VMc=c;
Histo.CHI=CHI;

ShE=-sum(h.*log(h+eps));
ShEspread=exp(ShE);

spread_ent_thres=prctile(h,100-ShEspread/180*100);
h_sp=min(h,spread_ent_thres);
Histo.AngularDistributionSE=h_sp;

spread_mean_thres=mean(h)/2;
Histo.HalfMeanThreshold=spread_mean_thres;
spread_mean=sum(h>=spread_mean_thres);

sigma=2;
g=fspecial('gaussian',[1 2*ceil(3*sigma)+1],sigma);
h2=imfilter(h,g,'circular');
Histo.AngularDistributionSmooth=h2;
[m,ind]=max(h2);
theta_max=bins(ind);
half_max_spread=sum(h2>(m/2));
Histo.Maximum=m;
Histo.MaximumAngle=theta_max;
Histo.HalfMaximumDispersion=half_max_spread;



%---------------------------------------------------------------------
function h=w_hist(bins,x,w)
h=zeros(size(bins));
step=bins(2)-bins(1);
x=round((x-bins(1))/step)+1;
[x,m]=sort(x);
w=w(m);
ind=find(x(2:end)>x(1:end-1))+1;
ind=[1; ind; numel(x)+1];

for i=1:numel(ind)-1
    x2=x(ind(i));
    if and(x2<=numel(h),x2>0)
        var=sum(w(ind(i):(ind(i+1)-1))); 
        h(x2)=var;
    elseif x2>numel(h)
        break
    end
end























