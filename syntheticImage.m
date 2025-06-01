function [I,Parameters,A]=syntheticImage(Parameters)

% Generate a single sinthetic slice (2D image)

% Parameters (input and output): structure with fields defining the parameters for
% image generation, as output some fields decribe propeties of the generation process 
% I (output): n x m  double matrix, pattern image
% A (output): n x m x 2 double matrix, 2D vectorial field defining the
% direction for fiber drawing in the slice


if strcmp(Parameters.HRGeneration,'on')  % generate the image at twice resolution and downscale it to obtain finer fibers
    Parameters2=Parameters;
    Parameters2.Size=Parameters.Size*2;
    Parameters2.Spacing=Parameters.Spacing/2;
    Parameters2.FiberLength=Parameters.FiberLength*2;
    if strcmp(Parameters.FiberNumber,'auto')
        if not(ischar(Parameters.FiberDensity))
            Parameters2.FiberDensity=Parameters.FiberDensity/2;
        end
    else
        Parameters2.FiberNumber=Parameters.FiberNumber/2;
    end

    [I,Parameters,A]=syntheticImage_pre(Parameters2);
    I=double(I);
    I=(I(1:2:end,:)+I(2:2:end,:))/2;
    I=(I(:,1:2:end)+I(:,2:2:end))/2;
    A=(A(1:2:end,:,:)+A(2:2:end,:,:))/2;
    A=(A(:,1:2:end,:)+A(:,2:2:end,:))/2;

    Parameters.Size=Parameters.Size/2;
    Parameters.Spacing=Parameters.Spacing*2;
    Parameters.FiberLength=Parameters.FiberLength/2;
    if strcmp(Parameters.FiberNumber,'auto')
        if not(ischar(Parameters.FiberDensity))
            Parameters.FiberDensity=Parameters.FiberDensity*2;
        end
    else
        Parameters.FiberNumber=Parameters.FiberNumber*2;
    end

else
    [I,Parameters,A]=syntheticImage_pre(Parameters);
end



%-------------------------------------------------------------------------------------------------------------------------
function [I,Parameters,A]=syntheticImage_pre(Parameters)


siz=Parameters.Size;

% generate modulation image
M=randn(siz,'double');
BackgroundSigma=Parameters.BackgroundSigma./Parameters.Spacing(1:2);
M=imgaussfilt(M,BackgroundSigma,'FilterDomain','spatial','Padding','circular','FilterSize',2*ceil(3*BackgroundSigma)+1);
M=M-prctile(M,Parameters.BackgroundPercentage,"all");
M=M/prctile(M,100-Parameters.TissuePercentage,"all");
M=min(max(M,0),1);
M=imgaussfilt(M,5,'FilterSize',2*ceil(4*5)+1,'Padding','circular');
% figure
% imshow(M,[])

if strcmp(Parameters.AngleMean,'random')
    AngleMean=(rand-0.5)*pi;
else
    AngleMean=Parameters.AngleMean/180*pi;
end
AngleSD=Parameters.AngleSD/180*pi;
AngleSigma=Parameters.AngleSigma./Parameters.Spacing(1:2);
AngleHighSigma=Parameters.AngleSigma./Parameters.Spacing(1:2);
if strcmp(Parameters.Wavelength,'auto')
    Wavelength=[];
else
    Wavelength=Parameters.Wavelength./Parameters.Spacing(1);
end
% generate the directional field trying to attain the target angular SD (AngleSD).
% Accuracy of the effective angular SD (RealAngleSD) is improved by a repeated attempts approach.
ang_deg_tol=0.1; % tolerance for field angular dispersion for multiple attempt at field generation 
var=Parameters.AngleFieldOctaves(end)/sum(Parameters.AngleFieldOctaves,"all");
var2=2^(numel(Parameters.AngleFieldOctaves)-1);
ang_tol=inf;
for i=1:20
    [A2,RealAngleMean2,RealAngleSD2,ActualWavelength2]=angleGeneration(siz,Parameters.SinWave,Wavelength,AngleMean,AngleSD*var,AngleSigma*var2,AngleHighSigma*var2,Parameters.AngleRangeLimiter,Parameters.FieldCoalescence,1,3);
    ang_tol2=abs(RealAngleSD2-AngleSD*var);
    if ang_tol2<ang_tol
        A=A2;
        RealAngleMean=RealAngleMean2;
        RealAngleSD=RealAngleSD2;
        ActualWavelength=ActualWavelength2;
        ang_tol=ang_tol2;
    end
    if ang_tol/pi*180<ang_deg_tol
        if i>1
            [num2str(i) ' attempts at field generation']
        end
        break
    end
    if i==20
        '20 attempts at field generation exausted at RealAngleSD error (deg):'
        (RealAngleSD-AngleSD*var)/pi*180
    end
end
% multiscale routine for field generation
for i=(numel(Parameters.AngleFieldOctaves)-1):-1:1
    var=sum(Parameters.AngleFieldOctaves(i:end),"all")/sum(Parameters.AngleFieldOctaves,"all");
    % var=Parameters.AngleFieldOctaves(i)/sum(Parameters.AngleFieldOctaves,"all");
    var2=2^(i-1);
    ang_tol=inf;
    for j=1:20
        [A2,RealAngleMean2,RealAngleSD2,ActualWavelength2]=angleGeneration(siz,Parameters.SinWave,Wavelength,A,AngleSD*var,AngleSigma*var2,AngleHighSigma*var2,Parameters.AngleRangeLimiter,Parameters.FieldCoalescence,1,1.5);
        ang_tol2=abs(RealAngleSD2-AngleSD*var);
        if ang_tol2<ang_tol
            Ab=A2;
            RealAngleMean=RealAngleMean2;
            RealAngleSD=RealAngleSD2;
            ActualWavelength=ActualWavelength2;
            ang_tol=ang_tol2;
        end
        if ang_tol/pi*180<ang_deg_tol
            if i>1
                [num2str(i) ' attempts at field generation']
            end
            break
        end
        if i==20
            '20 attempts at field generation exausted at RealAngleSD error (deg):'
            (RealAngleSD-AngleSD*var)/pi*180
        end
    end
    A=Ab;
end

% for the sinusoidal approach
ActualWavelength=ActualWavelength*mean(Parameters.Spacing(1:2));

ang=atan2(-A(:,:,1),A(:,:,2)); % angle in image reference
[PostModulationAngleMean,PostModulationAngleSD]=circular_stat_ang(ang(:),pi,M(:));
RealAngleMean=RealAngleMean/pi*180;
RealAngleSD=RealAngleSD/pi*180;
PostModulationAngleMean=PostModulationAngleMean/pi*180;
PostModulationAngleSD=PostModulationAngleSD/pi*180;
TheoreticalAngle=AngleMean/pi*180;


% draw the image using the directional field as guidance
if numel(Parameters.FiberLength)>1
    Parameters2=Parameters;
    Parameters2.FiberLength=Parameters.FiberLength(1);
    [I,Parameters2]=syntheticImage_pre2(A,Parameters2);
    I=I*Parameters.FiberLengthIntensity(1);
    Parameters.ActualFiberNumber=Parameters2.ActualFiberNumber;
    for i=2:numel(Parameters.FiberLength)
        Parameters2.FiberLength=Parameters.FiberLength(i);
        [I2,Parameters2]=syntheticImage_pre2(A,Parameters2);
        I=I+I2*Parameters.FiberLengthIntensity(i);
        Parameters.ActualFiberNumber=cat(2,Parameters.ActualFiberNumber,Parameters2.ActualFiberNumber);
    end
else
    [I,Parameters]=syntheticImage_pre2(A,Parameters);
end

% adjust contrast and average intensity to match desired input AverageSignal
p=prctile(I,[0.1 99.9],"all");
if p(2)>p(1)
    I=max(I-p(1),0);
end
if strcmp(Parameters.AverageSignalDefinition,'mean')
    I=I/mean(I,"all")*Parameters.AverageSignal;
else
    I=I/prctile(I,Parameters.AverageSignalDefinition,"all")*Parameters.AverageSignal;
end

% modulate the image
I=I.*M; 

% add baseline and noise
I=I+Parameters.AverageSignal/Parameters.SDR+randn(size(I))*(Parameters.AverageSignal/Parameters.SNRWhite) + randn(size(I)).*sqrt(I)*(sqrt(Parameters.AverageSignal)/Parameters.SNRShot) + randn(size(I)).*I/Parameters.SNRMult;

Parameters.ActualWavelength=ActualWavelength;
Parameters.DarkCurrent=Parameters.AverageSignal/Parameters.SDR;
Parameters.WhiteNoiseCoefficient=(Parameters.AverageSignal/Parameters.SNRWhite)^2;
Parameters.ShotNoiseCoefficient=Parameters.AverageSignal/(Parameters.SNRShot^2);
Parameters.MultNoiseCoefficient=1/(Parameters.SNRMult^2);
Parameters.RealAngleMean=RealAngleMean;
Parameters.RealAngleSD=RealAngleSD;
Parameters.PostModulationAngleMean=PostModulationAngleMean;
Parameters.PostModulationAngleSD=PostModulationAngleSD;
Parameters.TheoreticalAngle=TheoreticalAngle;

%------------------------------------------------------------------
function [I,Parameters]=syntheticImage_pre2(A,Parameters)
% draw the image I using the directional field A as guidance
% multiple scales of fiber thickness are controlled by OctaveIntensity

siz=Parameters.Size;
MeanSpacing=mean(Parameters.Spacing(1:2));
I=zeros(siz);
for i=1:numel(Parameters.OctaveIntensity)
    if Parameters.OctaveIntensity(i)>0
        Parameters2=Parameters;
        if i>1
            A2=imresize(A,2^(-i+1),"bilinear");
            if strcmp(Parameters.FiberLengthRescale,'on')
                Parameters2.FiberLength=Parameters.FiberLength*2^(-i+1);
            else
                Parameters2.FiberLength=Parameters.FiberLength;
            end
            if strcmp(Parameters.FiberNumber,'auto')
                Parameters2.FiberNumber=Parameters.FiberNumber;
            else
                Parameters2.FiberNumber=round((Parameters.FiberLength/Parameters2.FiberLength)*Parameters.FiberNumber*2^(-2*i+2));
            end
            Parameters2.SpatialJitter=(Parameters.AngleSpatialJitter/MeanSpacing)*2^(-i+1);
            
        else
            A2=A;
        end
        [I2,ActualFiberNumber2]=syntheticImage_aux(A2,Parameters2);
        ActualFiberNumber(i)=ActualFiberNumber2;
        if i>1
            I2=imresize(I2,size(I),"bicubic");
        end
        sig=0.5*2^i;
        I2=imgaussfilt(I2,sig,'FilterDomain','spatial','Padding','circular','FilterSize',2*ceil(4*sig)+1);

        p=prctile(I2,[0.1 99.9],"all");
        if p(2)>p(1)
            I2=max((I2-p(1))/(p(2)-p(1)),0);
        end
        I=I+Parameters.OctaveIntensity(i)*(I2/mean(I2,"all"));
    end
end
Parameters.ActualFiberNumber=ActualFiberNumber;

%-------------------------------------------------------------------------------------------------------------------------
function [I,FiberNumber]=syntheticImage_aux(A,Parameters)
% draw the image I using the directional field A as guidance at the minimal
% scale of fiber thickness

interp_method='linear';
n=1;
sp=Parameters.SeedSpacing;
siz=size(A);
siz=siz(1:2);
step=0.5;
num=round(Parameters.FiberLength/step/2);


if strcmp(Parameters.FiberModulationWL,'auto')
    frmax=max(floor(Parameters.FiberLength/8),1);
    Parameters.FiberModulationWL=Parameters.FiberLength./[1:frmax];
end
fiber_mod_num=Parameters.FiberModulationWL;
switch Parameters.FiberModulationSpectrum
    case 'white'
        fiber_mod_octaves=ones(size(fiber_mod_num));
    case 'pink'
        fiber_mod_octaves=sqrt(fiber_mod_num);
    case 'red'    % Brownian
        fiber_mod_octaves=fiber_mod_num;
end
cutoffmask=Parameters.FiberModulationWL>Parameters.FiberModulationCutoff;
if min(cutoffmask)==0
    fiber_mod_octaves(cutoffmask)=max(fiber_mod_octaves(not(cutoffmask)));
else
    fiber_mod_octaves=ones(size(fiber_mod_num));
end

num0=num;

if strcmp(Parameters.FiberLengthDistribution,'Poisson')
    pd=makedist('Poisson','lambda',num0);
end
% pd=makedist('Exponential','mu',num0);
% pd=makedist('Rayleigh','B',num0*sqrt(2/pi));


I=zeros(siz);
I0=zeros(siz);

switch Parameters.SeedMode
    case 'spaced'
        [Rs,Cs]=ndgrid([1:sp:siz(1)+sp/2],[1:sp:siz(2)+sp/2]);
        Rs=Rs(:);
        Cs=Cs(:);
        Rs=Rs+(rand(size(Rs))-0.5)*sp*Parameters.SeedRandomization;
        Cs=Cs+(rand(size(Rs))-0.5)*sp*Parameters.SeedRandomization;
        mask=or(Rs<2,Rs>(siz(1)-1));
        Rs(mask)=[];
        Cs(mask)=[];
        mask=or(Cs<2,Cs>(siz(2)-1));
        Rs(mask)=[];
        Cs(mask)=[];
        Parameters.FiberNumber=numel(Rs);
end

flag_zeropixel=true;

if strcmp(Parameters.FiberNumber,'auto')
    if strcmp(Parameters.FiberDensity,'packed')
        FiberNumber=round(siz(1)*siz(2)/Parameters.FiberLength*2);
    elseif strcmp(Parameters.FiberDensity,'loose')
        FiberNumber=round(siz(1)*siz(2)/Parameters.FiberLength*2);
    else
        FiberNumber=round(siz(1)*siz(2)/Parameters.FiberLength*Parameters.FiberDensity);
    end
else
    FiberNumber=Parameters.FiberNumber;
    Parameters.FiberDensity='not used';
end

SE=strel("rectangle",[3 3]);

% speedup_factor=10;
speedup_factor=max(round(FiberNumber/Parameters.BlendingSpeedupIterMax),1);

% preliminary  creation of stencils and corrsponding positions along a single fiber
current_fiber_number=0;
for j=1:round(FiberNumber/speedup_factor)
    Vq=[];
    r2=[];
    c2=[];
    for spk=1:speedup_factor
        current_fiber_number=current_fiber_number+1;
        ajit=randn*Parameters.AngleJitter/180*pi;
        sjit=randn(1,2)*Parameters.AngleSpatialJitter;
        switch Parameters.SeedMode
            case 'spaced'
                r0=Rs(current_fiber_number);
                c0=Cs(current_fiber_number);
            case 'random'
                if and(Parameters.SeedAvoidance,flag_zeropixel)
                    for i=1:1000
                        r0=1+rand*(siz(1)-1);
                        c0=1+rand*(siz(2)-1);
                        dr=(round(r0)-1):(round(r0)+1);
                        dc=(round(c0)-1):(round(c0)+1);
                        dr=1+mod(dr-1,siz(1));
                        dc=1+mod(dc-1,siz(2));
                        if max(I(dr,dc),[],"all")==0
                            break
                        end
                        if i==1000
                            if strcmp(Parameters.FiberDensity,'loose')
                                [row,col]=find(imerode(I==0,SE));
                                if isempty(row)
                                    flag_zeropixel=false;
                                else
                                    ind=round(1+rand*(numel(row)-1));
                                    r0=row(ind)+(rand-0.5);
                                    c0=col(ind)+(rand-0.5);
                                end
                            else
                                for k=1:1000
                                    r0=1+rand*(siz(1)-1);
                                    c0=1+rand*(siz(2)-1);
                                    if I(round(r0),round(c0))==0
                                        break
                                    end
                                    if k==1000
                                        [row,col]=find(I==0);
                                        if isempty(row)
                                            flag_zeropixel=false;
                                        else
                                            ind=round(1+rand*(numel(row)-1));
                                            r0=row(ind)+(rand-0.5);
                                            c0=col(ind)+(rand-0.5);
                                        end
                                    end
                                end
                            end
                        end
                    end
                else
                    if or(strcmp(Parameters.FiberDensity,'packed'),strcmp(Parameters.FiberDensity,'loose'))
                        FiberNumber=current_fiber_number;
                        break
                    else
                        [~,ind]=min(I,[],"all");
                        [r0,c0] = ind2sub(size(I),ind);
                        r0=r0+(rand-0.5);
                        c0=c0+(rand-0.5);
                    end
                end
        end
        
        r=r0;
        c=c0;
        v=1;

        
        if strcmp(Parameters.FiberLengthDistribution,'Poisson')
            num=max(round(random(pd)),5);
        elseif strcmp(Parameters.FiberLengthDistribution,'Chisquare6')  %  x^2*exp(-x)
            num=max(round(random('Chisquare',6)/6*num0),5);
        elseif strcmp(Parameters.FiberLengthDistribution,'Chisquare8') %  x^3*exp(-x)
            num=max(round(random('Chisquare',8)/8*num0),5);
        elseif strcmp(Parameters.FiberLengthDistribution,'Maxwell')  % == Chi3  x^2*exp(-x^2)
            num=max(round(sqrt(random('Chisquare',3)/3)*num0),5);
        elseif strcmp(Parameters.FiberLengthDistribution,'Chi4')   % x^3*exp(-x^2)
            num=max(round(sqrt(random('Chisquare',4)/4)*num0),5);
        elseif strcmp(Parameters.FiberLengthDistribution,'Chi5')     % x^4*exp(-x^2)
            num=max(round(sqrt(random('Chisquare',5)/5)*num0),5);
        end

        if Parameters.FiberModulationIntensity>0
            
            x1=[1:num]*step;
            x2=[0:-1:-num+1]*step;
            phas=rand(size(fiber_mod_octaves))*2*pi;
            fiber_mod1=fiber_mod_octaves(1)*cos(x1/fiber_mod_num(1)*2*pi+phas(1));
            fiber_mod2=fiber_mod_octaves(1)*cos(x2/fiber_mod_num(1)*2*pi+phas(1));
            for i=2:numel(fiber_mod_octaves)
                fiber_mod1=fiber_mod1+fiber_mod_octaves(i)*cos(x1/fiber_mod_num(i)*2*pi+phas(i));
                fiber_mod2=fiber_mod2+fiber_mod_octaves(i)*cos(x2/fiber_mod_num(i)*2*pi+phas(i));
            end
            fmax=max(max(fiber_mod1(:)),max(fiber_mod2(:)));
            fmin=min(min(fiber_mod1(:)),min(fiber_mod2(:)));
            fiber_mod1=(fiber_mod1-fmin)/(fmax-fmin);
            fiber_mod2=(fiber_mod2-fmin)/(fmax-fmin);
            fiber_mod1=1-Parameters.FiberModulationIntensity+Parameters.FiberModulationIntensity*fiber_mod1;
            fiber_mod2=1-Parameters.FiberModulationIntensity+Parameters.FiberModulationIntensity*fiber_mod2;
            
        end

        for i=1:num
            if Parameters.FiberModulationIntensity>0
                switch Parameters.FiberFading
                    case 'circle'
                        % v(i,1)=sqrt(1-(i/num)^2);
                        v(i,1)=sqrt(1-(i/num)^2)*fiber_mod1(i);
                    case 'sin'
                        v(i,1)=sin(pi/2+i/num*pi/2)*fiber_mod1(i);
                    otherwise
                        v(i,1)=1;
                end
            else
                switch Parameters.FiberFading
                    case 'circle'
                        v(i,1)=sqrt(1-(i/num)^2);
                    case 'sin'
                        v(i,1)=sin(pi/2+i/num*pi/2);
                    otherwise
                        v(i,1)=1;
                end
            end

            rf=floor(r(i)+sjit(1));
            cf=floor(c(i)+sjit(2));
            rd=r(i)+sjit(1)-rf;
            cd=c(i)+sjit(2)-cf;
            rf1=rf+1;
            cf1=cf+1;
            rf=1+mod(rf-1,siz(1));
            cf=1+mod(cf-1,siz(2));
            rf1=1+mod(rf1-1,siz(1));
            cf1=1+mod(cf1-1,siz(2));
            ao=(1-rd)*(1-cd)*A(rf,cf,:) + rd*(1-cd)*A(rf1,cf,:) + (1-rd)*cd*A(rf,cf1,:) + rd*cd*A(rf1,cf1,:);

            if size(ao,3)==2
                a=atan2(ao(1),ao(2))+ajit;
            else
                a=ao+ajit;
            end

            if i<num
                do=[sin(a) cos(a)];
                if i>1
                    if do(1)*d(1)+do(2)*d(2)>0
                        d=do;
                    else
                        d=-do;
                    end
                else
                    d=do;
                end
                r(i+1,1)=r(i)+d(1)*step;
                c(i+1,1)=c(i)+d(2)*step;
                
            end
        end

        [Vqb,r2b,c2b]=pencil_aux(r,c,v);
        Vq=cat(3,Vq,Vqb);
        r2=cat(1,r2,r2b);
        c2=cat(1,c2,c2b);

        r=r0;
        c=c0;
        v=1;
        for i=1:num
            if Parameters.FiberModulationIntensity>0
                switch Parameters.FiberFading
                    case 'circle'
                        v(i,1)=sqrt(1-(i/num)^2)*fiber_mod2(i);
                    case 'sin'
                        v(i,1)=sin(pi/2+i/num*pi/2)*fiber_mod2(i);
                    otherwise
                        v(i,1)=1;
                end
            else
                switch Parameters.FiberFading
                    case 'circle'
                        v(i,1)=sqrt(1-(i/num)^2);
                    case 'sin'
                        v(i,1)=sin(pi/2+i/num*pi/2);
                    otherwise
                        v(i,1)=1;
                end
            end

            
            rf=floor(r(i)+sjit(1));
            cf=floor(c(i)+sjit(2));
            rd=r(i)+sjit(1)-rf;
            cd=c(i)+sjit(2)-cf;
            rf1=rf+1;
            cf1=cf+1;
            rf=1+mod(rf-1,siz(1));
            cf=1+mod(cf-1,siz(2));
            rf1=1+mod(rf1-1,siz(1));
            cf1=1+mod(cf1-1,siz(2));
            ao=(1-rd)*(1-cd)*A(rf,cf,:) + rd*(1-cd)*A(rf1,cf,:) + (1-rd)*cd*A(rf,cf1,:) + rd*cd*A(rf1,cf1,:);
            if size(ao,3)==2
                a=atan2(ao(1),ao(2))+ajit;
            else
                a=ao+ajit;
            end
            if i<num
                do=[sin(a) cos(a)];
                if i>1
                    if do(1)*d(1)+do(2)*d(2)>0
                        d=do;
                    else
                        d=-do;
                    end
                else
                    d=do;
                end
                r(i+1,1)=r(i)-d(1)*step;
                c(i+1,1)=c(i)-d(2)*step;
                
            end
        end
        [Vqb,r2b,c2b]=pencil_aux(r,c,v);
        Vq=cat(3,Vq,Vqb);
        r2=cat(1,r2,r2b);
        c2=cat(1,c2,c2b);
    end

    % finalize an image containing just the single fiber
    I1=I0;
    I1=pencil(I1,r2,c2,Vq,n);

    % combine the single fiber image with the overall image I
    switch Parameters.FiberBlendingMode
        case 'max'
            I=max(I,I1);
        case 'add'
            if Parameters.FiberBlendingExponent==1
                I=I+I1;
            else
                I=I+I1.^Parameters.FiberBlendingExponent;
            end
        case 'exp'
            I=I+I1;
        case 'addsub'
            if rand>0.5
                I=I+I1;
            else
                I=I-I1;
            end
    end
end

switch Parameters.FiberBlendingMode
    case 'add'
        if not(Parameters.FiberBlendingExponent==1)
            I=I.^(1/Parameters.FiberBlendingExponent);
        end
    case 'exp'
        FiberDensity=FiberNumber*Parameters.FiberLength/(Parameters.Size(1)*Parameters.Size(2));
        I=1-exp(I/(mean(I,"all")/FiberDensity)*log(0.5));
end


%-------------------------------------------------------------------------------------------------------------------------
function I=pencil(I,r2,c2,Vq,n)
% place the strokes from the stencils in the image to form a single fiber

siz=size(I);

for i=1:numel(r2)
    indr=r2(i)-n:r2(i)+n;
    indc=c2(i)-n:c2(i)+n;
    indr=1+mod(indr-1,siz(1));
    indc=1+mod(indc-1,siz(2));
    I(indr,indc)=max(I(indr,indc),Vq(:,:,i));
end


%--------------------------------------------------------------------------------------------------
function [Vq,r2,c2]=pencil_aux(r,c,v)
% draw a single stroke in a stencil


r2=round(r);
c2=round(c);

dr=r-r2;
dc=c-c2;


h1=0.625;

h2=(1-h1)/2;

a=2*permute(dr,[3 2 1]);
a1=max(-a,0);
a2=abs(a);
a3=max(a,0);
s1=h2*(1-a2)+0.5*a1;
s2=h1*(1-a2)+0.5*a2;
s3=h2*(1-a2)+0.5*a3;
V0=cat(1,s1,s2,s3);

a=2*permute(dc,[3 2 1]);
a1=max(-a,0);
a2=abs(a);
a3=max(a,0);
s1=h2*(1-a2)+0.5*a1;
s2=h1*(1-a2)+0.5*a2;
s3=h2*(1-a2)+0.5*a3;
V=cat(2,s1.*V0,s2.*V0,s3.*V0);

v=permute(v,[3 2 1]);
Vq=V.*v;


%-------------------------------------------------------------------------------------------------------------------------
function [A,ang_mean,ang_sd,WL]=angleGeneration(siz,SinWave,WL,AngleMean,AngleSD,LowSigma,HighSigma,AngleRangeLimiter,FieldCoalescence,AngleFieldOctaves,alphaLimiter)
% generate the directional field A


if SinWave>0
    Xim=cos(AngleMean);
    Yim=sin(AngleMean);
    Xint=-Yim;
    Yint=Xim;
    if isempty(WL)
        WL=mean(LowSigma)*4;
        WL=siz(1)/max(round(siz(1)/WL),1);
    end
    
    [Xg,Yg]=ndgrid(1:siz(1),1:siz(2));
    
    K=[Xint Yint]/WL;
    V=sin((Xg*K(1)+Yg*K(2))*2*pi+rand*2*pi);
    

    sigma=WL/4*1;
    w1=tukeywin(siz(1),4*sigma/siz(1));
    w2=tukeywin(siz(2),4*sigma/siz(2));
    W=w1*w2';

    V2=imgaussfilt(V,sigma,'FilterSize',2*ceil(4*sigma)+1,'Padding','circular');
    
    V2=V2/std(V2,0,"all")*std(V,0,"all");
    V=W.*V+(1-W).*V2;
   
    fact=1;
    if AngleSD>eps
        for i=1:10
            Gx2=Xint-fact*Yint*V;
            Gy2=Yint+fact*Xint*V;
            ang=atan2(-Gx2,Gy2); % angle is in image reference system
            [~,ang_sd]=circular_stat_ang(ang(:),pi);
            fact=fact*(AngleSD/(ang_sd+eps));
        end
        
    end
    Gxwave=Xint-fact*Yint*V;
    Gywave=Yint+fact*Xint*V;
end



if SinWave<1

    Flag_Equalize=false;

    
    N=2*rand(siz)-1;
    A=imgaussfilt(N,LowSigma,'FilterSize',2*ceil(4*LowSigma)+1,'Padding','circular');
    if HighSigma<inf
        Ab=imgaussfilt(A,HighSigma,'FilterSize',2*ceil(4*HighSigma)+1,'Padding','circular');
        A=A-Ab;
    end
   
    for i=1:(AngleFieldOctaves-1)
        
        B=imgaussfilt(N,(2^i)*LowSigma,'FilterSize',2*ceil(4*(2^i)*LowSigma)+1,'Padding','circular');
        A=A/std(A,0,"all")/2+B/std(B,0,"all");
    end

    [Gx,Gy] = imgradientxy(A);
    G=sqrt(Gx.^2+Gy.^2);
   
    if Flag_Equalize
        G2=imgaussfilt(G,LowSigma*2,'FilterSize',2*ceil(4*LowSigma*2)+1,'Padding','circular');
        Gx=Gx./(G2+eps);
        Gy=Gy./(G2+eps);
        G=sqrt(Gx.^2+Gy.^2);
    end

    Gavg=mean(G,"all");

    % field of input angle (in image reference system)
    if size(AngleMean,3)>1
        Xim=AngleMean(:,:,2);
        Yim=-AngleMean(:,:,1);
        var=sqrt(Xim.^2+Yim.^2)+eps;
        Xim=Xim./var;
        Yim=Yim./var;
    else
        Xim=cos(AngleMean);
        Yim=sin(AngleMean);
    end

    % fine tuning of angular spread by iterated adjustments
    ang=atan2(Yim,Xim);
    [~,ang_sd0]=circular_stat_ang(ang(:),pi);
    
    
    if AngleSD>eps
        AngleSD2=max(AngleSD-ang_sd0,0.01);
        for i=1:10
            switch AngleRangeLimiter
                case 'exp'
                    var=(1-exp(AngleSD2/0.375*G/Gavg*log(0.5)))./(G+eps)*1.0;
                case 'none'
                    var=G/Gavg*0.5./(G+eps)*AngleSD2/0.4265;
            end
            Gx2=Gx.*var;
            Gy2=Gy.*var;
            % input angle is in image reference system
            Gx2=-Yim+Gx2;
            Gy2=Xim-Gy2;
            ang=atan2(-Gx2,Gy2); % angle is in image reference system
            [~,ang_sd]=circular_stat_ang(ang(:),pi);
            fact=max(AngleSD-ang_sd0,0)/max(ang_sd-ang_sd0,eps);
            if i==1
                Gx3=-Yim-Gx.*var;
                Gy3=Xim+Gy.*var;
                ang=atan2(-Gx3,Gy3); % angle is in image reference system
                [~,ang_sd]=circular_stat_ang(ang(:),pi);
                fact2=max(AngleSD-ang_sd0,0)/max(ang_sd-ang_sd0,eps);
                if abs(1-fact2)<abs(1-fact)
                    Gx=-Gx;
                    Gy=-Gy;
                    fact=fact2;
                end
            end
            fact=max(min(fact,2),0.5);
            AngleSD2=min(AngleSD2*fact,alphaLimiter);
            
        end
        
    else
        AngleSD2=0;
    end

    switch AngleRangeLimiter
        case 'exp'
            var=(1-exp(AngleSD2/0.375*G/Gavg*log(0.5)))./(G+eps)*1.0;
        case 'none'
            var=G/Gavg*0.5./(G+eps)*AngleSD2/0.4265;
    end
    Gx=Gx.*var;
    Gy=Gy.*var;
    G=sqrt(Gx.^2+Gy.^2);
    A=atan2(Gy,Gx);
    A=A*1+pi/2*FieldCoalescence;
    Gx=cos(A).*G;
    Gy=sin(A).*G;

    % angle is in image reference system
    Gxunc=-Yim+Gx;
    Gyunc=Xim-Gy;
end

if SinWave==0
    Gx=Gxunc;
    Gy=Gyunc;
elseif SinWave==1
    Gx=Gxwave;
    Gy=Gywave;
else
    a=1-SinWave;
    b=SinWave;
    Gx=(a*Gxunc+b*Gxwave)/sqrt(a^2+b^2);
    Gy=(a*Gyunc+b*Gywave)/sqrt(a^2+b^2);

end

ang=atan2(-Gx,Gy); % angle is in image reference system
[ang_mean,ang_sd]=circular_stat_ang(ang(:),pi);


A=cat(3,Gx,Gy);










