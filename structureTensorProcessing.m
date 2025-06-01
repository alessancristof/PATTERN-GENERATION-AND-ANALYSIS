function [FirstVariation,SecondVariation,Vmax1,Vmax2,StructureSNR,G1,G2,GradientSNR,Graw1,Graw2,G1var,G2var]=structureTensorProcessing(I,NoiseVariance,Spacing,STParametersIn)
% Computes gradient and structure tensor, diagonalizes structure tensor,
% propagates noise variance to obtain unbiased versions of the gradient and
% the structure tensor, computes their SNR.

% INPUTS
% I : n x m x s matrix, stack of s images of n x m size.
% NoiseVariance: n x m x s matrix, noise variance estimate of I.
% Spacing: 1x3 vector, spatial spacing of the voxels of I.
% STParametersIn: structure with custom parameters

% OUTPUTS
% FirstVariation: n x m x s matrix, maximal eigenvalue of the structure tensor.
% SecondVariation: n x m x s matrix, mninimal eigenvalue of the structure tensor.
% Vmax1: n x m x s matrix, first component of the eigenvector corresponding
%   to FirstVariation.
% Vmax2: n x m x s matrix, second component of the eigenvector
%   corresponding to FirstVariation.
% StructureSNR: n x m x s matrix, SNR of the sqare root of the structure
%   tensor trace.
% G1: n x m x s matrix, first component of the gradient after required processing (bias correction, normalization).
% G2: n x m x s matrix, second component of the gradient after required processing (bias correction, normalization).
% GradientSNR: n x m x s matrix, SNR of the gradient norm.
% Graw1: n x m x s matrix, first component of the original gradient.
% Graw2: n x m x s matrix, second component of the original gradient.
% G1var: n x m x s matrix, variance of the first component of the original gradient.
% G2var: n x m x s matrix, variance of the second component of the original gradient.

% defaults
STParameters.GradientDiscretization='Sharr Tap 5';
STParameters.GradientBiasRemoval=true;
STParameters.GradientNormalization=false;
STParameters.StructureTensorRadius=1;
STParameters.StructureTensorSmoothingDirection='2D';
STParameters.StructureBiasRemoval=true;

% overrides the defaults with the custom values
if not(isempty(STParametersIn))
    fl=fields(STParametersIn);
    for i=1:numel(fl)
        STParameters.(fl{i})=STParametersIn.(fl{i});
    end
end

sigma=STParameters.StructureTensorRadius./Spacing;
g1=fspecial('gaussian',[2*ceil(3*sigma(1))+1 1],sigma(1));
g2=fspecial('gaussian',[2*ceil(3*sigma(2))+1 1],sigma(2));
g2=g2';
if strcmp(STParameters.StructureTensorSmoothingDirection,'2D')
    g3=[];
else
    g3=fspecial('gaussian',[2*ceil(3*sigma(3))+1 1],sigma(3));
    g3=permute(g3,[2 3 1]);
end

% Produce the original 2D gradient components and the variance of their noises
[G1,G2,G1var,G2var,gmax]=firGradient2D(I,3,STParameters.GradientDiscretization,NoiseVariance);
Gvar=G1var+G2var; % noise variance of the sum of the components of the gradient
Gsqvar=4*(G1.^2).*G1var+4*(G2.^2).*G2var; % noise variance of the squared norm of the gradient

% Compute variable for bias correction of structure tensor
if STParameters.StructureBiasRemoval
    Mbias1=imfilter(G1var,g1,'replicate');
    Mbias1=imfilter(Mbias1,g2,'replicate');
    Mbias2=imfilter(G2var,g1,'replicate');
    Mbias2=imfilter(Mbias2,g2,'replicate');
    if not(isempty(g3))
        Mbias1=imfilter(Mbias1,g3,'replicate');
        Mbias2=imfilter(Mbias2,g3,'replicate');
    end
    % Mbias=imfilter(Gvar,g1,'replicate');
    % Mbias=imfilter(Mbias,g2,'replicate');
    % if not(isempty(g3))
    %     Mbias=imfilter(Mbias,g3,'replicate');
    % end
end

% Compute the noise variance of the trace of the structure tensor from the
% squared norm of the gradient
TRvar=imfilter(Gsqvar,g1.^2,'replicate');
TRvar=imfilter(TRvar,g2.^2,'replicate');
if not(isempty(g3))
    TRvar=imfilter(TRvar,g3.^2,'replicate');
end
clear Gsqvar

% Computation of the structure tensor from the gradient ------------------------------------------------------------------
M11=imfilter(G1.^2,g1,'replicate');
M11=imfilter(M11,g2,'replicate');
M22=imfilter(G2.^2,g1,'replicate');
M22=imfilter(M22,g2,'replicate');
M12=imfilter(G1.*G2,g1,'replicate');
M12=imfilter(M12,g2,'replicate');
if not(isempty(g3))
    M11=imfilter(M11,g3,'replicate');
    M22=imfilter(M22,g3,'replicate');
    M12=imfilter(M12,g3,'replicate');
end
%-------------------------------------------------------------------------

% structure tensor bias correction
if STParameters.StructureBiasRemoval
    M11=M11-Mbias1;
    M22=M22-Mbias2;
    clear Mbias1 Mbias2
    % M11=M11-Mbias/2;
    % M22=M22-Mbias/2;
    % clear Mbias
end

% diagonalization of unbiased structure tensor
[Lmax,Lmin,Vmax1,Vmax2]=matrixDiagonalize2D(M11,M22,M12);
clear M11 M12 M22

% computation of unbiased eigenvalues
if STParameters.StructureBiasRemoval
    Ltot=max(Lmin+Lmax,0); % nonnegative corrected trace of the structure tensor
    Lmin=max(Lmin,0);
    Lmax=Ltot-Lmin;
else
    Ltot=Lmin+Lmax; % trace of the structure tensor
end

% computation of variations (square root of eigenvalues), consistent with
% the scale of gradient
FirstVariation=sqrt(Lmax);
SecondVariation=sqrt(Lmin);
clear Lmax Lmin

% computation of structure SNR considering as signal the square root of the sum
% of eigenvalues (sqrt of the trace of the structure tensor), to be consistent to
% the scale of the gradient
TRsd=sqrt(TRvar); % noise SD of the structure tensor trace
Ltotrtsd=(sqrt(Ltot+TRsd)-sqrt(max(Ltot-TRsd,0)))/2; % noise SD of the square root of eigenvalues computed by direct propagation 
StructureSNR=sqrt(Ltot)./(Ltotrtsd+eps);
clear Ltot TRvar TRsd Ltotrtsd

% Processing of the original gradient to obtain a biased corrected Gradient
% and the gradient SNR ---------------------------------------------------------------
Gsq=G1.^2+G2.^2; % squared norm of the gradient
Gnorm=sqrt(Gsq); % norm of the gradient

Gnormvar=((G1.^2).*G1var+(G2.^2).*G2var)./(Gsq+eps); % noise variance of the norm of the gradient
Gnormsd=sqrt(Gnormvar); % noise SD of the norm of the gradient

% original unprocessed gradient
Graw1=G1;
Graw2=G2;

if STParameters.GradientBiasRemoval % perform bias correction
    Gnormun=sqrt(max((Gsq-Gvar),0)); % unbiased gradient norm
    GradientSNR=Gnormun./(Gnormsd+eps); % unbiased gradient SNR
    if STParameters.GradientNormalization
        var=1./(Gnorm+eps);
        var(Gnormun==0)=0; 
        G1=G1.*var; % unbiased normalized gradient components
        G2=G2.*var;
    else
        var=Gnormun./(Gnorm+eps);
        G1=G1.*var; % unbiased gradient components
        G2=G2.*var;
    end
else
    GradientSNR=Gnorm./(Gnormsd+eps); % gradient SNR
    if STParameters.GradientNormalization
        var=1./(Gnorm+eps);
        G1=G1.*var; % normalized gradient components
        G2=G2.*var;
    end
end
%---------------------------------------------------------------------
