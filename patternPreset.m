function Parameters=patternPreset(pattern_preset,CustomParameters)

% Create a structure with fields defining the parameters for image
% generation according to presets and custom user's inputs.

% pattern_preset (input): char defining a preset for parameters ('I Straight Lines','II Coherent Waves'
% ,'III Uncoherent Waves', 'IV Coherent Wavelets', 'V Uncoherent Wavelets',
% 'VI Entangled','VII Chaotic', 'VIII Composite','IX Multiscale Waves',
% 'X Fuzzy Multiscale Waves','Braids','Whirlpools','Nodes','Curls'). 
% CustomParameters (optional input): structure containing fields for custom parameter

% values to assign to the preset Parameter variable 
% Parameters (output): structure with fields containing the final generative parameters

% Explanation of some relevant field:

% Size: size in pixels of each slice.
% SliceNumber: number of slices.
% Spacing: 1x3 vector with pixel spacing along x and y and slice spacing
% along z.

% Vector field generation
% SinWave: number ranging [0 - 1], use a combination of the blurred noise (0) and the sinusoidal wave (1) to
% perturbe the directional field.
% HRGeneration: 'on' or 'off', use the high resolution generation
% AngleFieldOctaves: 1xn vector, weights for multiscale field generation approach, from the smallest scale to the largest. 
% AngleSigma: sigma of Gaussian low pass for directional variation in vector field (scale of waves).
% AngleHighSigma: sigma of Gaussian high pass for directional variation in vector field (scale of waves).
% AngleMean: target mean angle (deg) of main fiber firection, if char 'progressive' the slices span the 180 deg half circle;
% AngleSD: targed angular spread (SD, deg) for directional field. 
% AngleRangeLimiter: 'none' or 'exp', if when perturbing the directional
% vectorial field a limit is enforced to avoid "looping" paths.
% FieldCoalescence: [0 - 1], if > 0 rotate the angle of the directional
% field up to 90 deg (1), so that fibers form spiral and converging patterns.

% Drawing of the fibers
% OctaveIntensity: 1xn vector, intensity weights for generating fibers of
% n different thickness (by scales of 2), from thinnest to thickest
% FiberLength: scalar or 1xn vector, length of fibers in px, if a vector has same size of FiberLengthIntensity.
% FiberLengthIntensity: scalar or 1xn vector, weights for multiple length of fibers, if a vector has same size of FiberLength.
% AngleSpatialJitter: sigma of Gaussian distribution for randomly spatially shifting each fiber (use for chaotic patterns).

% Intensity of image
% AverageSignalDefinition: how the "average" value for normalizing the output image
% is computed (char 'mean' -> mean, number [0 - 100] -> percentile).
% AverageSignal: the "average" signal imposed to the output image;


% Modulation of intensity
% TissuePercentage: percentage of full intensity tissue in the image after
% modulation
% BackgroundPercentage: percentage of zero intensity tissue in the image after
% modulation;
% BackgroundSigma: sigma for the Gaussian blur that define the scale of intensity modulation 
% modulation

% Addition of noise
% SNRShot: SNR for "shot" noise (inf for no noise)
% SNRMult: SNR for "multiplicative" noise (inf for no noise)
% SNRWhite: SNR for base white noise (inf for no noise)
% SDR: "average" signal to "Dark Current" Ratio.




if nargin<2
    CustomParameters=[];
end

% define defaults
Parameters.Size=[512 512];
Parameters.SliceNumber=1;
Parameters.Spacing=[0.34 0.34 1];
Parameters.Class='uint16';
Parameters.HRGeneration='off';
Parameters.TissuePercentage=50;
Parameters.BackgroundPercentage=10;
Parameters.BackgroundSigma=20;
Parameters.AverageSignalDefinition='mean';
Parameters.AverageSignal=1000;
Parameters.SNRShot=10;
Parameters.SNRMult=10;
Parameters.SNRWhite=10;
Parameters.SDR=5; % Signal to Dark Current Ratio
Parameters.NoiseEstimateBypass='off';

Parameters.SinWave=0;
Parameters.Wavelength=29.01;
Parameters.AngleSigma=8;
Parameters.AngleHighSigma=8;
Parameters.AngleMean='progressive';
Parameters.AngleSD=20;
Parameters.AngleJitter=0;
Parameters.AngleSpatialJitter=0;
Parameters.AngleRangeLimiter='exp';
Parameters.FieldCoalescence=0;
Parameters.AngleFieldOctaves=1;

Parameters.OctaveIntensity=sqrt([1 1 1]);
Parameters.SeedMode='random';
Parameters.SeedAvoidance='true';
Parameters.SeedSpacing=10;
Parameters.SeedRandomization=0.5;
Parameters.FiberNumber='auto';
Parameters.FiberDensity='loose';
Parameters.FiberLength=128;
Parameters.FiberFading='sin';
Parameters.FiberBlendingMode='add';
Parameters.FiberBlendingExponent=2;
Parameters.FiberLengthDistribution='none';
Parameters.FiberLengthIntensity=1;
Parameters.FiberLengthRescale='on';

Parameters.FiberModulationIntensity=0;
Parameters.FiberModulationWL='auto';
Parameters.FiberModulationSpectrum='pink';
Parameters.FiberModulationCutoff=64;
Parameters.BlendingSpeedupIterMax=1024;



% define presets
switch pattern_preset
    case 'I Straight Lines'
        Parameters.OctaveIntensity=sqrt([1 1 1]);
        Parameters.SinWave=0;
        Parameters.FiberDensity='loose';
        Parameters.FiberLength=256;
        Parameters.AngleSigma=20;
        Parameters.AngleHighSigma=inf;
        Parameters.AngleMean='progressive';
        Parameters.AngleSD=0;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=0;
        Parameters.AngleRangeLimiter='exp';
        Parameters.FieldCoalescence=0;
        Parameters.FiberBlendingMode='add';
        Parameters.FiberBlendingExponent=2;
    case 'II Coherent Waves'
        Parameters.OctaveIntensity=sqrt([1 1 1]);
        Parameters.SinWave=1;
        Parameters.Wavelength=29.01;
        % Parameters.FiberDensity=0.7;
        Parameters.FiberDensity='loose';
        Parameters.FiberLength=256;
        % Parameters.AngleSigma=7.2533; % to obtain wavelength (4*sigma) an exact fraction of image size
        % Parameters.AngleHighSigma=7.2533;
        Parameters.AngleMean='progressive';
        Parameters.AngleSD=20;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=0;
        Parameters.AngleRangeLimiter='exp';
        Parameters.FieldCoalescence=0;
        Parameters.FiberBlendingMode='add';
        Parameters.FiberBlendingExponent=2;
    case 'III Uncoherent Waves'
        Parameters.OctaveIntensity=sqrt([1 1 1]);
        Parameters.SinWave=0;
        % Parameters.FiberDensity=0.7;
        Parameters.FiberDensity='loose';
        Parameters.FiberLength=128;
        Parameters.AngleSigma=8;
        Parameters.AngleHighSigma=8;
        Parameters.AngleMean='progressive';
        Parameters.AngleSD=20;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=0;
        Parameters.AngleRangeLimiter='exp';
        Parameters.FieldCoalescence=0;
        Parameters.FiberBlendingMode='add';
        Parameters.FiberBlendingExponent=2;
    case 'IV Coherent Wavelets'
        Parameters.OctaveIntensity=sqrt([1 1]);
        Parameters.SinWave=1;
        Parameters.Wavelength=6;
        % Parameters.FiberDensity=0.9;
        Parameters.FiberDensity='loose';
        Parameters.FiberLength=64;
        % Parameters.AngleSigma=1.5;
        % Parameters.AngleHighSigma=1.5;
        Parameters.AngleMean='progressive';
        Parameters.AngleSD=20;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=0;
        Parameters.AngleRangeLimiter='exp';
        Parameters.FieldCoalescence=0;
        Parameters.FiberBlendingMode='add';
        Parameters.FiberBlendingExponent=2;
    case 'V Uncoherent Wavelets'
        Parameters.OctaveIntensity=sqrt([1 1]);
        Parameters.SinWave=0;
        % Parameters.FiberDensity=0.9;
        Parameters.FiberDensity='loose';
        Parameters.FiberLength=64;
        Parameters.AngleSigma=1.5;
        Parameters.AngleHighSigma=1.5;
        Parameters.AngleMean='progressive';
        Parameters.AngleSD=20;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=0;
        Parameters.AngleRangeLimiter='exp';
        Parameters.FieldCoalescence=0;
        Parameters.FiberBlendingMode='add';
        Parameters.FiberBlendingExponent=2;
    case 'VI Entangled'
        Parameters.OctaveIntensity=sqrt([1 1 1]);
        Parameters.SinWave=1;
        Parameters.Wavelength=29.01;
        % Parameters.FiberDensity=0.7;
        Parameters.FiberDensity='loose';
        Parameters.FiberLength=256;
        % Parameters.AngleSigma=7.2533;
        % Parameters.AngleHighSigma=7.2533;
        Parameters.AngleMean='progressive';
        Parameters.AngleSD=20;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=100;
        Parameters.AngleRangeLimiter='exp';
        Parameters.FieldCoalescence=0;
        Parameters.FiberBlendingMode='add';
        Parameters.FiberBlendingExponent=2;
    case 'VII Chaotic'
        Parameters.OctaveIntensity=sqrt([1 1]);
        Parameters.SinWave=0;
        % Parameters.FiberDensity=1.5;
        Parameters.FiberDensity='loose';
        Parameters.FiberLength=256;
        Parameters.AngleSigma=10;
        Parameters.AngleHighSigma=inf;
        Parameters.AngleMean='progressive';
        Parameters.AngleSD=70;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=100;
        Parameters.AngleRangeLimiter='none';
        Parameters.FieldCoalescence=0;
        Parameters.FiberBlendingMode='add';
        Parameters.FiberBlendingExponent=2;
    case 'VIII Composite'
        Parameters.OctaveIntensity=sqrt([1 1 1]);
        Parameters.SinWave=0;
        % Parameters.FiberDensity=0.7;
        Parameters.FiberDensity='loose';
        Parameters.FiberLength=256;
        Parameters.AngleSigma=50;
        Parameters.AngleHighSigma=inf;
        Parameters.AngleMean='progressive';
        Parameters.AngleSD=50;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=0;
        Parameters.AngleRangeLimiter='exp';
        Parameters.FieldCoalescence=0;
        Parameters.FiberBlendingMode='add';
        Parameters.FiberBlendingExponent=2;
    case 'IX Multiscale Waves'
        Parameters.OctaveIntensity=sqrt([1 1 1]);
        Parameters.SinWave=0;
        % Parameters.FiberDensity=0.7;
        Parameters.FiberDensity='loose';
        Parameters.FiberLength=128;
        Parameters.FiberLengthDistribution='Poisson';
        Parameters.AngleSigma=4;
        Parameters.AngleHighSigma=4;
        Parameters.AngleMean='progressive';
        Parameters.AngleSD=40;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=0;
        Parameters.AngleRangeLimiter='exp';
        Parameters.FieldCoalescence=0;
        Parameters.FiberBlendingMode='add';
        Parameters.FiberBlendingExponent=2;
        Parameters.AngleFieldOctaves=[1 2 4 8];
    case 'X Fuzzy Multiscale Waves'
        Parameters.HRGeneration='on';
        Parameters.OctaveIntensity=sqrt([1 1 1 1]);
        Parameters.SinWave=0;
        % Parameters.FiberDensity=0.7;
        Parameters.FiberDensity='loose';
        Parameters.FiberLength=[128 32 16 8];
        Parameters.AngleSigma=0.5;
        Parameters.AngleHighSigma=0.5;
        Parameters.AngleMean='progressive';
        Parameters.AngleSD=30;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=0;
        Parameters.AngleRangeLimiter='exp';
        Parameters.FieldCoalescence=0;
        Parameters.FiberBlendingMode='add';
        Parameters.FiberBlendingExponent=2;
        Parameters.AngleFieldOctaves=[1 2 4 8 16 32];
        Parameters.FiberLengthDistribution='none';
        Parameters.FiberLengthIntensity=[3 1 1 1];
        Parameters.FiberLengthRescale='off';
        case 'Braids'
        Parameters.OctaveIntensity=sqrt([1 1 1]);
        Parameters.FiberDensity=0.7;
        Parameters.FiberLength=100;
        Parameters.AngleSigma=5;
        Parameters.AngleHighSigma=10;
        Parameters.AngleMean='random';
        Parameters.AngleSD=40;
        Parameters.AngleJitter=3;
        Parameters.AngleSpatialJitter=4;
        Parameters.AngleRangeLimiter='exp';
        Parameters.FieldCoalescence=0.7;
    case 'Whirlpools'
        Parameters.OctaveIntensity=sqrt([1]);
        Parameters.FiberDensity=1;
        Parameters.FiberLength=70;
        Parameters.AngleSigma=5;
        Parameters.AngleHighSigma=10;
        Parameters.AngleMean='random';
        Parameters.AngleSD=45;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=2;
        Parameters.AngleRangeLimiter='none';
        Parameters.FieldCoalescence=0.2;
    case 'Nodes'
        Parameters.OctaveIntensity=sqrt([1 1]);
        Parameters.FiberDensity=0.5;
        Parameters.FiberLength=100;
        Parameters.AngleSigma=10;
        Parameters.AngleHighSigma=15;
        Parameters.AngleMean='random';
        Parameters.AngleSD=45;
        Parameters.AngleJitter=0;
        Parameters.AngleSpatialJitter=1;
        Parameters.AngleRangeLimiter='none';
        Parameters.FieldCoalescence=0;
    case 'Curls'
        Parameters.OctaveIntensity=sqrt([1 1]);
        Parameters.FiberDensity=0.8;
        Parameters.FiberLength=100;
        Parameters.AngleSigma=10;
        Parameters.AngleHighSigma=10;
        Parameters.AngleMean='random';
        Parameters.AngleSD=50;
        Parameters.AngleJitter=5;
        Parameters.AngleSpatialJitter=1;
        Parameters.AngleRangeLimiter='none';
        Parameters.FieldCoalescence=0.2;
    otherwise
        Parameters=[];
        return
end

% force custom parameters on top of the presets
if not(isempty(CustomParameters))
    fl=fields(CustomParameters);
    for i=1:numel(fl)
        Parameters.(fl{i})=CustomParameters.(fl{i});
    end
end
