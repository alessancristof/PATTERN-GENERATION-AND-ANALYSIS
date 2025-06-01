function [I,GenerativeField,Parameters]=syntheticImageGeneration(Parameters,Size,SliceNumber,Spacing,Class,CustomParameters,flag_figure)
% Main function for generating synthetic patterns as 3D image stack (I) and
% the corresponding fiber direction field (A). Generation depends on parameters
% fed by the input structure variable Parameters or on a named preset. 

% INPUTS:
% If there is no input, presets and parameters are queried by modal dialog
% windows. Otherwise:
% Parameters (input): structure defining the parameters for image stack syntehsis
% (see function patternPreset for explanation)
% or char defining the name from some presets ('I Straight Lines','II Coherent Waves'
% ,'III Uncoherent Waves', 'IV Coherent Wavelets', 'V Uncoherent Wavelets',
% 'VI Entangled','VII Chaotic', 'VIII Composite','IX Multiscale Waves',
% 'X Fuzzy Multiscale Waves').
% Size (input): size of each slice
% SliceNumber (input): number of slices in the stack
% Spacing (optional input): 1x3 vector with pixel spacing along x and y and slice spacing
% along z, define spatial units for processing, default is [0.34 0.34 1].
% Class (optional input): char defining the numerical class of output,
% 'uint16' (unsigned 16 bit integer) or 'double' (64 bit floating point integer) 
% CustomParameters (optional input): structure with arbitrary few custom parameter values 
% to overwrite some produced by the preset. 
% flag_figure (optional input): binary flag for showing the resulting image

%OUTPUTS:
% I (output): n x m x s matrix, the synthetic image stack
% GenerativeField (output): n x m x s x 2 matrix, 2D vectorial field defining the
% direction for fiber drawing in the image stack
% Parameters (output): updated full parameters for image generation


% NOTES:

% Correspondence with article: : I. Straight Lines = 'I Straight Lines',
% II. Waves = 'III Uncoherent Waves', III. Multiscale Waves = 'IX Multiscale Waves',
% IV. Wavelets = 'V Uncoherent Wavelets', V. Incoherent waves = 'VI Entangled',
% VI. Chaotic Fibers = 'VII Chaotic', and VII. Realistic Fibers = 'X Fuzzy Multiscale Waves'. 

% How to modify the modulation and the no signal background by manual input
% Input structure CustomParameters with fields:
% TissuePercentage: percentage of full intensity tissue in the image after
% modulation (100 % for no modulation)
% BackgroundPercentage: percentage of zero intensity tissue in the image after
% modulation (0 % for no background);
% BackgroundSigma: sigma (in spatial units) for the Gaussian blur that defines the scale of intensity modulation 
% modulation

% How to add noise and an intensity baseline (dark current) by manual input
% Input structure CustomParameters with fields:
% SNRShot: SNR for "shot" noise (inf for no noise)
% SNRMult: SNR for "multiplicative" noise (inf for no noise)
% SNRWhite: SNR for base white noise (inf for no noise)
% SDR=10: "average" signal to "Dark Current" Ratio.

if nargin==0
    list={'I Straight Lines','II Coherent Waves','III Uncoherent Waves', 'IV Coherent Wavelets', 'V Uncoherent Wavelets','VI Entangled','VII Chaotic', 'VIII Composite','IX Multiscale Waves','X Fuzzy Multiscale Waves','Braids','Whirlpools','Nodes','Curls'};
    [indx,tf] = listdlg('Name','Pattern Presets','PromptString','Select a pattern preset.','ListSize',[250 150],'SelectionMode','single','ListString',list,'InitialValue',1);
    if not(tf)
        return
    end
    pattern_preset=list{indx};
    Parameters=patternPreset(pattern_preset);
    fieldnames=fields(Parameters);
    def={};
    for i=1:numel(fieldnames)
        var=Parameters.(fieldnames{i});
        flag_char(i)=ischar(var);
        if flag_char(i)
            def=cat(1,def,{var});
        else
            def=cat(1,def,{num2str(var)});
        end
    end
    fieldnames1=fieldnames(1:15);
    fieldnames2=fieldnames(16:26);
    fieldnames3=fieldnames(27:40);
    def1=def(1:15);
    def2=def(16:26);
    def3=def(27:40);
    
    answer1 = inputdlg(fieldnames1,'Stack and Noise Properties',ones(numel(fieldnames1),1)*[0.9 100],def1);
    answer2 = inputdlg(fieldnames2,'Directional Field Properties',ones(numel(fieldnames2),1)*[0.9 100],def2);
    answer3 = inputdlg(fieldnames3,'Fiber Drawing Properties',ones(numel(fieldnames3),1)*[0.9 100],def3);
    answer=cat(1,answer1,answer2,answer3);
    for i=1:numel(answer)
        var=str2num(answer{i});
        if isempty(var)
            Parameters.(fieldnames{i})=answer{i};
        else
            Parameters.(fieldnames{i})=var;
        end
        % if flag_char(i)
        %     Parameters.(fieldnames{i})=answer{i};
        % else
        %     Parameters.(fieldnames{i})=str2num(answer{i});
        % end
    end
    flag_figure=true;
    Class=Parameters.Class;
else
    if nargin<4
        Spacing=[];
    end
    if isempty(Spacing)
        Spacing=[0.34 0.34 1];
    end

    if nargin<5
        Class='uint16';
    end

    if nargin<6
        CustomParameters=[];
    end
    if nargin<7
        flag_figure=false;
    end

    if ischar(Parameters)
        Parameters=patternPreset(Parameters,CustomParameters);
        if isempty(Parameters)
            msgbox('Not valid preset name.')
            return
        end
    end
    Parameters.Size=Size;
    Parameters.SliceNumber=SliceNumber;
    Parameters.Spacing=Spacing;
    Parameters.Class=Class;
end

Parameters.Synthetic=true;

if strcmp(Parameters.AngleMean,'progressive') % mean angle increases with slices
    flag_progressive=true;
else
    flag_progressive=false;
end


I=zeros([Parameters.Size Parameters.SliceNumber]); % the syntethic image
GenerativeField=zeros([Parameters.Size Parameters.SliceNumber 2]); % the vectorial field that generates it
for i=1:Parameters.SliceNumber
    if flag_progressive
        Parameters.AngleMean=(i-1)*180/Parameters.SliceNumber-90;
    end

    % generate the slice
    [Ib,Parametersb,Ab]=syntheticImage(Parameters);
    I(:,:,i)=Ib; 
    GenerativeField(:,:,i,1)=Ab(:,:,1); 
    GenerativeField(:,:,i,2)=Ab(:,:,2);

    % add in Parameters extra descriptive values resulting from the
    % generation of each slice for accuracy evaluation purpose (the real angle mean and SD before and after the intensity modulation, the theoretical angle)
    if i==1
        Parameters.DarkCurrent=Parametersb.DarkCurrent;
        Parameters.WhiteNoiseCoefficient=Parametersb.WhiteNoiseCoefficient;
        Parameters.ShotNoiseCoefficient=Parametersb.ShotNoiseCoefficient;
        Parameters.MultNoiseCoefficient=Parametersb.MultNoiseCoefficient;

        Parameters.RealAngleMean=Parametersb.RealAngleMean;
        Parameters.RealAngleSD=Parametersb.RealAngleSD;
        Parameters.PostModulationAngleMean=Parametersb.PostModulationAngleMean;
        Parameters.PostModulationAngleSD=Parametersb.PostModulationAngleSD;
        Parameters.TheoreticalAngle=Parametersb.TheoreticalAngle;
    else
        Parameters.RealAngleMean=cat(1,Parameters.RealAngleMean,Parametersb.RealAngleMean);
        Parameters.RealAngleSD=cat(1,Parameters.RealAngleSD,Parametersb.RealAngleSD);
        Parameters.PostModulationAngleMean=cat(1,Parameters.PostModulationAngleMean,Parametersb.PostModulationAngleMean);
        Parameters.PostModulationAngleSD=cat(1,Parameters.PostModulationAngleSD,Parametersb.PostModulationAngleSD);
        Parameters.TheoreticalAngle=cat(1,Parameters.TheoreticalAngle,Parametersb.TheoreticalAngle);
    end
end



if flag_figure
    I2=uint8(I/prctile(I,99.9,"all")*255);
    if size(I,3)==1
        figure('Name','Synthetic Pattern');
        imshow(I2);
    else
        figure('Name','Synthetic Pattern Montage');
        montage(I2);
    end
end


% turn output in uint16
if strcmp(Class,'uint16')
    I=uint16(I);
    Asq=GenerativeField(:,:,:,1).^2+GenerativeField(:,:,:,2).^2;
    gmax=sqrt(max(Asq,[],"all"));
    GenerativeField=int16(GenerativeField/gmax*32767);
end

