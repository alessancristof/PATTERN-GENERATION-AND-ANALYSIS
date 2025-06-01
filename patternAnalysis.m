function [TImageProperties,TAnalysisSetting,TNoiseTot,TGradientTot,TGradientHistoTot,TStructureTot,TStructureHistoTot,TSyntheticAngle,TGenerativeFieldTot]=patternAnalysis(I,IParameters,NEParametersIn,BDParametersIn,STParametersIn,GenerativeField,FigureSliceNumber,FlagSave)

% INPUT
% 1) I : n x m x s matrix, stack of s images of n x m size, the pattern to be analyzed.
% 2) IParameters: parameters that define image I characteristics necessary for
% processing. If the image is synthetic, it can be a structure containing all the generative
% parameters, which triggers a validation. Otherwise it is a 1 x 3 vector defining Spacing directly. If empty Spacing is
% assumed unitary.
% If NEParametersIn, BDParametersIn, and STParametersIn are missing or empty, the processing
% settings are queried by modal dialog windows, otherwise they can be given to
% the function as arguments.
% 3) NEParametersIn (optional): structure containing the custom parameters for noise
% estimation, look in function noiseEstimator for an explanation of its fields.
% 4) BDParametersIn (optional): structure containing the custom parameters for
% background detection, in particular the field SNRThreshold to set the SNR cutoff. Look in function backgroundDetection for an
% explanation of its fields. 
% 4) STParametersIn (optional): structure containing the custom parameters for
% background detection, in particular the field GradientSNRThreshold to set
% the gradient SNR cutoff used in gradient analysis, look in function structureTensorProcessing for an explanation of its fields.
% 5) GenerativeField (optional): n x m x s x 2 matrix representing the
% 2D vectorial field used to generate the pattern
% 6) FigureSliceNumber: if not 0, generates figures with that slice index as
% example.
% 7) FlagSave: if true save the results in an excel file


if nargin<3
    NEParametersIn=[];
end
if nargin<4
    BDParametersIn=[];
end
if nargin<5
    STParametersIn=[];
end
if and(and(isempty(NEParametersIn),isempty(BDParametersIn)),isempty(STParametersIn))
    flag_inputdlg=true;
end

if nargin<6
    GenerativeField=[];
end
if nargin<7
    FigureSliceNumber=1;
end
if nargin<8
    FlagSave=false;
end


% determine if the image is synthetic
flag_synthetic=false;
if not(isempty(IParameters))
    if isfield(IParameters,'Synthetic')
        if IParameters.Synthetic
            flag_synthetic=true;
        end
    end
end

% retrieve voxel spacing if given in IParameters to define a spatial scale
% of the image
Spacing=[1 1 1];
if not(isempty(IParameters))
    if isstruct(IParameters)
        if isfield(IParameters,'Spacing')
            Spacing=IParameters.Spacing;
        end
    elseif numel(IParameters)==3
        Spacing=IParameters;
    end
end

% Default parameters values for noise estimation
NEParameters.HighpassWavelength=1.2;
NEParameters.PatchSize=10;
NEParameters.BinSize=100;
% if the image is tagged as synthetic performs a per slice noise estimation
% and compare the results with the generative values.
if flag_synthetic
    NEParameters.SeparateSlices=true; 
else
    NEParameters.SeparateSlices=false;
end
NEParameters.VartotCutoffPrc=50;
NEParameters.VartotTrimPrc=20;
NEParameters.VartotLowPrc=2;
NEParameters.DarkCurrentTrimPrc=NEParameters.VartotTrimPrc;
NEParameters.ICutoffPrc=99.5;
NEParameters.VarTrimPrc=NEParameters.VartotTrimPrc;
if flag_inputdlg
    NEParameters.ExtraFigureSliceNumber=FigureSliceNumber;
else
    NEParameters.ExtraFigureSliceNumber=0;
end

% defaults parameters for background detection (areas with low signal SNR, possibly not tissue)
BDParameters.SignalSmoothingRadius=1; % typically the same as StructureTensorRadius
BDParameters.Smoothing3D=true;
BDParameters.SNRThreshold=3; % Threshold for image intesity SNR.
BDParameters.LowSignalThreshold=1;
BDParameters.MaskSmoothingRadius=0.5;
BDParameters.AreaOpeningThreshold=100;

if flag_inputdlg % if input dialog windows are activated add these extras settings to ask for
    STParameters.FigureSliceNumber=FigureSliceNumber;
    STParameters.FlagSave=FlagSave;
end
% defaults parameters for gradient and structure generation
STParameters.GradientDiscretization='Sharr Tap 5';
STParameters.GradientBiasRemoval=true;
STParameters.GradientNormalization=false;
STParameters.StructureTensorRadius=1;
STParameters.StructureTensorSmoothingDirection='2D';
STParameters.StructureBiasRemoval=true;
% defaults parameters for gradient and structure analysis
STParameters.SliceNumberInterval=[]; % 1x2 vector, interval of slices to be analyzed, if empty they are detected automatically.
STParameters.SNIBackgroundThreshold=50; % threshold on percentage of signal background to automatically detect the slices to be analyzed.
STParameters.GradientSNRThreshold=2; % Threshold for gradient SNR
STParameters.StructureSNRThreshold=4; % Threshold for structure SNR
STParameters.AnisotropyDescriptor='DA'; % Type of anisotropy measure for analysis
STParameters.AnisotropyThreshold=0.8; % Threshold to isolate high anisotropy areas
STParameters.StructureWeightType='First Variation'; % What to be used as weights for structure analysis. 
STParameters.CropSize=3; % number of pixels to crop from the image border prior to analysis (to avoid border artifacts affecting the analysis).

% produces input dialog windows to interact with the custom parameters
if flag_inputdlg
    fieldnames=fields(NEParameters);
    def={};
    for i=1:numel(fieldnames)
        var=NEParameters.(fieldnames{i});
        flag_char(i)=ischar(var);
        if flag_char(i)
            def=cat(1,def,{var});
        else
            def=cat(1,def,{num2str(var)});
        end
    end
    answer = inputdlg(fieldnames,'Noise Estimation Settings',ones(numel(fieldnames),1)*[0.9 100],def);
    for i=1:numel(answer)
        if flag_char(i)
            NEParameters.(fieldnames{i})=answer{i};
        else
            NEParameters.(fieldnames{i})=str2num(answer{i});
        end
    end

    flag_char=[];
    fieldnames=fields(BDParameters);
    def={};
    for i=1:numel(fieldnames)
        var=BDParameters.(fieldnames{i});
        flag_char(i)=ischar(var);
        if flag_char(i)
            def=cat(1,def,{var});
        else
            def=cat(1,def,{num2str(var)});
        end
    end
    answer = inputdlg(fieldnames,'Background Detection Settings',ones(numel(fieldnames),1)*[0.9 100],def);
    for i=1:numel(answer)
        if flag_char(i)
            BDParameters.(fieldnames{i})=answer{i};
        else
            BDParameters.(fieldnames{i})=str2num(answer{i});
        end
    end

    flag_char=[];
    fieldnames=fields(STParameters);
    def={};
    for i=1:numel(fieldnames)
        var=STParameters.(fieldnames{i});
        flag_char(i)=ischar(var);
        if flag_char(i)
            def=cat(1,def,{var});
        else
            def=cat(1,def,{num2str(var)});
        end
    end
    answer = inputdlg(fieldnames,'Structure Tensor Analysis Settings',ones(numel(fieldnames),1)*[0.9 100],def);
    for i=1:numel(answer)
        if flag_char(i)
            STParameters.(fieldnames{i})=answer{i};
        else
            STParameters.(fieldnames{i})=str2num(answer{i});
        end
    end
    FigureSliceNumber=STParameters.FigureSliceNumber;
    FlagSave=STParameters.FlagSave;
end

% overrides the defaults parameters with the manual input values
if not(isempty(NEParametersIn))
    fl=fields(NEParametersIn);
    for i=1:numel(fl)
        NEParameters.(fl{i})=NEParametersIn.(fl{i});
    end
end
if not(isempty(BDParametersIn))
    fl=fields(BDParametersIn);
    for i=1:numel(fl)
        BDParameters.(fl{i})=BDParametersIn.(fl{i});
    end
end
if not(isempty(STParametersIn))
    fl=fields(STParametersIn);
    for i=1:numel(fl)
        STParameters.(fl{i})=STParametersIn.(fl{i});
    end
end

I=double(I);

[DarkCurrent,DarkVariance,ShotCoefficient,MultiplicativeCoefficient,NEParameters]=noiseEstimator(I,NEParameters,IParameters,FigureSliceNumber>0);
% if noise estimation was performed per slice, consider the average coefficients over
% all slices
if numel(DarkCurrent)>1 
    DarkCurrent=mean(DarkCurrent);
    DarkVariance=mean(DarkVariance);
    ShotCoefficient=mean(ShotCoefficient);
    MultiplicativeCoefficient=mean(MultiplicativeCoefficient);
end

Isig=max(I-DarkCurrent,0);
NoiseVariance=DarkVariance+ShotCoefficient*Isig+MultiplicativeCoefficient*Isig.^2;

% compute the signal background (areas with low intensity signal SNR).
[SignalBackgroundMask,BDParameters,SNR]=backgroundDetection(Isig,NoiseVariance,Spacing,BDParameters);

% compute gradient and diagonalized structure tensor 
[FirstVariation,SecondVariation,FirstEigenvector1,FirstEigenvector2,StructureSNR,G1,G2,GradientSNR,Graw1,Graw2,G1var,G2var]=structureTensorProcessing(I,NoiseVariance,Spacing,STParameters);


% detect a background for gradient (areas with low gradient SNR).
GradientBackgroundMask=GradientSNR<STParameters.GradientSNRThreshold;

% detect a background for structure, where there is no appreciable structure (area with low SNR for the
% square-rooted trace of the structure tensor).
StructureBackgroundMask=StructureSNR<STParameters.StructureSNRThreshold;

% detect areas with high pattern anisotropy by thresholding an
% anisotropy measure.
Aniso=structureDescriptors(FirstVariation,SecondVariation,FirstEigenvector1,FirstEigenvector2,{STParameters.AnisotropyDescriptor});
AnisotropyMask=Aniso>STParameters.AnisotropyThreshold;
AnisotropyMask(SignalBackgroundMask)=0;
AnisotropyMask(StructureBackgroundMask)=0;

% compute percentage of signal background per slice
for i=1:size(SignalBackgroundMask,3)
    var=SignalBackgroundMask(:,:,i);
    SignalBackgroundPercent(i)=sum(var,"all")/numel(var)*100;
end

% create a label image for the gradient (0 = signal background, 1 =
% gradient background, 2 = gradient is present).
GradientClass=ones(size(I),'uint8')*2;
GradientClass(GradientBackgroundMask)=1;
GradientClass(SignalBackgroundMask)=0;

% create a label image for the structure (0 = signal background, 1 =
% structure background, 2 = structure is present but isotropic, 3 = strucure is present and anisotropic).
StructureClass=ones(size(I),'uint8')*2;
StructureClass(AnisotropyMask)=3;
StructureClass(StructureBackgroundMask)=1;
StructureClass(SignalBackgroundMask)=0;

clear SignalBackgroundMask GradientBackgroundMask StructureBackgroundMask AnisotropyMask

% create figures for exemplative slice
if FigureSliceNumber>0
    overlay_fact=0.3;

    figure('Name',['Slice ' num2str(FigureSliceNumber) ': Gradient (red), Background (blue), and Gradient Background (cyan). Zoom in for detail.'])
    R=Isig(:,:,FigureSliceNumber);
    R=R/prctile(R,99.9,"all")*255;
    G=R;
    B=R;
    mask=GradientClass(:,:,FigureSliceNumber)==0;
    col=[0 0 255];
    R(mask)=(1-overlay_fact)*R(mask)+overlay_fact*col(1);
    G(mask)=(1-overlay_fact)*G(mask)+overlay_fact*col(2);
    B(mask)=(1-overlay_fact)*B(mask)+overlay_fact*col(3);
    mask=GradientClass(:,:,FigureSliceNumber)==1;
    col=[0 255 255];
    R(mask)=(1-overlay_fact)*R(mask)+overlay_fact*col(1);
    G(mask)=(1-overlay_fact)*G(mask)+overlay_fact*col(2);
    B(mask)=(1-overlay_fact)*B(mask)+overlay_fact*col(3);
    I=uint8(cat(3,R,G,B));
    imshow(I);
    hold on
    quiver(G2(:,:,FigureSliceNumber),G1(:,:,FigureSliceNumber),2,'r','Alignment','center');
    hold off

    [EVfirst, EVsecond]=structureDescriptors(FirstVariation,SecondVariation,FirstEigenvector1,FirstEigenvector2,{'First Eigenvector', 'Second Eigenvector'});
    R=Isig(:,:,FigureSliceNumber);
    R=R/prctile(R,99.9,"all")*255;
    G=R;
    B=R;
    mask=GradientClass(:,:,FigureSliceNumber)==0;
    col=[0 0 255];
    R(mask)=(1-overlay_fact)*R(mask)+overlay_fact*col(1);
    G(mask)=(1-overlay_fact)*G(mask)+overlay_fact*col(2);
    B(mask)=(1-overlay_fact)*B(mask)+overlay_fact*col(3);
    mask=StructureClass(:,:,FigureSliceNumber)==1;
    col=[255 0 0];
    R(mask)=(1-overlay_fact)*R(mask)+overlay_fact*col(1);
    G(mask)=(1-overlay_fact)*G(mask)+overlay_fact*col(2);
    B(mask)=(1-overlay_fact)*B(mask)+overlay_fact*col(3);
    mask=StructureClass(:,:,FigureSliceNumber)==3;
    col=[255 255 0];
    R(mask)=(1-overlay_fact)*R(mask)+overlay_fact*col(1);
    G(mask)=(1-overlay_fact)*G(mask)+overlay_fact*col(2);
    B(mask)=(1-overlay_fact)*B(mask)+overlay_fact*col(3);
    I=uint8(cat(3,R,G,B));
    var1=FirstVariation(:,:,FigureSliceNumber)-SecondVariation(:,:,FigureSliceNumber);
    figure('Name',['Slice ' num2str(FigureSliceNumber) ': First Eigenvector (red) weighted by variation difference, Background (blue), Low Structure (red), and ' STParameters.AnisotropyDescriptor '>' num2str(STParameters.AnisotropyThreshold) ' (yellow). Zoom in for detail.'])
    imshow(I);
    hold on
    quiver(var1.*EVfirst(:,:,FigureSliceNumber,2),var1.*EVfirst(:,:,FigureSliceNumber,1),2,'r','Alignment','center','ShowArrowHead','off');
    hold off
    figure('Name',['Slice ' num2str(FigureSliceNumber) ': Second Eigenvector (green, fiber direction) weighted by variation difference, Background (blue), Low Structure (red), and ' STParameters.AnisotropyDescriptor '>' num2str(STParameters.AnisotropyThreshold) ' (yellow). Zoom in for detail.'])
    imshow(I);
    hold on
    quiver(var1.*EVsecond(:,:,FigureSliceNumber,2),var1.*EVsecond(:,:,FigureSliceNumber,1),2,'g','Alignment','center','ShowArrowHead','off');
    hold off

    figure('Name',['Slice ' num2str(FigureSliceNumber) ': First (red channel) and Second (green channel) Variations'])
    R=FirstVariation(:,:,FigureSliceNumber);
    m=prctile(R,99.9,"all");
    R=R/m*255;
    G=SecondVariation(:,:,FigureSliceNumber)/m*255;
    B=zeros(size(R));
    I=uint8(cat(3,R,G,B));
    imshow(I);

    Aniso(GradientClass==0)=0;
    Aniso(StructureClass==1)=0;
    figure('Name',['Slice ' num2str(FigureSliceNumber) ': ' STParameters.AnisotropyDescriptor ' in foreground areas'])
    imshow(Aniso(:,:,FigureSliceNumber),[]);
    
end


% isolate the slices to be analyzed
if isempty(STParameters.SliceNumberInterval)
    start_slice=[];
    end_slice=[];
    for i=1:numel(SignalBackgroundPercent)
        if SignalBackgroundPercent(i)<STParameters.SNIBackgroundThreshold
            start_slice=i;
            break
        end
    end
    for i=numel(SignalBackgroundPercent):-1:1
        if SignalBackgroundPercent(i)<STParameters.SNIBackgroundThreshold
            end_slice=i;
            break
        end
    end
    if and(not(isempty(start_slice)),not(isempty(end_slice)))
        STParameters.SliceNumberInterval=[start_slice end_slice];
    end
end

if isempty(STParameters.SliceNumberInterval)
    return
end

% crop the border of images
if STParameters.CropSize>0
    GradientClass=StructureClass((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    StructureClass=StructureClass((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    G1=G1((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    G2=G2((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    FirstVariation=FirstVariation((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    SecondVariation=SecondVariation((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    FirstEigenvector1=FirstEigenvector1((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    FirstEigenvector2=FirstEigenvector2((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    Isig=Isig((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    NoiseVariance=NoiseVariance((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    if not(isempty(GenerativeField))
        GenerativeField=GenerativeField((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:,:);
    end
    Graw1=Graw1((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    Graw2=Graw2((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    G1var=G1var((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
    G2var=G2var((1+STParameters.CropSize):(end-STParameters.CropSize),(1+STParameters.CropSize):(end-STParameters.CropSize),:);
end

TGradientTot=[];
TGradientHistoTot=[];
TStructureTot=[];
TStructureHistoTot=[];
TNoiseTot=[];
TGenerativeFieldTot=[];
SliceNumberReal=[];
for i=STParameters.SliceNumberInterval(1):STParameters.SliceNumberInterval(2)
    if isempty(GenerativeField)
        var=[];
    else
        var=GenerativeField(:,:,i,:);
    end
    [TGradient, TGradientHisto, GradientHisto, TStructure, TStructureHisto, StructureHisto, TNoise, TGenerativeField]=oneSliceStructureAnalysis(i,GradientClass(:,:,i),StructureClass(:,:,i),G1(:,:,i),G2(:,:,i),FirstVariation(:,:,i),SecondVariation(:,:,i),FirstEigenvector1(:,:,i),FirstEigenvector2(:,:,i),Isig(:,:,i),NoiseVariance(:,:,i),var,STParameters.StructureWeightType,STParameters.StructureBiasRemoval,Graw1(:,:,i),Graw2(:,:,i),G1var(:,:,i),G2var(:,:,i));
    if not(isempty(TGradient))
        SliceNumberReal=cat(1,SliceNumberReal,i);
        TGradientTot=cat(1,TGradientTot,TGradient);
        TGradientHistoTot=cat(1,TGradientHistoTot,TGradientHisto);
        TStructureTot=cat(1,TStructureTot,TStructure);
        TStructureHistoTot=cat(1,TStructureHistoTot,TStructureHisto);
        TNoiseTot=cat(1,TNoiseTot,TNoise);
        if not(isempty(TGenerativeField))
            TGenerativeFieldTot=cat(1,TGenerativeFieldTot,TGenerativeField);
        end
        if i==FigureSliceNumber
            oneSliceGraphicalResults(TGradient, TGradientHisto, GradientHisto, TStructure, TStructureHisto, StructureHisto,STParameters,FigureSliceNumber);
        end
    end
end

% if the image is synthetic, compute discrepancies between theoretical
% values and mesured values for validation
TSyntheticAngle=[];
if flag_synthetic
    if isfield(IParameters,'RealAngleMean')
        TheoreticalAngle=IParameters.TheoreticalAngle;
        TheoreticalAngleSD=IParameters.AngleSD*ones(size(TheoreticalAngle));
        RealAngleMean=IParameters.RealAngleMean(SliceNumberReal);
        RealAngleSD=IParameters.RealAngleSD(SliceNumberReal);
        PostModulationAngleMean=IParameters.PostModulationAngleMean(SliceNumberReal);
        PostModulationAngleSD=IParameters.PostModulationAngleSD(SliceNumberReal);
        varnames={'Slice Number', 'Theoretical Angle', 'Theoretical Angle SD', 'Field Angle Mean', 'Field Angle SD', 'Modulated Field Angle Mean', 'Modulated Field Angle SD'};
        TSyntheticAngle=table(num2cell(SliceNumberReal),num2cell(TheoreticalAngle),num2cell(TheoreticalAngleSD),num2cell(RealAngleMean),num2cell(RealAngleSD),num2cell(PostModulationAngleMean),num2cell(PostModulationAngleSD),'VariableNames',varnames);

        % comparison with measured angles
        TheoreticalModulatedAngleDifference=cell2mat(TSyntheticAngle.('Modulated Field Angle Mean'))-cell2mat(TSyntheticAngle.('Theoretical Angle'));
        GradientGlobalThAngleDifference=cell2mat(TGradientTot.('Global Angle'))-cell2mat(TSyntheticAngle.('Theoretical Angle')); % difference from the global gradient diagonalization angle
        GradientThAngleDifference=cell2mat(TGradientHistoTot.('Angle Mean'))-cell2mat(TSyntheticAngle.('Theoretical Angle'));
        StructureThAngleDifference=cell2mat(TStructureHistoTot.('Angle Mean'))-cell2mat(TSyntheticAngle.('Theoretical Angle'));
        GradientGlobalAngleDifference=cell2mat(TGradientTot.('Global Angle'))-cell2mat(TSyntheticAngle.('Modulated Field Angle Mean')); % % difference from the global gradient diagonalization angle
        GradientAngleDifference=cell2mat(TGradientHistoTot.('Angle Mean'))-cell2mat(TSyntheticAngle.('Modulated Field Angle Mean'));
        StructureAngleDifference=cell2mat(TStructureHistoTot.('Angle Mean'))-cell2mat(TSyntheticAngle.('Modulated Field Angle Mean'));
        varnames={'Theoretical - Modulated Angle Difference', 'Gradient Global - Theoretical Angle Difference', 'Gradient - Theoretical Angle Difference', 'Structure - Theoretical Angle Difference', 'Gradient Global - Modulated Angle Difference', 'Gradient - Modulated Angle Difference', 'Structure - Modulated Angle Difference'};
        T2=table(num2cell(TheoreticalModulatedAngleDifference),num2cell(GradientGlobalThAngleDifference),num2cell(GradientThAngleDifference),num2cell(StructureThAngleDifference),num2cell(GradientGlobalAngleDifference),num2cell(GradientAngleDifference),num2cell(StructureAngleDifference),'VariableNames',varnames);
        TSyntheticAngle=cat(2,TSyntheticAngle,T2);

        % comparison with measured spread
        TheoreticalModulatedStdDifference=cell2mat(TSyntheticAngle.('Modulated Field Angle SD'))-cell2mat(TSyntheticAngle.('Theoretical Angle SD'));
        GradientThStdDifference=cell2mat(TGradientHistoTot.('Angle Std'))-cell2mat(TSyntheticAngle.('Theoretical Angle SD'));
        StructureThStdDifference=cell2mat(TStructureHistoTot.('Angle Std'))-cell2mat(TSyntheticAngle.('Theoretical Angle SD'));
        GradientStdDifference=cell2mat(TGradientHistoTot.('Angle Std'))-cell2mat(TSyntheticAngle.('Modulated Field Angle SD'));
        StructureStdDifference=cell2mat(TStructureHistoTot.('Angle Std'))-cell2mat(TSyntheticAngle.('Modulated Field Angle SD'));
        varnames={'Theoretical - Modulated Angle SD Difference', 'Gradient - Theoretical Angle SD Difference', 'Structure - Theoretical Angle SD Difference', 'Gradient - Modulated Angle SD Difference', 'Structure - Modulated Angle SD Difference'};
        T2=table(num2cell(TheoreticalModulatedStdDifference),num2cell(GradientThStdDifference),num2cell(StructureThStdDifference),num2cell(GradientStdDifference),num2cell(StructureStdDifference),'VariableNames',varnames);
        TSyntheticAngle=cat(2,TSyntheticAngle,T2);

        flag_ang=logical([0 1 0 1 0 1 0 1 1 1 1 1 1 1 0 0 0 0 0]);
        flag_secmom=logical([0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1]);
        T=table({'Mean over Slices'; 'Std over Slices'; 'Sqrt Second Moment over Slices'},'VariableNames',TSyntheticAngle.Properties.VariableNames(1));
        for i=2:size(TSyntheticAngle,2)
            if isa(TSyntheticAngle{:,i},'cell')
                var=cell2mat(TSyntheticAngle{:,i});
            else
                var=TSyntheticAngle{:,i};
            end
            if flag_ang(i)
                [me,sd]=circular_stat_ang(var,180);
            else
                me=mean(var);
                sd=std(var);
            end
            if flag_secmom(i)
                ssm=sqrt(me^2+sd^2);
            else
                ssm=nan;
            end
            T2=table({me; sd; ssm},'VariableNames',TSyntheticAngle.Properties.VariableNames(i));
            T=cat(2,T,T2);
        end
        TSyntheticAngle=cat(1,TSyntheticAngle,T);
    end
end


% compute the mean and the SD over the slices for all
% tables-----------------------------------------------------------------

T=table({'Mean over Slices'; 'Std over Slices'},'VariableNames',TNoiseTot.Properties.VariableNames(1));
for i=2:size(TNoiseTot,2)
    var=cell2mat(TNoiseTot{:,i});
    T2=table({mean(var); std(var)},'VariableNames',TNoiseTot.Properties.VariableNames(i));
    T=cat(2,T,T2);
end
TNoiseTot=cat(1,TNoiseTot,T);

T=table({'Mean over Slices'; 'Std over Slices'},'VariableNames',TGradientTot.Properties.VariableNames(1));
for i=2:size(TGradientTot,2)
    var=cell2mat(TGradientTot{:,i});
    if i==10
        [me,sd]=circular_stat_ang(var,180);
        T2=table({me; sd},'VariableNames',TGradientTot.Properties.VariableNames(i));
    else
        T2=table({mean(var); std(var)},'VariableNames',TGradientTot.Properties.VariableNames(i));
    end
    T=cat(2,T,T2);
end
TGradientTot=cat(1,TGradientTot,T);

T=table({'Mean over Slices'; 'Std over Slices'},'VariableNames',TStructureTot.Properties.VariableNames(1));
flag_stat=logical([0 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
flag_ang=logical([0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]);
for i=2:size(TStructureTot,2)
    var=cell2mat(TStructureTot{:,i});
    if flag_stat(i)
        if flag_ang(i)
            [me,sd]=circular_stat_ang(var,180);
            T2=table({me; sd},'VariableNames',TStructureTot.Properties.VariableNames(i));
        else
            T2=table({mean(var); std(var)},'VariableNames',TStructureTot.Properties.VariableNames(i));
        end
    else
        T2=table({[]; []},'VariableNames',TStructureTot.Properties.VariableNames(i));
    end
    T=cat(2,T,T2);
end
TStructureTot=cat(1,TStructureTot,T);

if not(isempty(TGenerativeFieldTot))
    T=table({'Mean over Slices'; 'Std over Slices'},'VariableNames',TGenerativeFieldTot.Properties.VariableNames(1));
    for i=2:size(TGenerativeFieldTot,2)
        var=cell2mat(TGenerativeFieldTot{:,i});
        T2=table({mean(var); std(var)},'VariableNames',TGenerativeFieldTot.Properties.VariableNames(i));
        T=cat(2,T,T2);
    end
    TGenerativeFieldTot=cat(1,TGenerativeFieldTot,T);
end

T=table({'Mean over Slices'; 'Std over Slices'},'VariableNames',TGradientHistoTot.Properties.VariableNames(1));
flag_ang=logical([0 1 0 1 0 1 0 0 0 0 0]);
for i=2:size(TGradientHistoTot,2)
    var=cell2mat(TGradientHistoTot{:,i});
    if flag_ang(i)
        [me,sd]=circular_stat_ang(var,180);
        T2=table({me; sd},'VariableNames',TGradientHistoTot.Properties.VariableNames(i));
    else
        T2=table({mean(var); std(var)},'VariableNames',TGradientHistoTot.Properties.VariableNames(i));
    end
    T=cat(2,T,T2);
end
TGradientHistoTot=cat(1,TGradientHistoTot,T);

T=table({'Mean over Slices'; 'Std over Slices'},'VariableNames',TStructureHistoTot.Properties.VariableNames(1));
flag_ang=logical([0 1 0 1 0 1 0 0 0 0 0]);
for i=2:size(TStructureHistoTot,2)
    var=cell2mat(TStructureHistoTot{:,i});
    if flag_ang(i)
        [me,sd]=circular_stat_ang(var,180);
        T2=table({me; sd},'VariableNames',TStructureHistoTot.Properties.VariableNames(i));
    else
        T2=table({mean(var); std(var)},'VariableNames',TStructureHistoTot.Properties.VariableNames(i));
    end
    T=cat(2,T,T2);
end
TStructureHistoTot=cat(1,TStructureHistoTot,T);

%------------------------------------------------------------------------


% Table with image properties
Tn={'Property', 'Value'};
if isstruct(IParameters)
    fl=fields(IParameters);
    for i=1:numel(fl)
        var={IParameters.(fl{i})};
        T2=table(fl(i),var,'VariableNames',Tn);
        if i==1
            TImageProperties=T2;
        else
            TImageProperties=cat(1,TImageProperties,T2);
        end
    end
else
    TImageProperties=table({'Size'; 'Spacing'},{size(I); Spacing},'VariableNames',Tn);
end

% Table with analysis settings
Tn={'Property', 'Value'};
% var={Spacing(1); Spacing(2); Spacing(3)};
% TAnalysisSetting=table({'Spacing 1'; 'Spacing 2'; 'Spacing 3'},var,'VariableNames',Tn);
TAnalysisSetting=table({'Spacing'},Spacing,'VariableNames',Tn);
fl=fields(NEParameters);
for i=1:numel(fl)
    var={NEParameters.(fl{i})};
    T2=table(fl(i),var,'VariableNames',Tn);
    TAnalysisSetting=cat(1,TAnalysisSetting,T2);
end
fl=fields(BDParameters);
for i=1:numel(fl)
    var={BDParameters.(fl{i})};
    T2=table(fl(i),var,'VariableNames',Tn);
    TAnalysisSetting=cat(1,TAnalysisSetting,T2);
end
fl=fields(STParameters);
for i=1:numel(fl)
    var={STParameters.(fl{i})};
    T2=table(fl(i),var,'VariableNames',Tn);
    TAnalysisSetting=cat(1,TAnalysisSetting,T2);
end


% save the results in an excel file
if FlagSave
    [file,path] = uiputfile('*.xlsx','Save the results to Excel file:','Analysis Results.xlsx');
else
    return
end
if file==0
    return
end

filename = fullfile(path,file);

writematrix("Image Properties",filename,'Sheet','Image Properties','Range','A1');
writetable(TImageProperties,filename,'Sheet','Image Properties','Range','A3');

writematrix("Processing Parameters",filename,'Sheet','Processing Parameters','Range','A1');
writetable(TAnalysisSetting,filename,'Sheet','Processing Parameters','Range','A3');

writematrix("Noise Descriptors",filename,'Sheet','Noise Descriptors','Range','A1');
writetable(TNoiseTot,filename,'Sheet','Noise Descriptors','Range','A3');

writematrix("Gradient Angular Distribution",filename,'Sheet','Gradient Angular Distribution','Range','A1');
writetable(TGradientHistoTot,filename,'Sheet','Gradient Angular Distribution','Range','A3');
writematrix("Structure Angular Distribution",filename,'Sheet','Structure Angular Distribution','Range','A1');
writetable(TStructureHistoTot,filename,'Sheet','Structure Angular Distribution','Range','A3');
writematrix("Signal and Global Structure Descriptors",filename,'Sheet','Global Structure Descriptors','Range','A1');
writetable(TGradientTot,filename,'Sheet','Global Structure Descriptors','Range','A3');
writematrix("Averaged Local Structure Descriptors",filename,'Sheet','Averaged Local Structure','Range','A1');
writetable(TStructureTot,filename,'Sheet','Averaged Local Structure','Range','A3');

if not(isempty(TSyntheticAngle))
    writematrix("Generative Angle Descriptors",filename,'Sheet','Generative Angle Descriptors','Range','A1');
    writetable(TSyntheticAngle,filename,'Sheet','Generative Angle Descriptors','Range','A3');
end
if not(isempty(TGenerativeFieldTot))
    writematrix("Generative Field Pixelwise RMS Errors",filename,'Sheet','Generative Field','Range','A1');
    writetable(TGenerativeFieldTot,filename,'Sheet','Generative Field','Range','A3');
end



%-----------------------------------------------------------------------------------------
function [TGradient, TGradientHisto, GradientHisto, TStructure, TStructureHisto, StructureHisto, TNoise, TGenerativeField]=oneSliceStructureAnalysis(SliceNumber,GradientClass,StructureClass,G1,G2,FirstVariation,SecondVariation,FirstEigenvector1,FirstEigenvector2,Signal,SignalNoiseVariance,GenerativeField,StructureWeightType,GlobalStructureBiasRemoval,Graw1,Graw2,G1var,G2var)


TGradient=[];
TGradientHisto=[];
GradientHisto=[];
TStructure=[];
TStructureHisto=[];
StructureHisto=[];
TNoise=[];
% TParentNoise=[];

TGenerativeField=[];

% compute percentage of signal background
SignalBackgroundMask=GradientClass==0;
SliceBackgroundPercentage=sum(SignalBackgroundMask,"all")/numel(SignalBackgroundMask)*100;
Signal2=Signal(not(SignalBackgroundMask));

if numel(Signal2)<50
    return
end

% compute statistics on intensity and noise
SignalNoiseVariance2=SignalNoiseVariance(not(SignalBackgroundMask));
SliceSignalMean=mean(Signal2,"all");
SliceSignalStd=std(Signal2,0,"all");
SliceNoiseSDMean=mean(sqrt(SignalNoiseVariance2),"all");
SliceSNR=SliceSignalMean/SliceNoiseSDMean;

mask=Signal2>=prctile(Signal2,50,"all");
SliceSignalMean50=mean(Signal2(mask));
SliceSignalNoiseSD50=mean(sqrt(SignalNoiseVariance2(mask)));
SliceSNR50=SliceSignalMean50/SliceSignalNoiseSD50;


% TParentNoise=[];
% if not(isempty(ParentNoise))
%     if not(isempty(ParentNoise.DarkVariance))
%         DarkCurrent=ParentNoise.DarkCurrent;
%         DarkVariance=ParentNoise.DarkVariance;
%         ShotCoefficient=ParentNoise.ShotCoefficient;
%         MultiplicativeCoefficient=ParentNoise.MultiplicativeCoefficient;
%         ParentPercentage=[];
%         ParentPercentage50=[];
%         if isempty(ParentSignal)
%             WhiteNoiseSD=sqrt(DarkVariance);
%             ShotNoiseSD=sqrt(ShotCoefficient*SliceSignalMean);
%             MultNoiseSD=sqrt(MultiplicativeCoefficient*(SliceSignalMean^2));
%             NoiseSD=sqrt(WhiteNoiseSD^2+ShotNoiseSD^2+MultNoiseSD^2);
%             WhiteNoiseSD50=sqrt(DarkVariance);
%             ShotNoiseSD50=sqrt(ShotCoefficient*SliceSignalMean50);
%             MultNoiseSD50=sqrt(MultiplicativeCoefficient*(SliceSignalMean50^2));
%             NoiseSD50=sqrt(WhiteNoiseSD50^2+ShotNoiseSD50^2+MultNoiseSD50^2);
%         else
%             for i=1:size(ParentSignal,4)
%                 var=ParentSignal(:,:,:,i);
%                 var=var(not(SignalBackgroundMask));
%                 var2=mean(var,"all");
%                 WhiteNoiseSD(i)=sqrt(DarkVariance(i));
%                 ShotNoiseSD(i)=sqrt(ShotCoefficient(i)*var2);
%                 MultNoiseSD(i)=sqrt(MultiplicativeCoefficient(i)*(var2^2));
%                 NoiseSD(i)=sqrt(WhiteNoiseSD(i)^2+ShotNoiseSD(i)^2+MultNoiseSD(i)^2);
%                 mask=var>=prctile(var,50,"all");
%                 var=var(mask);
%                 var2=mean(var,"all");
%                 WhiteNoiseSD50(i)=sqrt(DarkVariance(i));
%                 ShotNoiseSD50(i)=sqrt(ShotCoefficient(i)*var2);
%                 MultNoiseSD50(i)=sqrt(MultiplicativeCoefficient(i)*(var2^2));
%                 NoiseSD50(i)=sqrt(WhiteNoiseSD50(i)^2+ShotNoiseSD50(i)^2+MultNoiseSD50(i)^2);
%             end
%             NoiseSDtot=sqrt(sum(NoiseSD.^2,"all"));
%             NoiseSD50tot=sqrt(sum(NoiseSD50.^2,"all"));
%             for i=1:numel(NoiseSD)
%                 ParentPercentage(i)=(NoiseSD(i)/NoiseSDtot)^2*100;
%                 ParentPercentage50(i)=(NoiseSD50(i)/NoiseSD50tot)^2*100;
%             end
%         end
%         WhiteNoisePercentage=(WhiteNoiseSD./NoiseSD).^2*100;
%         ShotNoisePercentage=(ShotNoiseSD./NoiseSD).^2*100;
%         MultNoisePercentage=(MultNoiseSD./NoiseSD).^2*100;
%         WhiteNoisePercentage50=(WhiteNoiseSD50./NoiseSD50).^2*100;
%         ShotNoisePercentage50=(ShotNoiseSD50./NoiseSD50).^2*100;
%         MultNoisePercentage50=(MultNoiseSD50./NoiseSD50).^2*100;
%         VariableNames={'Slice Number'};
%         TParentNoise=table({SliceNumber},'VariableNames',VariableNames);
%         if not(isempty(ParentPercentage))
%             for i=1:numel(ParentPercentage)
%                 VariableNames={['Ch' num2str(i) ' Parent Noise Percentage'], ['Ch' num2str(i) ' Parent Noise Percentage Full']};
%                 T2=table({ParentPercentage(i)},{ParentPercentage50(i)},'VariableNames',VariableNames);
%                 TParentNoise=cat(2,TParentNoise,T2);
%             end
%         end
%         for i=1:numel(WhiteNoisePercentage)
%              VariableNames={['White Noise Percentage in Ch' num2str(i)], ['Shot Noise Percentage in Ch' num2str(i)], ['Multiplicative Noise Percentage in Ch' num2str(i)], ['White Noise Percentage in Full Ch' num2str(i)], ['Shot Noise Percentage in Full Ch' num2str(i)], ['Multiplicative Noise Percentage in Full Ch' num2str(i)]};
%              T2=table({WhiteNoisePercentage(i)},{ShotNoisePercentage(i)},{MultNoisePercentage(i)},{WhiteNoisePercentage50(i)},{ShotNoisePercentage50(i)},{MultNoisePercentage50(i)},'VariableNames',VariableNames);
%              TParentNoise=cat(2,TParentNoise,T2);
%         end
%     end
% end

% write table with statistics about slice intensity and noise
VariableNames={'Slice Number','Background %','Foreground Signal Mean','Foreground Noise SD','Foreground SNR','Full Signal Area %','Full Signal Mean','Full Signal Noise SD','Full Signal SNR'};
TNoise=table({SliceNumber},{SliceBackgroundPercentage},{SliceSignalMean},{SliceNoiseSDMean},{SliceSNR},{(100-SliceBackgroundPercentage)/2},{SliceSignalMean50},{SliceSignalNoiseSD50},{SliceSNR50},'VariableNames',VariableNames);

% gradient analysis
GradientBackgroundMask=GradientClass==1;
SliceNoGradientPercentage=sum(GradientBackgroundMask,"all")/numel(mask)*100;
mask=GradientClass==2;
SliceGradientPercentage=sum(mask,"all")/numel(mask)*100;

G1=G1(mask); % gradient analysis is performed only in foreground and significant gradient SNR
G2=G2(mask);
Graw1=Graw1(mask);
Graw2=Graw2(mask);
G1var=G1var(mask);
G2var=G2var(mask);


% note that g1 and g2 are swapped and changed sign to have 0 angle along x axis of image
[theta_mean,theta_std,SE,SEspread,spread_mean,theta_exkurt,GradientHisto]=circularHistogramStatistics(G1,G2,[],false);
VariableNames={'Slice Number','Angle Mean','Angle Std','Von Mises Angle Mean','Von Mises Angle Std','Angle at Maximum','Half Maximum Dispersion','Shannon Entropy','Entropic Dispersion','Half Mean Dispersion','Angle Excess Kurtosis'};
TGradientHisto=table({SliceNumber},{theta_mean},{theta_std},{GradientHisto.VMmu},{GradientHisto.VMsigma},{GradientHisto.MaximumAngle},{GradientHisto.HalfMaximumDispersion},{SE},{SEspread},{spread_mean},{theta_exkurt},'VariableNames',VariableNames);

% comparison of gradient direction with generative field, if present
if not(isempty(GenerativeField))
    gf1=GenerativeField(:,:,:,1);
    gf2=GenerativeField(:,:,:,2);
    gf1=gf1(mask);
    gf2=gf2(mask);
    gfn=sqrt(gf1.^2+gf2.^2)+eps;
    gf1=gf1./gfn;
    gf2=gf2./gfn;
    gn=sqrt(G1.^2+G2.^2)+eps;
    g1b=G1./gn;
    g2b=G2./gn;
    var=min(max(-g1b.*gf2+g2b.*gf1,-1),1);
    delta_ang=acos(var);

    delta_ang=delta_ang/pi*180;
    [me,sd]=circular_stat_ang(delta_ang,180,gn);
    GF_Gradient_RMS=sqrt(me^2+sd^2);
end

% global structure analysis on all the slices
if GlobalStructureBiasRemoval
    M11=sum(Graw1.^2-G1var,"all");
    M22=sum(Graw2.^2-G2var,"all");
    % Mbias=sum(gvar,"all");
    % M11=sum(Graw1.^2,"all")-Mbias/2;
    % M22=sum(Graw2.^2,"all")-Mbias/2;
    M12=sum(Graw1.*Graw2,"all");
    [Lmax,Lmin,Vmax1,Vmax2]=matrixDiagonalize2D(M11,M22,M12);
    Ltot=max(Lmin+Lmax,0);
    Lmin=max(Lmin,0);
    Lmax=Ltot-Lmin;
elseif not(GlobalStructureBiasRemoval) % if the global structure is computed without bias removal (false), the processed gradient is used (which can be bias corrected)
    M11=sum(G1.^2,"all");
    M22=sum(G2.^2,"all");
    M12=sum(G1.*G2,"all");
    [Lmax,Lmin,Vmax1,Vmax2]=matrixDiagonalize2D(M11,M22,M12);
else % if GlobalStructureBiasRemoval is neither true nor false (for instance empty), use raw gradient
    M11=sum(Graw1.^2,"all");
    M22=sum(Graw2.^2,"all");
    M12=sum(Graw1.*Graw2,"all");
    [Lmax,Lmin,Vmax1,Vmax2]=matrixDiagonalize2D(M11,M22,M12);
end

[SliceGlobalCoherence,SliceGlobalDA,SliceGlobalEccentricity,SliceGlobalAngle] = structureDescriptors(Lmax,Lmin,Vmax1,Vmax2,{'Coherence','DA','Eccentricity','Second Angle'});

VariableNames={'Slice Number','Background %','Foreground Signal Mean','Foreground Signal Std','No Gradient %','Gradient %','Coherence','Degree of Anisotropy','Eccentricity','Global Angle'};
TGradient=table({SliceNumber},{SliceBackgroundPercentage},{SliceSignalMean},{SliceSignalStd},{SliceNoGradientPercentage},{SliceGradientPercentage},{SliceGlobalCoherence},{SliceGlobalDA},{SliceGlobalEccentricity},{SliceGlobalAngle},'VariableNames',VariableNames);
         
[TotalVariation,Coherence,DA,Eccentricity,SecondAngle,SecondEigenvector] = structureDescriptors(FirstVariation,SecondVariation,FirstEigenvector1,FirstEigenvector2,{'Total Variation','Coherence','DA','Eccentricity','Second Angle','Second Eigenvector'});

mask=StructureClass==1;
SliceNoStructurePercentage=sum(mask,"all")/numel(mask)*100;
MaskStructure=StructureClass>=2;
SliceStructurePercentage=sum(MaskStructure,"all")/numel(mask)*100;

% compute weights for statistics
if strcmp(StructureWeightType,'Variation Difference')
   Weights=(FirstVariation-SecondVariation);
elseif strcmp(StructureWeightType,'Total Variation')
    Weights=TotalVariation;
elseif strcmp(StructureWeightType,'Signal')
    Weights=Signal;
elseif strcmp(StructureWeightType,'First Variation')
    Weights=FirstVariation;
else
    Weights=[];
end

if not(isempty(Weights))
    SecondEigenvector=SecondEigenvector.*Weights;
end
e1=SecondEigenvector(:,:,:,1);
e2=SecondEigenvector(:,:,:,2);
e1=e1(MaskStructure);
e2=e2(MaskStructure);
mask=StructureClass==2;
SliceAnisotropyPercentage=sum(mask,"all")/numel(mask)*100;
mask=StructureClass==3;
SliceIsotropyPercentage=sum(mask,"all")/numel(mask)*100;

% compute circular statistics on the structure derived direction
% note that e1 and e2 are swapped and changed sign to have 0 angle along x axis of image
[theta_mean,theta_std,SE,SEspread,spread_mean,theta_exkurt,StructureHisto]=circularHistogramStatistics(e2,-e1,[],false);

VariableNames={'Slice Number','Angle Mean','Angle Std','Von Mises Angle Mean','Von Mises Angle Std','Angle at Maximum','Half Maximum Dispersion','Shannon Entropy','Entropic Dispersion','Half Mean Dispersion','Angle Excess Kurtosis'};
TStructureHisto=table({SliceNumber},{theta_mean},{theta_std},{StructureHisto.VMmu},{StructureHisto.VMsigma},{StructureHisto.MaximumAngle},{StructureHisto.HalfMaximumDispersion},{SE},{SEspread},{spread_mean},{theta_exkurt},'VariableNames',VariableNames);

% comparison of structure direction with generative field 
if not(isempty(GenerativeField))
    gf1=GenerativeField(:,:,:,1);
    gf2=GenerativeField(:,:,:,2);
    gf1=gf1(MaskStructure);
    gf2=gf2(MaskStructure);
    gfn=sqrt(gf1.^2+gf2.^2)+eps;
    gf1=gf1./gfn;
    gf2=gf2./gfn;
    en=sqrt(e1.^2+e2.^2)+eps;
    e1b=e1./en;
    e2b=e2./en;
    var=max(min(e2b.*gf2+e1b.*gf1,1),-1);
    delta_ang=acos(var);

    delta_ang=delta_ang/pi*180;
    [me,sd]=circular_stat_ang(delta_ang,180,en);
    GF_Structure_RMS=sqrt(me^2+sd^2);

    VariableNames={'Slice Number','Gradient Angle RMS Error','Structure Angle RMS Error'};
    TGenerativeField=table({SliceNumber},{GF_Gradient_RMS},{GF_Structure_RMS},'VariableNames',VariableNames);
end


if isempty(Weights)
    SliceCoherenceMean=mean(Coherence(MaskStructure),"all");
    SliceCoherenceStd=std(Coherence(MaskStructure),0,"all");
    SliceDAMean=mean(DA(MaskStructure),"all");
    SliceDAStd=std(DA(MaskStructure),0,"all");
    SliceEccentricityMean=mean(Eccentricity(MaskStructure),"all");
    SliceEccentricityStd=std(Eccentricity(MaskStructure),0,"all");
else
    Weights=Weights(MaskStructure);
    var=sum(Weights,"all");
    SliceCoherenceMean=sum(Coherence(MaskStructure).*Weights,"all")/var;
    SliceCoherenceStd=sqrt(sum(((Coherence(MaskStructure)-SliceCoherenceMean).^2).*Weights,"all")/var);
    SliceDAMean=sum(DA(MaskStructure).*Weights,"all")/var;
    SliceDAStd=sqrt(sum(((DA(MaskStructure)-SliceCoherenceMean).^2).*Weights,"all")/var);
    SliceEccentricityMean=sum(Eccentricity(MaskStructure).*Weights,"all")/var;
    SliceEccentricityStd=sqrt(sum(((Eccentricity(MaskStructure)-SliceCoherenceMean).^2).*Weights,"all")/var);
end

[me,sd]=circular_stat_ang(SecondAngle(MaskStructure),180); % non weighted angle
SliceAngleMean=me;
SliceAngleStd=sd;
var1=mean(FirstVariation(MaskStructure),"all");
SliceDAIndirectMean=(var1-mean(SecondVariation(MaskStructure),"all"))/(var1+eps);

VariableNames={'Slice Number','Background %','No Strucure %','Structure %','Anisotropy %','Isotropy %','Coherence Mean','Coherence Std','Degree of Anisotropy Mean','Degree of Anisotropy Std','Degree of Anisotropy Indirect Mean','Eccentricity Mean','Eccentricity Std','Angle Mean','Angle Std'};
TStructure=table({SliceNumber},{SliceBackgroundPercentage},{SliceNoStructurePercentage},{SliceStructurePercentage},{SliceAnisotropyPercentage},{SliceIsotropyPercentage},{SliceCoherenceMean},{SliceCoherenceStd},{SliceDAMean},{SliceDAStd},{SliceDAIndirectMean},{SliceEccentricityMean},{SliceEccentricityStd},{SliceAngleMean},{SliceAngleStd},'VariableNames',VariableNames);
      

%-------------------------------------------------------------------------------
function oneSliceGraphicalResults(TGradient, TGradientHisto, GradientHisto, TStructure, TStructureHisto, StructureHisto,STParameters,SliceNumber)
% create figure that show results in the exemplative slice

flag_full=true;

sangle=uistyle("BackgroundColor",[1 0.8 0.8]);
sangle2=uistyle("FontWeight","bold");
saniso=uistyle("BackgroundColor",[1 1 0.8]);
sperc=uistyle("BackgroundColor",[0.8 0.8 1]);
ssig=uistyle("BackgroundColor",[0.8 1 0.8]);

hf=uifigure('Name',['Slice ' num2str(SliceNumber) ' Statistics'],'Position',[200 200 1200 420]);
if STParameters.GradientNormalization
    var="Normalized Gradient Angular Distribution Analysis";
else
    var="Gradient Angular Distribution Analysis";
end
uilabel(hf,"Text",var,"Position",[20 390 1200-40 20]);
ht=uitable(hf,'Data',TGradientHisto,'Position',[20 320 1200-40 70]);
addStyle(ht,sangle,'cell',[1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8; 1 9; 1 10]);
addStyle(ht,sangle2,'cell',[1 2; 1 4; 1 6]);
ht.Tooltip='Non-directional angles of the vectors normal to the gradient (along the fibers) in the foreground area. Different average and dispersion mesures are computed by circular statistics.';

var=['Structure Angular Distribution Analysis (weighted by ' char(STParameters.StructureWeightType) ')'];
uilabel(hf,"Text",var,"Position",[20 290 1200-40 20]);
ht=uitable(hf,'Data',TStructureHisto,'Position',[20 220 1200-40 70]);
addStyle(ht,sangle,'cell',[1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8; 1 9; 1 10]);
addStyle(ht,sangle2,'cell',[1 2; 1 4; 1 6]);
ht.Tooltip='Non-directional angles of the second eigenvector vectors (along the fibers) in the foreground area. Different average and dispersion mesures are computed by circular statistics.';

uilabel(hf,"Text","Signal and Gobal Structure Descriptors","Position",[20 190 1200-40 20]);
ht=uitable(hf,'Data',TGradient,'Position',[20 120 1200-40 70]);
addStyle(ht,sangle,'cell',[1 10]);
addStyle(ht,sangle2,'cell',[1 10]);
addStyle(ht,saniso,'cell',[1 7; 1 8; 1 9]);
addStyle(ht,sperc,'cell',[1 2; 1 5; 1 6]);
addStyle(ht,ssig,'cell',[1 3; 1 4]);
ht.Tooltip='Background percentage area and signal mean  and dispersion in the foreground area. Descriptor of the general structure tensor obtained by diagonalizing the 2D covariance matrix of the gradient components averaged all over the foreground area.';

uilabel(hf,"Text","Average of Local Structure Descriptors ","Position",[20 90 1200-40 20]);
ht=uitable(hf,'Data',TStructure,'Position',[20 20 1200-40 70]);
addStyle(ht,sangle,'cell',[1 14; 1 15]);
addStyle(ht,sangle2,'cell',[1 14]);
addStyle(ht,saniso,'cell',[1 7; 1 8; 1 9; 1 10; 1 11; 1 12; 1 13]);
addStyle(ht,sperc,'cell',[1 2; 1 3; 1 4; 1 5; 1 6]);
ht.Tooltip='Area percentages of the classification labels. Average and dispersion of local structure descriptors. Local structure obtained from the local average of the covariance matrix of the gradient components. The structure scale defines the size of the local averanging.';


hf=figure('Name',['Slice ' num2str(SliceNumber) ' Angular Distributions']);
tiledlayout(hf,"vertical")

nexttile
Bins=GradientHisto.Bins;
AngularDistribution=GradientHisto.AngularDistribution;
AngularDistributionSE=GradientHisto.AngularDistributionSE;
ang_off2=0;
if TGradientHisto.('Angle Mean'){1}>0
    mask=Bins<(TGradientHisto.('Angle Mean'){1}-90);
    Bins(mask)=Bins(mask)+180;
    if GradientHisto.MaximumAngle<(TGradientHisto.('Angle Mean'){1}-90)
        ang_off2=180;
    end
else
    mask=Bins>(TGradientHisto.('Angle Mean'){1}+90);
    Bins(mask)=Bins(mask)-180;
    if GradientHisto.MaximumAngle>(TGradientHisto.('Angle Mean'){1}+90)
        ang_off2=-180;
    end
end
had=bar(Bins,AngularDistribution);
set(had,'FaceColor',[0.5 0 0],'EdgeColor',[0.5 0 0])
hold on
if flag_full
    had=bar(Bins,AngularDistributionSE);
    set(had,'FaceColor',[0.3 0 0],'EdgeColor',[0.3 0 0])
    if numel(GradientHisto.VMc)==2
        y=GradientHisto.VMc(1)+GradientHisto.VMc(2)*exp(cos((Bins-GradientHisto.VMmu)/180*2*pi)/(GradientHisto.VMsigma/180*2*pi)^2);
    else
        y=vonMisesDistribution(Bins/180*2*pi,GradientHisto.VMmu/180*2*pi,GradientHisto.VMsigma/180*2*pi);
    end
    plot(Bins,y,'Color','k','LineStyle','none','Marker','*','LineWidth',2);
    plot(Bins,GradientHisto.AngularDistributionSmooth,'Color',[0.5 0 0],'LineStyle','none','Marker','o','LineWidth',2);
    plot(xlim,[GradientHisto.HalfMeanThreshold GradientHisto.HalfMeanThreshold],'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',1.5)
    plot(xlim,[GradientHisto.Maximum/2 GradientHisto.Maximum/2],'Color','m','LineStyle','--','LineWidth',0.5)
    plot([TGradientHisto.('Angle Mean'){1} TGradientHisto.('Angle Mean'){1}],ylim,'r','LineWidth',2)
    plot([GradientHisto.MaximumAngle GradientHisto.MaximumAngle]+ang_off2,ylim,'r','LineWidth',1)
    plot([TGradientHisto.('Angle Mean'){1} TGradientHisto.('Angle Mean'){1}]-TGradientHisto.('Angle Std'){1},ylim,'r','LineStyle',':','LineWidth',1.5)
    plot([TGradientHisto.('Angle Mean'){1} TGradientHisto.('Angle Mean'){1}]+TGradientHisto.('Angle Std'){1},ylim,'r','LineStyle',':','LineWidth',1.5)
    legend('Angular distr.','Ang. distr. under entropic dispersion','von Mises Interpolation','Smooth Angular distr.','Half mean threshold', 'Half Maximum Threshold')
    text(0.01,0.95,['Percentage of excluded null vectors : ' num2str(GradientHisto.ExclusionPercentage)],'Units','normalized');
    text(0.01,0.90,['Angle mean: ' num2str(TGradientHisto.('Angle Mean'){1}) ' +- ' num2str(TGradientHisto.('Angle Std'){1})],'Units','normalized','Color','r');
    text(0.01,0.85,['Von Mises Angle mean: ' num2str(GradientHisto.VMmu) ' +- ' num2str(GradientHisto.VMsigma)],'Units','normalized','Color','r');
    text(0.01,0.80,['Angle at maximum +- half maximum dispersion/2: ' num2str(GradientHisto.MaximumAngle) ' +- ' num2str(GradientHisto.HalfMaximumDispersion/2)],'Units','normalized','Color','r');
    text(0.01,0.75,['Shannon entropy: ' num2str(TGradientHisto.('Shannon Entropy'){1})],'Units','normalized');
    text(0.01,0.71,['Entropic dispersion: ' num2str(TGradientHisto.('Entropic Dispersion'){1})],'Units','normalized','Color',[0.6 0 0]);
    text(0.01,0.65,['Half mean dispersion: ' num2str(TGradientHisto.('Half Mean Dispersion'){1})],'Units','normalized','Color',[0.4 0.4 0.4]);
    text(0.01,0.60,['Excess Kurtosis (Mardia): ' num2str(TGradientHisto.('Angle Excess Kurtosis'){1})],'Units','normalized','Color','r');
else
    plot([TGradientHisto.('Angle Mean'){1} TGradientHisto.('Angle Mean'){1}],ylim,'r','LineWidth',2)
    plot([TGradientHisto.('Angle Mean'){1} TGradientHisto.('Angle Mean'){1}]-TGradientHisto.('Angle Std'){1},ylim,'r','LineStyle',':','LineWidth',1.5)
    plot([TGradientHisto.('Angle Mean'){1} TGradientHisto.('Angle Mean'){1}]+TGradientHisto.('Angle Std'){1},ylim,'r','LineStyle',':','LineWidth',1.5)
    % legend('Angular distr.','Angle Mean','Angle SD Confidence Interval')
    text(0.01,0.90,['Angle Mean +- SD: ' num2str(TGradientHisto.('Angle Mean'){1}) ' +- ' num2str(TGradientHisto.('Angle Std'){1})],'Units','normalized','Color','r');
end
ylimits=ylim;
ax=gca;
if STParameters.GradientNormalization
    title("Normalized Gradient Angular Distribution");
else
    title("Gradient Angular Distribution");
end

nexttile
Bins=StructureHisto.Bins;
AngularDistribution=StructureHisto.AngularDistribution;
AngularDistributionSE=StructureHisto.AngularDistributionSE;
ang_off=0;
ang_off2=0;
if TGradientHisto.('Angle Mean'){1}>0
    mask=Bins<(TGradientHisto.('Angle Mean'){1}-90);
    Bins(mask)=Bins(mask)+180;
    if TStructureHisto.('Angle Mean'){1}<(TGradientHisto.('Angle Mean'){1}-90)
        ang_off=180;
    end
    if StructureHisto.MaximumAngle<(TGradientHisto.('Angle Mean'){1}-90)
        ang_off2=180;
    end
else
    mask=Bins>(TGradientHisto.('Angle Mean'){1}+90);
    Bins(mask)=Bins(mask)-180;
    if TStructureHisto.('Angle Mean'){1}>(TGradientHisto.('Angle Mean'){1}+90)
        ang_off=-180;
    end
    if StructureHisto.MaximumAngle>(TGradientHisto.('Angle Mean'){1}+90)
        ang_off2=-180;
    end
end
had=bar(Bins,AngularDistribution);
set(had,'FaceColor',[0.5 0 0.5],'EdgeColor',[0.5 0 0.5])
hold on
if flag_full
    had=bar(Bins,AngularDistributionSE);
    set(had,'FaceColor',[0.3 0 0.3],'EdgeColor',[0.3 0.4 0.3])
    if numel(StructureHisto.VMc)==2
        y=StructureHisto.VMc(1)+StructureHisto.VMc(2)*exp(cos((Bins-StructureHisto.VMmu)/180*2*pi)/(StructureHisto.VMsigma/180*2*pi)^2);
    else
        % y=StructureHisto.VMc*exp(cos((Bins-StructureHisto.VMmu)/180*2*pi)/(StructureHisto.VMsigma/180*2*pi)^2);
        y=vonMisesDistribution(Bins/180*2*pi,StructureHisto.VMmu/180*2*pi,StructureHisto.VMsigma/180*2*pi);
    end
    plot(Bins,y,'Color','k','LineStyle','none','Marker','*','LineWidth',2);
    plot(Bins,StructureHisto.AngularDistributionSmooth,'Color',[0.5 0 0.5],'LineStyle','none','Marker','o','LineWidth',2);
    plot(xlim,[StructureHisto.HalfMeanThreshold StructureHisto.HalfMeanThreshold],'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',1.5)
    plot(xlim,[StructureHisto.Maximum/2 StructureHisto.Maximum/2],'Color','m','LineStyle','--','LineWidth',0.5)
    plot([TStructureHisto.('Angle Mean'){1} TStructureHisto.('Angle Mean'){1}]+ang_off,ylim,'m','LineWidth',2)
    plot([StructureHisto.MaximumAngle StructureHisto.MaximumAngle]+ang_off2,ylim,'m','LineWidth',1)
    plot([TStructureHisto.('Angle Mean'){1} TStructureHisto.('Angle Mean'){1}]-TStructureHisto.('Angle Std'){1}+ang_off,ylim,'m','LineStyle',':','LineWidth',1.5)
    plot([TStructureHisto.('Angle Mean'){1} TStructureHisto.('Angle Mean'){1}]+TStructureHisto.('Angle Std'){1}+ang_off,ylim,'m','LineStyle',':','LineWidth',1.5)
    legend('Angular distr.','Ang. distr. under entropic dispersion','von Mises Interpolation','Smooth Angular distr.', 'Half mean threshold', 'Half Maximum Threshold')
    text(0.01,0.95,['Percentage of excluded null vectors : ' num2str(StructureHisto.ExclusionPercentage)],'Units','normalized');
    text(0.01,0.90,['Angle mean: ' num2str(TStructureHisto.('Angle Mean'){1}) ' +- ' num2str(TStructureHisto.('Angle Std'){1})],'Units','normalized','Color',[1 0 1]);
    text(0.01,0.85,['Von Mises Angle mean: ' num2str(StructureHisto.VMmu) ' +- ' num2str(StructureHisto.VMsigma)],'Units','normalized','Color',[1 0 1]);
    text(0.01,0.80,['Angle at maximum +- half maximum dispersion/2: ' num2str(StructureHisto.MaximumAngle) ' +- ' num2str(StructureHisto.HalfMaximumDispersion/2)],'Units','normalized','Color','m');
    text(0.01,0.75,['Shannon entropy: ' num2str(TStructureHisto.('Shannon Entropy'){1})],'Units','normalized');
    text(0.01,0.71,['Entropic dispersion: ' num2str(TStructureHisto.('Entropic Dispersion'){1})],'Units','normalized','Color',[0.6 0 0.6]);
    text(0.01,0.65,['Half mean dispersion: ' num2str(TStructureHisto.('Half Mean Dispersion'){1})],'Units','normalized','Color',[0.4 0.4 0.4]);
    text(0.01,0.60,['Excess Kurtosis (Mardia): ' num2str(TStructureHisto.('Angle Excess Kurtosis'){1})],'Units','normalized','Color','m');
else
    plot([TStructureHisto.('Angle Mean'){1} TStructureHisto.('Angle Mean'){1}]+ang_off,ylim,'m','LineWidth',2)
    plot([TStructureHisto.('Angle Mean'){1} TStructureHisto.('Angle Mean'){1}]-TStructureHisto.('Angle Std'){1}+ang_off,ylim,'m','LineStyle',':','LineWidth',1.5)
    plot([TStructureHisto.('Angle Mean'){1} TStructureHisto.('Angle Mean'){1}]+TStructureHisto.('Angle Std'){1}+ang_off,ylim,'m','LineStyle',':','LineWidth',1.5)
    text(0.01,0.90,['Angle Mean +- SD: ' num2str(TStructureHisto.('Angle Mean'){1}) ' +- ' num2str(TStructureHisto.('Angle Std'){1})],'Units','normalized','Color',[1 0 1]);
end
title(['Structure Angular Distribution (weighted by ' char(STParameters.StructureWeightType) ')']);

ylim(ax,ylimits);
