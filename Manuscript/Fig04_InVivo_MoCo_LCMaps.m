
clc; clear all; close all;
currentDir = mfilename('fullpath');
cd(fileparts(currentDir));

%% SET PATHS
addpath(genpath('..'));
pathInfo; 
figureSpecs;

experimentNames = {'Low','Medium','High'};

%% ITERATE OVER HEALTHY VOLUNTEERS
for HV = [1 2]
    
    %% LOAD PARAMETERS FOR IN-VIVO
    InViVoParams;
    
    %% LOAD AND REGISTER RECONSTRUCTIONS
    Xall = [];%comment
    for f=1:length(experimentNames)

        fileName = sprintf('%s/Results/Exp2_Motion/HV%d/An-Ve/%s',mainDir,HV, experimentNames{f});%comment
        suff={'Dc_1_MotB0Corr','Dc_2_MotB0Corr','Di_MotB0Corr'};
        [xTemp,MSTemp,MTTemp] = readNII(fileName,suff);

%         %%% Match xGT image
%         matchHistFlag=0;pctAlign = [5:5:95];
%         if HV==1 && ~exist('matchingCompleteXGT','var')%Only to be performed once
%             xGT = matchHist(xGT,xTemp{3},pctAlign, 0);
%             matchingCompleteXGT=1;
%         end

        %%% Align
        if HV==1
             x = {xTemp{3},xGT,xDGT; ...
                 cat(4,xTemp{1}, xTemp{2}), [], dynInd(DGT,[2 1],4)};

             res = {MSTemp{3},  MSxGT, MSDGT};
             orient ={MTTemp{3}, MTxGT, MTDGT};
        elseif HV==2
             x = {xTemp{3},xGT; ...
                 cat(4,xTemp{1}, xTemp{2}),[] };%No GT LC maps (see manuscript)

             res = {MSTemp{3},  MSxGT};
             orient = {MTTemp{3}, MTxGT};                  
        end

        ref = 1;%Take xTemp as a reference so that indexing is correct
        exclDim = [3, 0.3];
        padding = [0 0 10]; 
        if gpu; x = gatherStruct(x,0);end
        xAligned = alignVolumes(x, res, orient, ref, exclDim, 0, padding);

        if isempty(Xall); Xall = xAligned; else Xall = cat(6, Xall, xAligned);end
        %4th dim: x and LCmaps / 5th: recon, GT and dGT / 6th: motion experiments
    end

    %% CREATE MASK FOR PLOTTING
    Mmagn = abs(dynInd(Xall,{1,2,1},4:6))>th*rrX(2);
     
    %% PLOT LC MAPS
    if HV==1
        xPlot = cat(5, ...
           dynInd(Xall, {2,1,1},4:6),...%LCM_pitch - from recon - low motion
           dynInd(Xall, {2,1,2},4:6),...%LCM_pitch - from recon - medium motion
           dynInd(Xall, {2,1,3},4:6),...%LCM_pitch - from recon - high motion
           dynInd(Xall, {2,3,1},4:6),...%LCM_pitch - ground truth 
           dynInd(Xall, {3,1,1},4:6),...%LCM_roll - from recon - low motion
           dynInd(Xall, {3,1,2},4:6),...%LCM_roll - from recon - medium motion
           dynInd(Xall, {3,1,3},4:6),...%LCM_roll - from recon - high motion
           dynInd(Xall, {3,3,1},4:6));%LCM_roll - ground truth 
    elseif HV==2
        xPlot = cat(5, ...
           dynInd(Xall, {2,1,1},4:6),...
           dynInd(Xall, {2,1,2},4:6),...
           dynInd(Xall, {2,1,3},4:6),...
           dynInd(Xall, {3,1,1},4:6),...
           dynInd(Xall, {3,1,2},4:6),...
           dynInd(Xall, {3,1,3},4:6));
    end

    xPlot = pi/180*xPlot;%In Hz/degree
    N = size(Xall);
    
    %%% Extract FOV and do NIFTI index handling
    for i=1:3; if isempty(supFOV{i}); supFOV{i} = 1:N(i);end;end
    if isempty(idxNIFTI_LCM); idxNIFTI_LCM = ceil((N(1:3)+1)/2);end
    idxPlotND = idxNIFTI_LCM -  ([supFOV{1}(1) supFOV{2}(1) supFOV{3}(1)]-1);
    xPlot = dynInd(xPlot,supFOV,1:3);
    Mmagn = dynInd(Mmagn,supFOV,1:3);
    
    Text = [];    
    plotND([], real(xPlot), rrD, idxPlotND, 1, {1:3,1}, MTxGT, Text, Mmagn, {3}, 101,[],'[$Hz/ ^\circ$]');

    saveFig(fullfile(saveDir,sprintf('Figure4_LCM_InVivo_HV%d',HV)),[],saveFlag,[],1,imageResol,[],imageType);

end
