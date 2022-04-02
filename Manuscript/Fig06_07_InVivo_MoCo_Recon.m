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
    %%% Load in-vivo parameters
    InViVoParams;%comment

    %%% Load reconstructions and register
    Xall = [];
    for f=1:length(experimentNames)
        
        fileName=fullfile(mainDir, 'Results','Exp2_Motion',sprintf('HV%d',HV),'An-Ve',experimentNames{f});
        suff={'Aq_MotCorr','Di_MotCorr', 'Di_MotB0Corr'};%Uncorrected - Motion corrected - Motion+B0 corrected
        [xTemp,MSTemp,MTTemp] = readNII(fileName,suff);

        %%% Align
        x = []; x = {xTemp{1}, xTemp{2}, xTemp{3}, xGT }; %Uncorrected - Motion corrected - Motion+B0 corrected - Ground Truth
        res = {MSTemp{1}, MSTemp{2}, MSTemp{3},MSxGT};%Resolution
        orient = {MTTemp{1}, MTTemp{2}, MTTemp{3},MTxGT};%Orientation
        ref= 3;%Take dephCorr as a reference
        exclDim=[3 , .3];
        padding = [0 0 0];
        if gpu; x = gatherStruct(x,0);end%Load into the gpu to accelerate regisrtation
        xAligned = alignVolumes(x, res, orient, ref, exclDim, 0, padding);

        if isempty(Xall); Xall = xAligned; else Xall = cat(4, Xall, xAligned);end
        %4th dim:  motion experiments/ 5th: recons and xGT
    end
    
    %%% Create mask
    Mmagn = abs(dynInd(Xall,{1,4},4:5))>th*rrX(2);%GT image
      
    N = size(Xall);
    Text = [];
    
    %%% Extract FOV and do NIFTI index handling
    for i=1:3; if isempty(supFOV{i}); supFOV{i} = 1:N(i);end;end
    if isempty(idxNIFTI); idxNIFTI = ceil((N(1:3)+1)/2);end
    idxPlotND = idxNIFTI - ([supFOV{1}(1) supFOV{2}(1) supFOV{3}(1)]-1);
    Xall = dynInd(Xall,supFOV,1:3);
    Mmagn = dynInd(Mmagn,supFOV,1:3);

    %% Figure 6
    N = size(Xall); %4th dimension are the different recons/ 5th the different experiment
    xPlot = Xall;%4th dimension are the different recons/ 5th the different experiment
    plotND([],abs(xPlot),rrX,idxPlotND,0,{2,1},MTxGT,Text, Mmagn, {0}, 101);
    saveFig(fullfile(saveDir,sprintf('Figure6_HV%d',HV)),[],saveFlag,[],1,imageResol,[],imageType);

    %% Figure 7
    xPlot = dynInd(Xall,3,4);%4th dimension are the different recons/ 5th the different experiment

    plotND([],abs(xPlot),rrX,idxPlotND,0,{[],1},MTxGT,Text, Mmagn, {0}, 101);
    %title(sprintf('Linear maps estimated from different motion-free scan combinations.'),'Interpreter','latex','FontSize',16,'Linewidth',15)
    saveFig(fullfile(saveDir,sprintf('Figure7_HV%d',HV)),[],saveFlag,[],1,imageResol,[],imageType);
    
    %% Table 2 - similarity metrics
    Xall = permute(Xall,[1:3 5 4]); %4th dimension are the different recons/ 5th the different experiment
    Xall  = gather(Xall);

    SNR = zeros( (size(Xall,4)-1), size(Xall,5));
    SSIM = zeros( (size(Xall,4)-1), size(Xall,5));
    MI = zeros( (size(Xall,4)-1), size(Xall,5));
    AP = zeros( (size(Xall,4)-1), size(Xall,5));

    for i =1:size(Xall,4)-1%Reconstructions - don't include the GT
        for j = 1:size(Xall, 5) %Experiments

            im = abs( dynInd(Xall,[i,j],4:5));%Use magnitude since phase is not guaranteed to be the same
            ref = abs(dynInd(Xall , [size(Xall,4), 1],4:5));

            SNR(i,j) = gather ( getSNR(im,ref) );
            SSIM(i,j) = gather ( getSSI(im,ref,[],[],0) );
            MI(i,j) = gather ( getMI(im,ref) );
            AP(i,j) = gather(getAP(im,ref));

        end  
    end

    %Compute where metrics have improved
    SSIM_Improved = inf*onesL(SSIM); SSIM_Improved(2:end,:) = (SSIM(2:end,:) - SSIM(1:end-1,:))>0;
    MI_Improved = inf*onesL(MI); MI_Improved(2:end,:) = (MI(2:end,:) - MI(1:end-1,:))>0;
    SNR_Improved = inf*onesL(SNR); SNR_Improved(2:end,:) = (SNR(2:end,:) - SNR(1:end-1,:))>0;
    AP_Improved = inf*onesL(AP); AP_Improved(2:end,:) = (AP(2:end,:) - AP(1:end-1,:))<0;
    
    save(fullfile(saveDir, sprintf('Table2_Metrics_HV%d.mat',HV)),'SNR','SSIM','MI','AP','SNR_Improved','SSIM_Improved','MI_Improved','AP_Improved');

end

