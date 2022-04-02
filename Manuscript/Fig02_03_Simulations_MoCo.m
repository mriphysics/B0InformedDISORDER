
clc; clear all; close all;
currentDir = mfilename('fullpath');
cd(fileparts(currentDir));

%% SET PATHS
addpath(genpath('..'));

pathInfo;
figureSpecs;

rrx = [0  150000]; 
xTh=rrx(2)*.2;

%% ITERATE OVER SNR LEVELS
for SNRdB = [30 inf] %Respectively 30dB and no noise.
    
    if SNRdB==30;suff = '_Noise';else;suff='';end
    
    %% SET PARAMETERS USED IN SIMULATIONS
    rec.Sim.resSynth = 1;
    rec.Sim.motionType = 'interleavedExponential';
    rec.Sim.NinterleavesT = 15;
    rec.Sim.tran = 10;
    rec.Sim.rot = 15;
    rec.Sim.snrdB = SNRdB; 

    rec.Sim.groupSweepsSynth = 2; 
    rec.Alg.parXT.groupSweeps = 2;
    rec.Sim.synthX = 1;     rec.Sim.synthT = 1;     rec.Sim.synthD = 1;   
    rec.Sim.resRecon = [1 ];

    %% LOAD SIMULATIONS RESULTS
    %Uncorrected
    rec.Sim.provideX = 0;   rec.Sim.provideT = 0;   rec.Sim.provideD = 0;
    rec.Sim.estX = 1;       rec.Sim.estT = 0;       rec.Sim.estD = 0;
    rec.Sim.nameSave = fullfile(mainDir,'Results','Simulations',simulationName(rec));
    ssUnCorr = load(rec.Sim.nameSave); 

    %Motion corrected
    rec.Sim.estX = 1;       rec.Sim.estT = 1;       rec.Sim.estD = 0;
    rec.Sim.nameSave = fullfile(mainDir,'Results','Simulations',simulationName(rec));
    ssMotCorr = load(rec.Sim.nameSave); 

    %Motion + B0 corrected
    rec.Sim.estX = 1;       rec.Sim.estT = 1;       rec.Sim.estD = 1;
    rec.Sim.nameSave = fullfile(mainDir,'Results','Simulations',simulationName(rec));
    ssMotB0Corr = load(rec.Sim.nameSave); 

    %(Motion + B0) provided
    rec.Sim.provideX = 0;   rec.Sim.provideT = 1;   rec.Sim.provideD = 1;
    rec.Sim.estX = 1;       rec.Sim.estT = 0;       rec.Sim.estD = 0;
    rec.Sim.nameSave = fullfile(mainDir,'Results','Simulations',simulationName(rec));
    ssMotB0Prov = load(rec.Sim.nameSave); 

    %Motion-free
    rec.Sim.synthX = 1;     rec.Sim.synthT = 0;     rec.Sim.synthD = 0;
    rec.Sim.provideX = 0;   rec.Sim.provideT = 0;   rec.Sim.provideD = 0;
    rec.Sim.estX = 1;       rec.Sim.estT = 0;       rec.Sim.estD = 0;
    rec.Sim.nameSave = fullfile(mainDir,'Results','Simulations',simulationName(rec));
    ssMF =  load(rec.Sim.nameSave); 

    %% LOAD MASK
    [M, ~, MTM ] = readNII(fullfile(mainDir,'Data','Exp1_Pose','HV1','7T','NIFTI','Pose1'),{'Mask'});
    %M = mapVolume(M{1},ssUnCorr.xGTRes, MTM{1},ssUnCorr.MT); 
    M = resampling(M{1}, size(ssUnCorr.xGTRes));
    M = single(M>.5);
    plotND([], abs(ssUnCorr.xGTRes), rrx,[],0,{[],2},ssUnCorr.MT,[], M, {2,[],3},100);

    %% LOAD IMAGES AND COMPUTE SNR
    xGT = ssUnCorr.xGTRes;
    xUnCorr = ssUnCorr.xRes;
    xMotB0Prov = ssMotB0Prov.xRes;
    xMotCorr = ssMotCorr.xRes; 
    xMotB0Corr = ssMotB0Corr.xRes; 
    xMF = ssMF.xRes; 

    %%% Compute SNR
    MSNR=onesL(M);%Disabling mask for metric calculation
    SNRUnCorr = getSNR(xUnCorr,xGT, MSNR);
    SNRMotCorr = getSNR(xMotCorr,xGT, MSNR);
    SNRMotB0Corr = getSNR(xMotB0Corr,xGT, MSNR);
    SNRMotB0Prov = getSNR(xMotB0Prov,xGT, MSNR);
    SNRMF = getSNR(xMF,xGT, MSNR);

    MSNR=M;
    SNRbrainUnCorr = getSNR(xUnCorr,xGT, MSNR);
    SNRbrainMotCorr = getSNR(xMotCorr,xGT, MSNR);
    SNRbrainMotB0Corr = getSNR(xMotB0Corr,xGT, MSNR);
    SNRbrainMotB0Prov = getSNR(xMotB0Prov,xGT, MSNR);
    SNRbrainMF = getSNR(xMF,xGT, MSNR);

    %% PLOT RECONSTRUCTIONS
    rrError = [0 .1* rrx(2)];
    X = abs( cat(5,xUnCorr,xMotCorr,xMotB0Corr,xMotB0Prov, xMF));
    X = cat(4,X, rescaleND( abs(xGT)-abs(X),rrx, rrError) );

    %%% Plot with contour
    plotND([], X, rrx,[],0,{1,1}, ssMF.MT, [], M, {2,[],3}, 100);
    saveFig(fullfile(saveDir,sprintf('Figure3_ImageRecon%s',suff)),[],saveFlag ,[],1,imageResol,[],imageType);

    %%% Plot without contour
    plotND([], X, rrx,[],0,{1,1}, ssMF.MT, [], M, {0}, 100);
    saveFig(fullfile(saveDir,sprintf('Figure3_ImageRecon_NoContour%s',suff)),[],saveFlag,[],1,imageResol,[],imageType);

    %%% Plot with annotated SNR
    Text = {sprintf('SNR %.2f dB', SNRUnCorr) , sprintf('SNR %.2f dB', SNRMotCorr), sprintf('SNR %.2f dB', SNRMotB0Corr), sprintf('SNR %.2f dB', SNRMotB0Prov),sprintf('SNR %.2f dB', SNRMF);...
            sprintf('SNRbrain %.2f dB', SNRbrainUnCorr) , sprintf('SNRbrain %.2f dB', SNRbrainMotCorr), sprintf('SNRbrain %.2f dB', SNRbrainMotB0Corr), sprintf('SNRbrain %.2f dB', SNRbrainMotB0Prov),sprintf('SNRbrain %.2f dB', SNRbrainMF)};
    plotND([], X, rrx,[],0,{1,1},ssMF.MT, Text, M, {2,[],3}, 100,[],[],[1 .5 1]);
    saveFig(fullfile(saveDir,sprintf('Figure3_ImageRecon_Annotated%s',suff)),[],saveFlag,[],1,imageResol,[],imageType);

    %%  LOAD LC MAPS
    DGT = pi/180*ssMotB0Corr.DGTRes;%In Hz/deg
    Mmagn = abs(xGT)>xTh;
    DMotB0Corr = pi/180*ssMotB0Corr.DRes;
    
    perm = [1:3 4 5:6];
    DGT = permute(DGT,perm);
    DMotB0Corr = permute(DMotB0Corr,perm);

    temp1 = cat(5, dynInd(DGT,2,6) , dynInd(DMotB0Corr,2,6) , dynInd(DGT,2,6)-dynInd(DMotB0Corr,2,6));
    temp2 = cat(5, dynInd(DGT,1,6) , dynInd(DMotB0Corr,1,6) , dynInd(DGT,1,6)-dynInd(DMotB0Corr,1,6));

    %%% Plot without annotation
    plotND([], cat(5, temp1,temp2),[-10 10],[],1,{[],1},ssMF.MT,[],Mmagn,{3,[0 0 0]},200,[],'[$Hz/ ^\circ$]');
    saveFig(fullfile(saveDir,sprintf('Figure2_LCMs%s',suff)),[],saveFlag,[],1,imageResol,[],imageType);
   

    %%% Plot with annotation
    plotND([], cat(5, temp1,temp2) , [-10 10], [],1, {[],1},ssMF.MT,...
            {'$d_{roll,GT}$','$d_{roll}$','$d_{roll,GT}-d_{roll}$','$d_{pitch,GT}$','$d_{pitch}$','$d_{pitch,GT}-d_{pitch}$'},...
             Mmagn, {3,[0 0 0]},200,[],'[$Hz/ ^\circ$]');
    saveFig(fullfile(saveDir,sprintf('Figure2_LCMs_Annotated%s',suff)),[],saveFlag,[],1,imageResol,[],imageType);


    %% PLOT MOTION TRACES
    TGT = ssMotCorr.TGT;
    TMotCorr = ssMotCorr.TRes;

    TMotB0Corr = ssMotB0Corr.TRes;

    limPlot = 1.1*[repmat( multDimMax(abs(dynInd(TGT,1:3,6))),[1 2]) ;
               repmat(convertRotation(multDimMax(abs(dynInd(TGT,4:6,6))),'rad','deg'),[1 2]) ] .* repmat([-1 1],[2 1]);

    limPlot= dynInd(limPlot,1,1, 1.5*dynInd(limPlot,1,1));

    import java.awt.Robot %Tool for generating mouse click as visMotion requires this to continue the script
    import java.awt.event.*
    keys = Robot;
 
    visMotion(TGT,1.5,[],2,[],[],limPlot,[],1,456);
    saveFig(fullfile(saveDir,sprintf('FigureS1_Tr_GT%s',suff)),[],saveFlag,[],1,imageResol,[],imageType);
    
    keys.keyPress(java.awt.event.KeyEvent.VK_SPACE ); keys.keyRelease(java.awt.event.KeyEvent.VK_SPACE );

    visMotion(TMotCorr,1.5 ,[],2,[],[],limPlot,[],1,457);
    saveFig(fullfile(saveDir,sprintf('FigureS1_Tr_MotCorr%s',suff)),[],saveFlag,[],1,imageResol,[],imageType);
    
    keys.keyPress(java.awt.event.KeyEvent.VK_SPACE ); keys.keyRelease(java.awt.event.KeyEvent.VK_SPACE );

    visMotion(TMotB0Corr,1.5,[],2,[],[],limPlot,[],1,458);
    saveFig(fullfile(saveDir,sprintf('FigureS1_Tr_MotB0Corr%s',suff)),[],saveFlag,[],1,imageResol,[],imageType);

    keys.keyPress(java.awt.event.KeyEvent.VK_SPACE ); keys.keyRelease(java.awt.event.KeyEvent.VK_SPACE );
           
end
