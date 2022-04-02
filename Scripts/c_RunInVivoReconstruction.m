
clc; clear all; close all;
currentDir = mfilename('fullpath');
cd(fileparts(currentDir));

%% SET PATHS
addpath(genpath('..'));
pathInfo; 

%% RUN RECONSTRUCTIONS
for HV = [1 2] %Run each Healthy Volunteer (HV)

    if HV==1
       scansToRecon = {'Low','Medium','High'};%No motion free scan in this experiment
    else
       scansToRecon = {'Low','Medium','High','MotionFree'}; 
    end
    
    for i=1:length(scansToRecon) %Run each experiment
    
        %%% LOAD DATA
        fileName = fullfile(mainDir,'Data','Exp2_Motion',sprintf('HV%d',HV), sprintf('%s.mat',scansToRecon{i}));
        load(fileName);%This is a prepared reconstruction structure with all the data needed to reconstruct
        rec.Names.pathOu = fullfile(mainDir,'Results','Exp2_Motion',sprintf('HV%d',HV));%In this folder, a subfolder An-Ve (with all the reconstructions) and An-Ve_Log (with the logging information) will be created
        
        %%% SET PARAMETERS 
        rec.Alg.WriteSnapshots=0;
        rec.Alg.parXT.writeInter=0;
        rec.Alg.parXT.perc=[0.9 0.9 0.9];%Energy to be preserved after coil compression (if strictly lower than one, otherwise number of components preserved) / Redundant components used for motion estimation / Energy preserved when estimating motion
        rec.Alg.parXT.meanT=0;
        rec.Alg.UseSoftMasking=0;
        rec.Dyn.Debug=0;%Set flag to 2 to see results in every iteration.

        rec.Dyn.GPUbS=[6 7];
        rec.Dyn.MaxMem=[6e6 10e6 1e6];%[6e6 2e6 1e6];%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing

        %%% MOTION CORRECTION
        rec.Plan.Suff='_MotCorr';rec.Plan.SuffOu='';
        rec.Alg.AlignedRec=2;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections

        rec.Alg.parXT.groupSweeps=2;% Factor to group sweeps
        rec.Alg.parXT.redFOV=.4;
        rec.Alg.parXT.traLimXT=[0.03 0.01];
        rec.Alg.parXT.disableGrouping = 1;
        rec.Alg.parXT.convTransformJoint = 1;
        
        rec.Alg.resPyr = [.25 0.5 1];
        rec.Alg.nExtern = [10 40 15];
        
        solveXTD(rec); 

        %%% (MOTION + B0) CORRECTION
        rec.Alg.parXB.useTaylor = 1;%First order Taylor model
        rec.Plan.Suff='_MotB0Corr';rec.Plan.SuffOu='';
        
        rec.Alg.parXB.nGD_Taylor = 2;
        rec.Alg.parXB.C.filterTaylor.sp = 1.5;%In mm
        rec.Alg.parXB.C.filterTaylor.gibbsRinging = 0.1;
        rec.Alg.parXB.C.unwrapTaylor=1;
        
        rec.Alg.parXB.taylorDelay=6;
        rec.Alg.parXB.redFOV = 1/(3);

        rec.Alg.parXB.Optim.alphaList=10.^(-1* [-3:0.7:6] + 4.5);
        rec.Alg.parXB.Optim.weightedEst=1;

        solveXTD(rec);close all;
    end
        
end




