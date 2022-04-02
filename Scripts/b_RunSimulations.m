
clc; clear all; close all;
currentDir = mfilename('fullpath');
cd(fileparts(currentDir));

%% SET PATHS
addpath(genpath('..'));
pathInfo; 

%% LOAD DATA
%%% Load reconstruction structure
recName = fullfile(mainDir,'Data','Exp1_Pose','HV1','7T','Pose1.mat');%Name where an example reconstruction structure is stored which has the necessary structure to be called by the reconstruction
ssRec = load(recName,'rec');rec= ssRec.rec;clear ssRec;%ss is the naming convention used for an intermediate matlab structure when loading files
rrx = [0  150000]; th = .2;%Display range and relative threshold value for plotting
xTh=rrx(2)*th;

%%% Load both image and LCMs that are anatomically aligned
dName = fullfile(mainDir,'Results','Exp1_Pose','HV1','7T','B0Model.mat');%Name where the fitted LC maps are stored
ssD = load(dName);
xGT = ssD.x;%Anatomical image
DGT = ssD.D; %LC maps in Hz/rad
DGT = bsxfun(@times, DGT, abs(xGT)>xTh);%Filter LC maps as the fit outside tissue areas is noisy (and meaningless)
MT = ssD.MT;%s_form
clear ssD;

%%% Plot to inspect
M = abs(xGT)>xTh;
plotND({abs(xGT),1,0}, DGT/180*pi, [[-10 10] rrx],[],1,[],MT,{'Anatomical reference';'$d_{pitch}$';'$d_{roll}$'},M,{3},[],'Image and LC maps used for simulations.','[$Hz/ ^\circ$]');

%% SET SIMULATION PARAMETERS
%%% General
rec.Plan.Suff='_Sim';rec.Plan.SuffOu='Sim';
rec.Alg.WriteSnapshots=0;
rec.Alg.AlignedRec=1;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
rec.Alg.UseSoftMasking= 0;
rec.Dyn.GPUbS = [6 7];%GPU block sizes used in the recsontruction. Depending on the machine used to run these simulations, this could require other values.
rec.Dyn.MaxMem=[10e6 10e6 10e6];%%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing
rec.Dyn.Log=0;
rec.Dyn.Debug=0;%Set flag to 2 to see results in every iteration.

%%% Parameters for motion correction Reconstruction
rec.Alg.parXT.writeInter=0;
rec.Alg.parXT.meanT=1;%Enforce mean transformation to be zero, as we will also generate zero-mean motion traces.
rec.Alg.parXT.redFOV=0.0; %Disable FOV reduction when estimating motion as pure rigid motion is simulated.
rec.Alg.parXT.convTransformJoint = 1;%To estimated all motion parameters in every iteration

%%% Parameters for the B0 model used
rec.Alg.parXB.useTaylor=1;%To indicate the use of the Taylor B0 model
rec.Alg.parXB.taylorDelay=15;%Number of iterations after which to activate the LC maps estimation
rec.Alg.parXB.nGD_Taylor=2;%Number iterations of the Gradient-Desccent (GD) solver within each outer iteration
rec.Alg.parXB.Optim.alphaList=10.^(-1* [-10:.7:6]);%Step size for the GD based LC map estimation
rec.Alg.parXB.C.unwrapTaylor = 1;%Whether to use the unwrapping strategy to filter LC maps.

%%% Simulations: Data synthesis
rec.Sim.rrx=rrx;%Used for plotting only
rec.Sim.xGT = xGT;
rec.Sim.DGT = DGT;
rec.Sim.disableGrouping=1;
rec.Sim.motionType = 'interleavedExponential';%See generateTransformTrace.m for other motion types
rec.Sim.expFact = 0.1;
rec.Sim.NinterleavesT = 15;%Time of each interleave (in seconds)
rec.Sim.tran = 10;%Range of translation of the simulated motion trace (in mm)
rec.Sim.rot = 15;%%Range of rotations of the simulated motion trace (in degrees)
rec.Sim.meanT = 1;%Simulate motion trace with zero mean
rec.Sim.useMask = 0;

%%% Simulations: reconstruction
rec.Sim.blockConvergence=1;%Avoid early convergence
rec.Sim.nExtern=[45 15];%Number of outer iterations at every resolution level
rec.Sim.resRecon=[.5 1];%Resolution levels for reconstruction
rec.Sim.resSynth = rec.Sim.resRecon(end);%Resolution at which to simulate k-space - set it to the finest resolution

rec.Alg.parXT.groupSweeps = 2;%Grouping of sweeps for data reconstruction. Note that the #sweeps is determined by the sampling trajectory as provided in ../Data/Exp1_Pose/HV1/7T/Pose1.mat
rec.Sim.groupSweepsSynth = rec.Alg.parXT.groupSweeps;%Grouping of sweeps for data synthesis. Set it to the same number of reconstruction sweeps to avoid intra-shot motion.

%% RUN SIMULATIONS

for SNRdB = [30 inf] %Respectively 30dB and no noise.
    
    rec.Sim.snrdB = SNRdB; %[inf 30 20 10] ==> No noise / High SNR / Acceptable SNR / Poor SNR

    %Uncorrected
    rec.Sim.synthX = 1;     rec.Sim.synthT = 1;     rec.Sim.synthD = 1; %Which parts of the signal model to activate when simulating k-space
    rec.Sim.provideX = 0;   rec.Sim.provideT = 0;   rec.Sim.provideD = 0; %Which parts of the signal model to provide the GT values for (and therefore not estimating)
    rec.Sim.estX = 1;       rec.Sim.estT = 0;       rec.Sim.estD = 0; %Which parts of the signal model to estimate when reconstructing    
    rec.Sim.nameSave = fullfile(mainDir,'Results','Simulations',simulationName(rec));%Generate a filename that incorporates the different simulation options
    controlledSimulation(rec);%Run the simulations and save results

    %Motion corrected
    rec.Sim.estX = 1;       rec.Sim.estT = 1;       rec.Sim.estD = 0;     
    rec.Sim.nameSave = fullfile(mainDir,'Results','Simulations',simulationName(rec));
    controlledSimulation(rec);

    %(Motion + B0) corrected
    rec.Sim.estX = 1;       rec.Sim.estT = 1;      rec.Sim.estD = 1;     
    rec.Sim.nameSave = fullfile(mainDir,'Results','Simulations',simulationName(rec));
    controlledSimulation(rec);

    %(Motion + B0) provided
    rec.Sim.provideX = 0;   rec.Sim.provideT = 1;   rec.Sim.provideD = 1; 
    rec.Sim.estX = 1;       rec.Sim.estT = 0;       rec.Sim.estD = 0;     
    rec.Sim.nameSave = fullfile(mainDir,'Results','Simulations',simulationName(rec));
    controlledSimulation(rec);

    %Motion-free
    rec.Sim.synthX = 1;     rec.Sim.synthT = 0;     rec.Sim.synthD = 0;   
    rec.Sim.provideX = 0;   rec.Sim.provideT = 0;   rec.Sim.provideD = 0; 
    rec.Sim.estX = 1;       rec.Sim.estT = 0;       rec.Sim.estD = 0;     
    rec.Sim.nameSave = fullfile(mainDir,'Results','Simulations',simulationName(rec));
    controlledSimulation(rec);

end

