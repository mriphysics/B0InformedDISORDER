
%This script contains plotting parameters used for the In-Vivo results.

rrD = [-10 10]; %Dynamic range of the LC maps in Hz/degree

if HV==1
    nameB0MF = strcat(mainDir, '/Results/Exp1_Pose/HV1/7T/B0Model.mat');ss = load(nameB0MF);%This is from the pose experiment (=different acquisition)
    %Load LC maps and associated image and orientation information
    DGT = ss.D;%LC maps
    xDGT = ss.x;%Anatomical image
    MSDGT = ss.MS ;%Resolution
    MTDGT = ss.MT;%sform
    
    %Load GT image and orientation information
    xGT = ss.x;
    MSxGT = ss.MS ;
    MTxGT = ss.MT;
    clear ss
    
    %Range of the motion parameters
    limTraces = {[-1.2 1.2; -1.7 1.7],[-1.7 1.7; -3.5 3.5],[-7.5 7.5; -14 14]};
    %Dynamic range of the image
    rrX= [5 170]*1e3;
    th=0.2;%Relative threshold to extract tissue
    
    %ROI to extract
    supFOV={ 12:144 ,[],[]};%Support of the FOV to display
    idxNIFTI = [68 68 75]; %Index of the image slices to plot in logical units
    idxNIFTI_LCM = [90 64 70];  %Index of the LC map slices to plot in logical units
    
elseif HV==2

    nameMF = strcat(mainDir, '/Results/Exp2_Motion/HV2/An-Ve/MotionFree');
    [xGT,MSxGT,MTxGT] = readNII(nameMF, {'Di_MotCorr'});
    xGT=xGT{1};MSxGT=MSxGT{1};MTxGT=MTxGT{1};
    %Motion range
    limTraces = {[-.7 .7; -1.2 1.2],[-1.4 1.4; -2.1 2.1],[-4.9 4.9; -11 11]};
    %Dynamic range
    rrX= [15 210]*1e3;
    th=0.2;
    %ROI to extract
    supFOV={ 10:141 ,5:169, 14:134};
    idxNIFTI = [90    88    74]; 
    idxNIFTI_LCM = [76 78 82]; 

end
