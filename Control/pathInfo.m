
%This script contains the location of the main folder on the system and creates the necessary sub-folders.

mainDir = '/home/ybr19/Publications/Journal/2021_07_MRM_B0InformedDISORDER/Repository/';
%example: mainDir = '/home/user/Code/B0InformedDISORDER/'

dirToCreate  = [];
dirToCreate{1} = fullfile(mainDir,'Results','Exp1_Pose','HV1','3T'); 
dirToCreate{2} = fullfile(mainDir,'Results','Exp1_Pose','HV1','7T'); 
dirToCreate{3} = fullfile(mainDir,'Results','Exp2_Motion','HV1'); 
dirToCreate{4} = fullfile(mainDir,'Results','Exp2_Motion','HV2'); 
dirToCreate{5} = fullfile(mainDir,'Results','Simulations'); 
dirToCreate{6} = fullfile(mainDir,'Manuscript','Figures'); 

for i = 1:length(dirToCreate)
    if ~exist(dirToCreate{i},'dir'); mkdir( dirToCreate{i}) ;end
end

gpu = gpuDeviceCount();%Variable to detect gpu functionality. This will speed up generating the figures.

%%
dirToRemove = [];
dirToRemove{1} = fullfile(mainDir,'Data_NotForPublication'); 
dirToRemove{2} = fullfile(mainDir,'Scripts_NotForPublication'); 
dirToRemove{3} = fullfile(mainDir,'Methods_NotForPublication'); 
dirToRemove{4} = fullfile(mainDir,'Manuscript_NotForPublication'); 

for i = 1:length(dirToRemove)
    if exist(dirToRemove{i},'dir'); rmpath( genpath(dirToRemove{i})) ;end
end