clc; clear all; close all;
currentDir = mfilename('fullpath');
cd(fileparts(currentDir));

%% SET PATHS
addpath(genpath('..'));
pathInfo; 
figureSpecs;

experimentNames = {'Low','Medium','High'};
suff={'MotCorr', 'MotB0Corr'};
removeLegend = 1;

%% ITERATE OVER HEALTHY VOLUNTEERS
for sameRange=[1]
    for HV=[1  2]
    %%% Load in-vivo parameters
    InViVoParams;

    %%% Plot motion traces
        for f = 1:length(experimentNames)
            for s = 1:length(suff)
                fileNameMotion =fullfile(mainDir,'Results','Exp2_Motion',sprintf('HV%d',HV),'An-Ve',sprintf('%s_Tr_%s.mat',experimentNames{f},suff{s}) );

                load(fileNameMotion)
                time = multDimMea(MotionInfo.timeState,2:ndims(MotionInfo.timeState));
                if sameRange
                    limTracesToUse=limTraces{end};
                else
                    limTracesToUse=limTraces{f};
                end
                MS = sqrt(sum(MotionInfo.Par.Mine.APhiRec(1:3,1:3).^2,1));
                visMotion(MotionInfo,MS,time,2,[],[],limTracesToUse,[],1)
                set(gcf, 'Position', get(0,'Screensize'),'Color',[1 1 1]);

                %Remove legend
                if removeLegend; hLeg=findobj(gcf,'type','legend');set(hLeg,'visible','off');end

                if sameRange
                    saveFig(fullfile(saveDir,sprintf('Figure5_MotionTrace_HV%d_%s_%s_SameRange',HV,experimentNames{f}, suff{s})),[],saveFlag,[],1,imageResol,[],imageType);
                else
                    saveFig(fullfile(saveDir,sprintf('Figure5_MotionTrace_HV%d_%s_%s',HV,experimentNames{f}, suff{s})),[],saveFlag,[],1,imageResol,[],imageType);
                end
                close 
            end
        end
    end
end
