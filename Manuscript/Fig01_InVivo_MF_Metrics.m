
clc; clear all; close all;
currentDir = mfilename('fullpath');
cd(fileparts(currentDir));

%% SET PATHS
addpath(genpath('..'));
pathInfo; 
figureSpecs;

%% LOAD 3T AND 7T DATA
%%% 7T
pathData = fullfile(mainDir,'Results','Exp1_Pose','HV1','7T');
H0 = 7; rrx_7T = [0 200000];
ss = load(fullfile(pathData ,'B0Model.mat'));
R2_7T = ss.R2;
RMSE_7T = ss.RMSE;
MAE_7T = ss.MAE;
x_7T= rescaleND( ss.x, [0 1], rrx_7T);
MT_7T = ss.MT;MS_7T = ss.MS;

%%% 3T
pathData = fullfile(mainDir,'Results','Exp1_Pose','HV1','3T');
rrx_3T = [0 200];
ss = load(fullfile(pathData ,'B0Model.mat'));
R2_3T = ss.R2;
RMSE_3T = ss.RMSE;
MAE_3T = ss.MAE;
x_3T= rescaleND( ss.x, [0 1], rrx_3T);
MT_3T = ss.MT;MS_3T = ss.MS;

%% REGISTER BOTH LC MAPS
xi =  {x_3T,                         x_7T;...
       cat(4,R2_3T,RMSE_3T,MAE_3T) , cat(4,R2_7T,RMSE_7T,MAE_7T)};
orient = {MT_3T,MT_7T};
res = {MS_3T,MS_7T};
if gpu;xi = gatherStruct(xi,0);end
xo = alignVolumes(xi, res, orient, 1, 3);

%% PLOT
%%% R2-SCORE
rr= [0 1];
Mmagn = abs(dynInd(abs(xo),{1,1},4:5))>0.08;
plotND([],permute(dynInd(real(xo),{2,1},4:5),[1:3 5 4]),rr ,[],0,{[],2},MT_3T,[],Mmagn,{1},[]);
hColorbar = findobj(gcf,'Type','Colorbar'); set(get(hColorbar,'label'),'string','[a.u.]','FontSize',20,'Interpreter','latex','Linewidth',15);
saveFig(fullfile(saveDir,'Figure1_R2_3T'),[],saveFlag,[],1,imageResol );

Mmagn = abs(dynInd(abs(xo),{1,2},4:5))>0.2;
plotND([],permute(dynInd(real(xo),{2,2},4:5),[1:3 5 4]),rr ,[],0,{[],2},MT_3T,[],Mmagn,{1},[]);
hColorbar = findobj(gcf,'Type','Colorbar'); set(get(hColorbar,'label'),'string','[a.u.]','FontSize',20,'Interpreter','latex','Linewidth',15);
saveFig(fullfile(saveDir,'Figure1_R2_7T'),[],saveFlag,[],1,imageResol );

%%% RMSE
rr_RMSE = [0 15];
rr= rr_RMSE*(3/7);
Mmagn = abs(dynInd(abs(xo),{1,1},4:5))>0.08;
plotND([],permute(dynInd(abs(xo),{3,1},4:5),[1:3 5 4]),rr ,[],0,{[],2},MT_3T,[],Mmagn,{1},[]);
hColorbar = findobj(gcf,'Type','Colorbar'); set(get(hColorbar,'label'),'string','[Hz]','FontSize',20,'Interpreter','latex','Linewidth',15);
saveFig(fullfile(saveDir,'Figure1_RMSE_3T'),[],saveFlag,[],1,imageResol );

rr= rr_RMSE*(7/7);
Mmagn = abs(dynInd(abs(xo),{1,2},4:5))>0.2;
plotND([],permute(dynInd(abs(xo),{3,2},4:5),[1:3 5 4]),rr ,[],0,{[],2},MT_3T,[],Mmagn,{1},[]);
hColorbar = findobj(gcf,'Type','Colorbar'); set(get(hColorbar,'label'),'string','[Hz]','FontSize',20,'Interpreter','latex','Linewidth',15);
saveFig(fullfile(saveDir,'Figure1_RMSE_7T'),[],saveFlag,[],1,imageResol );

%%% MAE
rr_MAE = [0 15];
rr= rr_MAE*(3/7);
Mmagn = abs(dynInd(abs(xo),{1,1},4:5))>0.08;
plotND([],permute(dynInd(abs(xo),{4,1},4:5),[1:3 5 4]),rr ,[],0,{[],2},MT_3T,[],Mmagn,{1},[]);
hColorbar = findobj(gcf,'Type','Colorbar'); set(get(hColorbar,'label'),'string','[Hz]','FontSize',20,'Interpreter','latex','Linewidth',15);
saveFig(fullfile(saveDir,'Figure1_MAE_3T'),[],saveFlag,[],1,imageResol );

rr= rr_MAE*(7/7);
Mmagn = abs(dynInd(abs(xo),{1,2},4:5))>0.2;
plotND([],permute(dynInd(abs(xo),{4,2},4:5),[1:3 5 4]),rr ,[],0,{[],2},MT_3T,[],Mmagn,{1},[]);
hColorbar = findobj(gcf,'Type','Colorbar'); set(get(hColorbar,'label'),'string','[Hz]','FontSize',20,'Interpreter','latex','Linewidth',15);
saveFig(fullfile(saveDir,'Figure1_MAE_7T'),[],saveFlag,[],1,imageResol );

