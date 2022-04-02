
clc; clear all; close all;
currentDir = mfilename('fullpath');
cd(fileparts(currentDir));

%% SET PATHS
addpath(genpath('..'));
pathInfo; 
figureSpecs;

%% LOAD 3T AND 7T DATA
%%% 7T
pathData = fullfile( mainDir, 'Results','Exp1_Pose','HV1','7T');
rrx_7T = [0 200000];%Display range of anatomical images
ss = load(fullfile(pathData ,'B0Model.mat'));
D_7T = ss.D;%LC maps
x_7T= rescaleND( ss.x, [0 1], rrx_7T);%Rescale so that registration can be performed between 3T and 7T
MT_7T = ss.MT;MS_7T = ss.MS;%sform and resolution

%%% 3T
pathData =  fullfile( mainDir, 'Results','Exp1_Pose','HV1','3T');
rrx_3T = [0 200];
ss = load(fullfile(pathData ,'B0Model.mat'));
D_3T = ss.D;
x_3T= rescaleND( ss.x, [0 1], rrx_3T);
MT_3T = ss.MT;MS_3T = ss.MS;

%% REGISTER BOTH LC MAPS
xi =  {x_3T, x_7T;...
       D_3T, D_7T};
orient = {MT_3T,MT_7T};
res = {MS_3T,MS_7T};
if gpu;xi = gatherStruct(xi,0);end
xo = alignVolumes(xi, res, orient, 1, 3);

%% PLOT
%%% 3T LC maps 
rr = [-10 10]*(3/7);%Scale LC maps display range to account for field strength
Mmagn = abs(dynInd(abs(xo),{1,1},4:5))>0.08;
plotND([],permute(dynInd(real(xo),{2:3,1},4:5),[1:5])/(180/pi),rr ,[],1,{[],2},MT_3T,[],Mmagn,{3},[]);
hColorbar = findobj(gcf,'Type','Colorbar'); set(get(hColorbar,'label'),'string','[$Hz/ ^\circ$]','FontSize',20,'Interpreter','latex','Linewidth',15);
saveFig(fullfile(saveDir,'Figure1_LCM3T'),[],saveFlag,[],1,imageResol,[],imageType);

%%% 7T LC maps 
rr = [-10 10]*(7/7);
Mmagn = abs(dynInd(abs(xo),{1,2},4:5))>0.15;
plotND([],permute(dynInd(real(xo),{2:3,2},4:5),[1:5])/(180/pi),rr ,[],1,{[],2},MT_3T,[],Mmagn,{3},[]);
hColorbar = findobj(gcf,'Type','Colorbar'); set(get(hColorbar,'label'),'string','[$Hz/ ^\circ$]','FontSize',20,'Interpreter','latex','Linewidth',15);
saveFig(fullfile(saveDir,'Figure1_LCM7T'),[],saveFlag,[],1,imageResol );

