
clc; clear all; close all;
currentDir = mfilename('fullpath');
cd(fileparts(currentDir));

%% SET PATHS
addpath(genpath('..'));
pathInfo; 

%% RUN 
for H0 = [3 7] %Field strength
    
    fprintf('\nFitting LC maps on the %dT pose experiment.\n',H0);
    
    if H0==7
        rrx = [0 150000];th=.2;%Dynamic range of the images and the relative threshold to use for tissue extraction (see figure captions)
    elseif H0==3
        rrx = [0 250];th=.2;   
    end
    
    pathData = fullfile( mainDir, 'Data','Exp1_Pose','HV1',sprintf('%dT',H0));
    pathResults = fullfile( mainDir, 'Results','Exp1_Pose','HV1',sprintf('%dT',H0)); 
    
    xThr = rrx(2)*th;%Threshold to extract tissue 
    
    %% LOAD REGISTRATION INFORMATIO
    load(fullfile(pathData, 'RegistrationInfo.mat'));%Load registration parameters needed to fit the B0 model
    %This load the variable "eulerAngles" for all the poses
    
    %% LOAD REFERENCE IMAGE
    fileName = fullfile(pathData,'NIFTI', sprintf('Pose%d',1));%Load the reference pose
    [x,MS,MT] = readNII(fileName,{'x_Registered'});
    x=x{1}; %Image
    MS=MS{1};%Resolution
    MT=MT{1};%s_form
    N=size(x);%Array size

    plotND([],abs(x),rrx,[],0,{[],2},MT);

    %% LOAD ALL POSES
    NPoses = size(eulerAngles,2);%Number of scans
    for pose = 1:NPoses
        %%% Read registered nifti image
        fprintf('Reading induced B0 for pose %d \n', pose)        
        fileName = fullfile(pathData,'NIFTI', sprintf('Pose%d',pose));
        [xx] = readNII(fileName,{'dB0_Registered'}); %Induced B0 in the head frame in Hz
        %%% Store in 4 dimensional array
        if pose == 1; f = xx{1}; else f = cat(4,f,xx{1});end
    end

    %% FIT LC MAPS
    %%% Fit LC maps using the model F = A*D with A the rotaton angles,D the Linear Coefficient maps (LC maps) and F the induced fields
    
    %%% Get pitch and roll rotation angles from the registration parameters
    pitch = convertRotation( eulerAngles(1,:).', 'deg','rad');%Convert to radians
    roll = convertRotation( eulerAngles(2,:).', 'deg','rad');

    A = zeros(8,2); % First order Taylor model
    A(:,1) = pitch ;
    A(:,2) = roll ;

    %%% Flatten to do fit for every voxel simultaneously
    F = resSub(f,1:3)';
    
    %%% Least-Squares fit
    D = A\F; 
    FResidual = F - A*D;%Residuals 
    
    %%% Reshape to an image size
    D = resSub(D',1,N);%These are the fitted LC maps    
    fResidual = resSub(FResidual',1,N);
    
    %% COMPUTE METRICS
    RMSE = sqrt(multDimMea( fResidual.^2 , 4)) ;
    MAE = multDimMea( abs(fResidual) , 4) ;
    R2 = 1 - multDimSum( fResidual.^2 , 4)./ multDimSum(   bsxfun(@minus,f, multDimMea(f,4)).^2   , 4  ) ;
   
    %% SAVE
    D = single(gather(D));%LC maps
    x = single(gather(x));%Image from the reference pose (useful for simulations as this is anatomically aligned with the LC maps)
    RMSE =  single(gather(RMSE));
    MAE =  single(gather(MAE));
    R2 = single(gather(R2));
   
    save(fullfile(pathResults ,'B0Model.mat') , 'D','x','MS','MT','R2','RMSE','MAE', '-v7.3')

end

