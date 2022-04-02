
function [y,TGT,TGTFullTempRes] = synthesizeKSpace(rec)

%SYNTHESISEYDISORDER  synthesises k-space for a set of synthesis options.
%The code is based on the function controlledSimulation.m to have the same
%data structures needed for the reconstruction.
%
disp('=======================  START DATA SYNTHESIS ===============================');

%% Set grouping differently
rec.Alg.parXT.groupSweeps = rec.Sim.groupSweepsSynth;
rec.Dyn.Debug = 0;

%% Same code as in SOLVEXTB.m
tsta=tic;

%SET DEFAULT PARAMETERS
rec=disorderAlgorithm(rec);
rec=B0Algorithm(rec);

% %FILENAME
folderSnapshots=strcat(rec.Names.pathOu, filesep, 'An-Ve_Sn');
[~,fileName]=fileparts(generateNIIFileName(rec));

c = clock;

%SHORTCUTS FOR COMMONLY USED STRUCTURE FIELDS AND INITIALIZERS
voxSiz=rec.Enc.AcqVoxelSize;%Acquired voxel size
parXT=rec.Alg.parXT;parXB=rec.Alg.parXB;
gpuIn=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);if gpuIn;gpuFIn=2;else gpuFIn=0;end 
gpu=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);

if size(rec.y,5)>=parXT.maximumDynamics;fprintf('Problem too big to fit in memory (%d dynamics for a limit of %d)\n',size(rec.y,5),parXT.maximumDynamics);rec.Fail=1;return;end
on=cell(1,5);for n=1:5;on{n}=ones(1,n);end
typ2Rec=rec.Dyn.Typ2Rec;

%RECONSTRUCTION PLAN
%[resPyr,~,~,resIso,~]=pyramidPlan(voxSiz,parXT.resolMax,parXT.NWend,parXT.accel);
resPyr = rec.Sim.resSynth; %YB: Change the resolution to the one for data-synthesis: rec.Sim.resSynth
resIso=sqrt(2)*((prod(voxSiz).^(1/3))./resPyr);%resIso = [resIso(end-(length(resPyr)-1):end)]
L = length(resIso);
estT=ones([1 L]);
estB = ones([1 L]);

fprintf('   Resolution:         %s \n',sprintf(' %g ',resPyr))

%ROI COMPUTATION AND EXTRACTION IN THE READOUT DIRECTION
rec.Enc.ROI=computeROI(onesL(rec.M),[],[0 0 0],[1 0 0]); %YB: At this point non-logical values
if rec.Dyn.Debug>=2
    fprintf('ROI:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));
    fprintf('Acquisition directions: %s\n',rec.Par.Scan.MPS);
end
for n=typ2Rec';datTyp=rec.Plan.Types{n};
    if ~ismember(n,[5 12]);rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,1);end
end

%PERMUTE PHASE ENCODES TO THE FIRST DIMENSIONS
perm=1:rec.Plan.NDims;perm(1:3)=[3 2 1]; 
for n=typ2Rec'        
    if ~ismember(n,[5 12]);datTyp=rec.Plan.Types{n};
        rec.(datTyp)=permute(rec.(datTyp),perm);
    end
end
voxSiz=voxSiz(perm(1:3));parXT.apod=parXT.apod(perm(1:3));
if ~isfield(rec.Par.Mine,'permuteHist'); rec.Par.Mine.permuteHist = [];end
rec.Par.Mine.permuteHist{end+1} = perm(1:4);

%COIL ARRAY COMPRESSION AND RECONSTRUCTED DATA GENERATION
NX=size(rec.M);
[S,y,eivaS]=compressCoils(rec.S,parXT.perc,rec.y);

NS=size(S);
if gpuIn;y=gpuArray(y);end
if rec.Dyn.Debug>=2 && ~isempty(parXT.perc)
    if parXT.perc(1)<1;fprintf('Number of compressed coil elements at%s%%: %d (%s )\n',sprintf(' %0.2f',parXT.perc*100),NS(4),sprintf(' %d',eivaS));else fprintf('Number of compressed coil elements: %d (%s )\n',NS(4),sprintf(' %d',eivaS));end
    fprintf('Number of coils used for intermediate image estimation: %d\n', eivaS(1));
    fprintf('Number of coils used for motion estimation: %d\n', eivaS(2));
    fprintf('Number of coils used for final image estimation: %d\n', eivaS(3));
end

%APODIZE + REARRANGE DATA + CORRECT FOR INTERLEAVED REPEATS 
NY=size(y);NY(end+1:rec.Plan.NDims)=1;
y=bsxfun(@times,y,ifftshift(buildFilter(NY(1:3),'tukey',[10 10 1],gpuIn,parXT.apod))); %Apodize 
[y,NY]=resSub(y,5:rec.Plan.NDims);NY(end+1:5)=1;%Different averages along the 5th dimension
y=gather(y);

%ACQUISITION STRUCTURE
NEchos=max(rec.Par.Labels.TFEfactor,rec.Par.Labels.ZReconLength);
NProfs=numel(rec.Assign.z{2});
if NProfs<=1 || NProfs~=numel(rec.Assign.z{3});fprintf('SEPARABLE Y-Z TRAJECTORIES. ALIGNED RECONSTRUCTION IS NOT PERFORMED\n');rec.Fail=1;return;end
if sum(cdfFilt(abs(diff(rec.Assign.z{2}(:))),'med'))<sum(cdfFilt(abs(diff(rec.Assign.z{3}(:))),'med'));fprintf('Slow direction is: %s\n',rec.Par.Scan.MPS(4:5));else fprintf('Slow direction is: %s\n',rec.Par.Scan.MPS(7:8));end
%YB: z{2} = AP and z{3} = LR
kTraj=zeros([NProfs 2],'single');
for n=1:2;kTraj(:,n)=rec.Assign.z{perm(n)}(:);end

if NEchos==1;NEchos=NProfs;end
NShots=NProfs/NEchos;
if mod(NShots,1)~=0;fprintf('NUMBER OF SHOTS: %.2f, PROBABLY AN INTERRUPTED SCAN\n',NShots);rec.Fail=1;return;end
NRepeats=NY(5);
if NShots<NRepeats
    NShots=NRepeats;
    NProfs=NShots*NEchos;
    kTraj=repmat(kTraj,[NRepeats 1]);
    NEchos=NProfs/NShots;
end
if rec.Dyn.Debug>=2
    fprintf('Number of echoes: %d\n',NEchos);
    fprintf('Number of shots: %d\n',NShots);
end
kTrajSS=reshape(kTraj,[NEchos NShots 2]);
%PLOT TRAJECTORIES
if rec.Alg.WriteSnapshots;visTrajectory(kTrajSS,2,strcat(folderSnapshots,filesep,'Trajectories'),strcat(fileName,rec.Plan.SuffOu));end
kRange=single(zeros(2,2));
for n=1:2;kRange(n,:)=rec.Enc.kRange{perm(n)};end
kShift=(floor((diff(kRange,1,2)+1)/2)+1)';
kSize=(diff(kRange,1,2)+1)';
kIndex=bsxfun(@plus,kTraj,kShift);
rec.Par.Mine.DISORDER.kTraj=gather(kTraj);
rec.Par.Mine.DISORDER.kRange=gather(kRange);
rec.Par.Mine.DISORDER.kShift=gather(kShift);
rec.Par.Mine.DISORDER.kSize=gather(kSize);
rec.Par.Mine.DISORDER.kIndex=gather(kIndex);

%STEADY-STATE CORRECTIONS
isSteadyState=0;
if mod(NShots,NRepeats)==0
    NShotsUnique=NShots/NRepeats;  
    isSteadyState=(NShotsUnique==1);%No shots have been identified
end
fprintf('Steady-state: %d\n',isSteadyState);
if isSteadyState
    kTrajs=abs(fct(kTraj)); 
    N=size(kTrajs);
    [~,iAs]=max(dynInd(kTrajs,1:floor(N(1)/2),1),[],1);
    iM=min(iAs,[],2);
    NSamples=NProfs/(iM-1);
    NSweeps=round(NProfs/(NSamples*parXT.groupSweeps));
    NSamples=NProfs/NSweeps;
    fprintf('Number of sweeps: %d\n',NSweeps);
else
    NSamples=NEchos;
    NSweeps=NShots;
end
NStates=NSweeps;
%SWEEP SUBDIVISION AND TIME INDEX
sweepSample=ceil(NSweeps*(((1:NProfs)-0.5)/NProfs));
stateSample=sweepSample;
timeIndex=zeros([NY(1:2) NY(5)]);
if gpu;timeIndex=gpuArray(timeIndex);end
for n=1:size(kIndex,1)
    for s=1:NY(5) 
        if timeIndex(kIndex(n,1),kIndex(n,2),s)==0
            timeIndex(kIndex(n,1),kIndex(n,2),s)=n;
            break;
        end
    end
end
for m=1:2;timeIndex=ifftshift(timeIndex,m);end
if isSteadyState
    timeSample=(0:NProfs-1)*rec.Par.Labels.RepetitionTime(1)/1000;
else
    TPerShot=rec.Par.Labels.ScanDuration/NShots;    
    timeSample=(sweepSample-1)*TPerShot;
    if strcmp(rec.Par.Scan.Technique,'TSE') || strcmp(rec.Par.Scan.Technique,'TIR')
        timeSample=timeSample+repmat(0:NEchos-1,[1 NShots])*(rec.Par.Labels.TE(1)/(1000*NEchos/2));
    else
        timeSample=timeSample+repmat(0:NEchos-1,[1 NShots])*rec.Par.Labels.RepetitionTime(1)/1000;
    end
end
sweepToSample=cell(1,NSweeps);
sampleToSweep=1:NSweeps;
for n=1:NSweeps;sweepToSample{n}=n;end

%CREATE RECONSTRUCTION DATA ARRAY
NX=size(rec.M);NX(end+1:3)=1;
rec.x=zeros(NX,'single');%Uncorrected non-regularized
rec.d=zeros([NX L-sum(resPyr~=1)],'single');%Corrected non-regularized 
if parXT.computeCSRecon;rec.r=zeros([NX L-sum(resPyr~=1)],'single');end %Corrected regularized
x=zeros(NX,'single');if gpuIn;x=gpuArray(x);end
typ2RecI=[12;16]; % 12-> Reconstruction 16-> Volumetric alignment
if parXT.computeCSRecon;typ2RecI=[typ2RecI;18];end %18-> Filtered reconstruction
if parXB.useTaylor;typ2RecF=29;else typ2RecF=[];end %29-> Linear coefficients of Bo fields with motion
for n=typ2RecI'
    if ~any(rec.Dyn.Typ2Rec==n);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,n);rec.Dyn.Typ2Wri(n)=1;end
end
typ2Rec=rec.Dyn.Typ2Rec;

%INITIALIZE TRANSFORM, NUMBER OF EXTERN ITERATIONS, TOLERANCES, IDENTIFY 
%STEADY STATE SEQUENCES 
rec.T=zeros([on{4} NSweeps 6],'single');
subT=ones(1,NStates);%To perform subdivisions of T
nExtern=99;

tolType={'Energy','RelativeResidual2Error'};%For CG--without vs with motion correction
tol=[parXT.tolerSolve 0];%For CG
nIt=[300 1];%For CG
nowithin=0;%To stop estimation of within motion

%ENERGY
En=cell(1,L);Re=cell(1,L);EnX=[];
mSt=cell(1,L);iSt=cell(1,L);
tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time arranging: %.3f s\n',tend);end

%SOLVE WITHOUT MOTION
tsta=tic;
resAni=single(on{3});
[~,indMinRes]=min(voxSiz);        
resAni(indMinRes)=voxSiz(indMinRes)/resPyr(end);
indNoMinRes=1:3;indNoMinRes(indMinRes)=[];
resAni(indNoMinRes)=max(resAni(indMinRes),voxSiz(indNoMinRes));   
resAni=voxSiz./resAni;
BlSz=max(sqrt(prod(voxSiz(1:2)./resAni(1:2))),1);%Block sizes baseline
if ~isfield(rec.Dyn,'GPUbS'); rec.Dyn.GPUbS = [2 4];end
E.bS=round([rec.Dyn.GPUbS].^log2(BlSz));E.oS=on{2};%Block sizes for our GPU
if rec.Dyn.Debug>=2;fprintf('Block sizes:%s\n',sprintf(' %d',E.bS));end
E.Je=1;%To use joint encoding/decoding
E.Uf=cell(1,3);for n=1:3;E.Uf{n}.NY=NY(n);E.Uf{n}.NX=NX(n);end%Folding
EH.Ub=cell(1,3);for n=1:3;EH.Ub{n}.NY=NY(n);EH.Ub{n}.NX=NX(n);end%Unfolding
[E.UAf,EH.UAb]=buildFoldM(NX(1:3),NY(1:3),gpuIn,1);%Folding/Unfolding matrix (empty is quicker!)
R=[];%Regularizer
ncx=1:eivaS(1);

%MASKING
CX.Ma = rec.M;
if rec.Alg.UseSoftMasking==0 && ~all(ismember(CX.Ma(:),[0 1])) % Need to convert rec.M into discrete mask
    %maxint = max(abs(rec.M(:)));  level = graythresh(abs(rec.M ) / maxint);
    maxint =1;  level = 0.45 * mean(abs(rec.M(:)));
    CX.Ma = single( abs(rec.M)> level*maxint);
    CX.Ma = imdilate(CX.Ma, strel('disk', 7));
    rec.M = CX.Ma;
end
if ~rec.Sim.useMask;rec.M = onesL(rec.M);disp('Mask set to ones for data synthesis.');end
CX.Ma = rec.M;

if gpuIn;CX.Ma=gpuArray(CX.Ma);end
if rec.Alg.UseSoftMasking==2;CX.Ma(:)=1;end %YB: 2 means not masking
discMa=all(ismember(CX.Ma(:),[0 1])); %YB: discrete masking - if every element is rather 0 or 1
if ~discMa %Soft masking
    fprintf('Using soft masking\n');
    M=buildFilter(2*[E.Uf{1}.NX E.Uf{2}.NX E.Uf{3}.NX],'tukeyIso',1,gpuIn,1,1);
    CX.Ma=abs(filtering(CX.Ma,M,1));M=[];
    if size(S,6)>1 || isfield(rec,'W')%We do not unfold in the PE direction
        NMa=[E.Uf{1}.NX E.Uf{2}.NX E.Uf{3}.NX];
        for m=1:3
            if ~isempty(E.Uf{m});CX.Ma=fold(CX.Ma,n,E.Uf{m}.NX,E.Uf{m}.NY,E.UAf{m});end
        end         
        CX.Ma=resampling(CX.Ma,NMa,2);
    end
    if isfield(rec,'W');CX.Ma=bsxfun(@times,CX.Ma,rec.W);end
    Ti.la=CX.Ma;
    R.Ti.la=1./(abs(Ti.la).^2+0.01);
    CX=[];
end

%SENSITIVITIES
E.Sf=dynInd(S,ncx,4);E.dS=[1 ncx(end)];
if gpuIn;E.Sf=gpuArray(E.Sf);end

%PRECONDITION
if ~discMa;P.Se=(normm(E.Sf,[],4)+R.Ti.la).^(-1);%Precondition
else P.Se=(normm(E.Sf,[],4)+1e-9).^(-1);
end

gibbsRing=parXT.UseGiRingi*parXT.GibbsRingi;
tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing non-corrected reconstruction: %.3f s\n',tend);end

%ITERATING THROUGH THE DIFFERENT LEVELS
for l=1:L    
    %WE RESET MOTION CONVERGENCE INSIDE
    nReset=0;THist=[]; % YB: setting nReset to zero ensures resetting motion convergence
    leff=l-sum(resPyr~=1);

    %UPSAMPLE THE TRANSFORM AND GENERATE MOTION PARAMETERS   
    rec.T=repelem(rec.T,1,1,1,1,subT,1); %YB: so intra-shot motion is not an extra dimension but are repeated elements in the 5th dimension!
    subTT=repelem(subT,1,subT);
    sampleToSweep=repelem(sampleToSweep,1,subT);
    sweepToSample=cell(1,NSweeps);
    c=0;
    for n=1:length(subT)
        if subT(n)==2
            indState=find(stateSample==(n+c));
            Ni=length(indState);
            indState=indState(floor(Ni/2)+1); % YB: halfway since you split the segment in half
            stateSample(indState:end)=stateSample(indState:end)+1; % YB: so max(stateSample) can change over resolutions
            sweepToSample{sampleToSweep(n+c)}=[sweepToSample{sampleToSweep(n+c)} n+c:n+c+1];%YB: Not exactly similar to sweepSample (not used anymore) and relates the states to the sweeps(stays constant)
            c=c+1;
        else
            sweepToSample{sampleToSweep(n+c)}=[sweepToSample{sampleToSweep(n+c)} n+c];
        end
    end 
    timeState=regionprops(stateSample,timeSample,'MeanIntensity'); %YB: if find(subT(:)==1) is not empty, max( stateSample) > NSweeps
    timeStateLimitInf=regionprops(stateSample,timeSample,'MinIntensity');
    timeStateLimitSup=regionprops(stateSample,timeSample,'MaxIntensity');
    timeState={timeState.MeanIntensity};
    timeStateLimitInf={timeStateLimitInf.MinIntensity};
    timeStateLimitSup={timeStateLimitSup.MaxIntensity};
    timeState=cat(1,timeState{:});
    timeStateLimitInf=cat(1,timeStateLimitInf{:});
    timeStateLimitSup=cat(1,timeStateLimitSup{:});
    if leff<=1;timeSweep=timeState;end
    %NOTE THIS ONLY WORKS IF THERE IS NO INTERSHOT MOTION
    if leff<=1
        if ~isSteadyState;timeState=cat(3,timeStateLimitInf,timeStateLimitSup);else timeState=cat(3,timeStateLimitInf,timeState,timeStateLimitSup);end
    end
    
    %RELATIVE SPATIAL RESOLUTION CALCULATIONS AT THIS LEVEL
    resAni=single(on{3});
    [~,indMinRes]=min(voxSiz);        
    resAni(indMinRes)=voxSiz(indMinRes)/resPyr(l);% If no rounding issues, this should be the 4mm of the coarses resolution
    indNoMinRes=1:3;indNoMinRes(indMinRes)=[];
    resAni(indNoMinRes)=max(resAni(indMinRes),voxSiz(indNoMinRes)); 
    resAni=voxSiz./resAni;
    
    NT=size(rec.T);    
    NStates=NT(5);
    subT=ones(1,NStates);%To perform subdivisions of T 

    cT=zeros(NT(1:5),'single');%Flag to indicate whether a given motion state can be assumed to have converged     
    w=parXT.winit*ones(NT(1:5),'single');%Weights of LM iteration            
    NSamplesPerState=accumarray(stateSample(:),ones(NProfs,1))';    
    voxSizRes=voxSiz .* ( NX ./  ceil(NX(1:3)./(2.^(-log2(resAni)))) );
    if rec.Dyn.Debug>=2
        fprintf('\n========================== Resolution level: %d ==========================\n',l);
        fprintf('Voxel size: %.2fx%.2fx%.2f mm\n',voxSizRes );
        fprintf('Number of motion states: %d\n',NStates);
        fprintf('Minimum/Maximum number of echoes: %d/%d\n',min(NSamplesPerState),max(NSamplesPerState));
        if any(NSamplesPerState(:)<parXT.echoesMin);fprintf('Trying to operate below the minimum required number of reads, motion estimates will not be obtained for stability\n');end
    end
    if estT(l)
        cT(NSamplesPerState(:)<parXT.echoesMin)=1;
        cT=reshape(cT,NT(1:5));
    end

    %ENCODING STRUCTURES AND BLOCK SIZES FOR GPU COMPUTATION AT THIS LEVEL
    if prod(round(NX(1:3).*resAni))<=rec.Dyn.MaxMem(2) && gpuIn;S=gpuArray(S);rec.M=gpuArray(rec.M);end
    if ~discMa%Soft masking
        [xRes,SRes,yRes,timeIndexRes,R.Ti.la,E.rG,E.kG,E.rkG,E.Fof,E.Fob]=downsampleOperators(-log2(resAni),x,S,y,timeIndex,Ti.la,gibbsRing,0,gpuIn); 
        R.Ti.la=1./(abs(R.Ti.la).^2+0.01);
    else
        [xRes,SRes,yRes,timeIndexRes,CX.Ma,E.rG,E.kG,E.rkG,E.Fof,E.Fob]=downsampleOperators(-log2(resAni),x,S,y,timeIndex,rec.M,gibbsRing,0,gpuIn);
    end
    
    NXres=size(xRes);NXres(end+1:3)=1;NYres=size(yRes);NYres(end+1:5)=1;
    [MSRes,MTRes] = mapNIIGeom(voxSiz,rec.Par.Mine.APhiRec*Tpermute(rec.Par.Mine.permuteHist{end}) , 'resampling',[],NX,NXres);
    
    %B0 FIELD INITIALISATION
    %%% Taylor model
    if ~parXT.exploreMemory && any(estB) && parXB.useTaylor
        E = updateTaylor (E, xRes, NXres, MSRes, parXB, rec.Par.Labels);
        E.Dl.Geom.APhiRecOrig = rec.Par.Mine.APhiRec;
        E.Dl.Geom.permuteHist = {rec.Par.Mine.permuteHist{end}};
        if isfield(parXB,'blockConversion');E.Dl.Geom.block = parXB.blockConversion;end
    end

    %if l==L;S=[];y=[];end
    for n=1:2;yRes=fftGPU(yRes,n)/sqrt(NYres(n));end 
    %WEIGHTING FUNCTION FOR MOTION ESTIMATION (ROI)
    if strcmp(rec.Par.Labels.FatShiftDir,'F');E.nF=1:floor((1-parXT.redFOV)*NXres(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'H');E.nF=floor(1+parXT.redFOV*NXres(3)):NXres(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.nF=1:NXres(3);end
    EH.Mbe=zeros([1 1 NYres(3)],'like',real(xRes));EH.Mbe(E.nF)=1;%To use the extracted region for energy computation

    %HARMONIC SAMPLING
    E.mSt=timeIndexRes;  
    E.mSt(timeIndexRes~=0)=stateSample(timeIndexRes(timeIndexRes~=0)); 

    if estT(l) || ~isfield(EH,'We');EH.We=[];end    
    if parXT.fractionOrder~=0;GRes=buildFilter([NXres(1:2) length(E.nF)],'FractionalFiniteDiscreteIso',NX./NXres,gpuIn,parXT.fractionOrder);else GRes=[];end%For fractional finite difference motion estimation           
    [ySt,FSt,GRes,E.iSt,E.nSt,E.mSt,nSa]=buildHarmonicSampling(yRes,GRes,E.mSt,NStates,NXres(1:2));iSt{l}=cat(1,E.iSt{:});yRes=gather(yRes);%Samples as cell arrays
    if any(nSa(1:NStates)==0);nowithin=1;fprintf('Probably not DISORDER, returning without performing corrections\n');return;end          
    
    if leff<=1;mStSweeps=E.mSt;end
    BlSz=max(sqrt(prod(voxSiz(1:2)./resAni(1:2))),1);%Block sizes baseline
    E.bS=round(rec.Dyn.GPUbS.^log2(BlSz));E.oS=on{2};%Block sizes for our GPU
    if rec.Dyn.Debug>=2;fprintf('Block sizes:%s\n',sprintf(' %d',E.bS));end
    E.Je=1;%To use joint encoding/decoding
    E.Uf=cell(1,3);for n=1:3;E.Uf{n}.NY=NYres(n);E.Uf{n}.NX=NXres(n);end%Folding
    EH.Ub=cell(1,3);for n=1:3;EH.Ub{n}.NY=NYres(n);EH.Ub{n}.NX=NXres(n);end%Unfolding
    [E.UAf,EH.UAb]=buildFoldM(NXres(1:3),NYres(1:3),gpuIn,1);%Folding/Unfolding matrix (empty is quicker!)

    %CONTROL OF TRANSFORM CONVERGENCE / CONSTRAIN THE TRANSFORM / ROI IN
    %READOUT FOR MOTION ESTIMATION
    E.w=w;E.cT=cT;E.winit=parXT.winit;
    if parXT.meanT;CT.mV=5;else CT=[];end
    E.tL=parXT.traLimXT*resIso(l);%Minimum update to keep on doing motion estimation

    %% Actual data synthesis
    
    rng(22)
    EGT=E;
    
    %%% Motion trace generation
    if isfield(rec.Sim , 'expFact'); interleaveInfo = [rec.Sim.NinterleavesT, rec.Sim.expFact] ; else interleaveInfo = [rec.Sim.Ninterleaves 0];end
    TGT = generateTransformTrace(NT(5),rec.Sim.motionType, interleaveInfo, rec.Sim.tran, rec.Sim.rot, rec.Sim.meanT, timeState, voxSiz );
    if ~rec.Sim.synthT; TGT = 0 * TGT; end

    DGT = permute(rec.Sim.DGT,[3 2 1 5 6 4]); DGT = single(DGT); if gpu;DGT= gpuArray(DGT);end 
    xGT = permute(rec.Sim.xGT,[3 2 1 5 6 4]); xGT = single(xGT); if gpu;xGT= gpuArray(xGT);end 
    
    %%% Taylor maps generation
    if rec.Sim.synthD       
        DRes = resampling(DGT,[NXres, on{2}, 2],0,2*ones(1,6));
        EGT.Dl.D = DRes;
    end
    
    %%% Image generation
    xGTRes = resampling(xGT,size(xRes));
    if ~rec.Sim.synthX; xGTRes = 0 * xGTRes;end  
  
    [EGT.Tr,EGT.Fs,EGT.nEc]=compressMotion(TGT,FSt, resIso(end), parXT, [],rec.Sim.disableGrouping);
    if ~isempty(EGT.Fs);EGT.NSe=length(EGT.Fs{1});EGT.NMs=size(EGT.Tr,5);else EGT.NSe=1;EGT.NMs=1;end%Number of segments (including out of elliptic shutter) / Number of motion states            
    EGT.NRe=EGT.NSe-EGT.NMs; %Number of outern segments
    
    EGT = transformationConversion(EGT, [], [], 'Dl');%Make motion parameters for Taylor model 

    EGT.dS(1)=EGT.NSe;
    ncx=1:eivaS(1+2*estT(l));
    EGT.Sf=dynInd(SRes,ncx,4);EGT.dS(2)=ncx(end);
    if gpuIn;EGT.Sf=gpuArray(EGT.Sf);end
    
    %Encode to generate k-space
    ySt = encode(xGTRes,EGT);

    %%% Add noise
    sigma=SNRToLevels(rec.Sim.snrdB,xGTRes,SRes);
    for n=3; ySt=fftGPU(ySt,n)/sqrt(NYres(n));end%noise has to be added in full k-space, so still need to make 3rd dimension in k-space
    no = sigma*plugNoise(ySt)/sqrt(prod(NYres(1:3)) ); %Need to rescale since using normalised DFT
    ySt = ySt + no;
    for n=3; ySt=ifftGPU(ySt,n)*sqrt(NYres(n));end
    
    %%% Reshape back
    ySt = dynInd(ySt,iSt{1},1 , dynInd(ySt,1:length(iSt{1}),1) ); 
    ySt = reshape(ySt, NYres([ 1:2 5 3:4])); 
    ySt = permute(ySt,[1 2  4 5 3]);

    %%% Back into image domain in first 2 dimensions
    for n=1:2;ySt=ifftGPU(ySt,n)*sqrt(NYres(n));end
    
    %%%Return
    y=ySt;
    disp('=======================  END DATA SYNTHESIS  ===============================');

    %%%make TGT on temporal resolution of each repetition time
    TGTFullTempRes = zeros( [ on{4}, length(stateSample),size(TGT,6)] );
    for i=1:max(stateSample) 
        idx = find(stateSample==i);
        TGTFullTempRes = dynInd(TGTFullTempRes, idx, 5, repmat( dynInd( TGT,i,5),[on{4} length(idx) 1 ]) );
    end
    
    TGTFullTempResSyntRes = zeros( [ on{4}, length(find(EGT.mSt<EGT.NMs)) ,size(TGT,6)] );
    for i=1:max(EGT.NMs)
        idx = find(EGT.mSt==i);
        TGTFullTempResSyntRes = dynInd(TGTFullTempResSyntRes, idx, 5, repmat( dynInd( TGT,i,5),[on{4} length(idx) 1 ]) );
    end

    return
    
end
end