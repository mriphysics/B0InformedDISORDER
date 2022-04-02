
function rec=solveXTD(rec)

%SOLVEXT  Performs an aligned and B0-corrected reconstruction according to the manuscript
% "Data-driven motion-corrected brain MRI incorporating pose
% dependent B0 fields", Y Brackenier, L Cordero-Grande, R Tomi-Tricot, T
% Wilkinson, P Bridgen, A Price, S J Malik, E De Vita and J V Hajnal, 2022.
%   REC=SOLVEXTB(REC)
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the encoding information (rec.Enc)
%   ** REC is a reconstruction structure with reconstructed data and
%   surrogate information (rec.(rec.Plan.Types))
%

tsta=tic;

%SET DEFAULT PARAMETERS
rec=disorderAlgorithm(rec);
rec=B0Algorithm(rec);

%FILENAME
folderSnapshots=strcat(rec.Names.pathOu, filesep,'An-Ve_Sn');%YB: Anatomical - Volumetric
[~,fileName]=fileparts(generateNIIFileName(rec));

%LOGGING INFO
if rec.Dyn.Log
    logDir = strcat(rec.Names.pathOu,filesep,'An-Ve_Log', filesep);
    logName =  strcat( logDir, fileName,rec.Plan.Suff,rec.Plan.SuffOu, '.txt'); 
    if exist(logName,'file'); delete(logName) ;end; if ~exist(logDir,'dir'); mkdir(logDir);end
    diary(logName) %eval( sprintf('diary %s ', logName))
end
fprintf('Running file: %s \n', strcat( fileName, rec.Plan.Suff, rec.Plan.SuffOu) )
c=clock; fprintf('Date of reconstruction: %d/%d/%d \n \n', c(1:3));

%SHORTCUTS FOR COMMONLY USED STRUCTURE FIELDS AND INITIALIZERS
voxSiz=rec.Enc.AcqVoxelSize;%Acquired voxel size
parXT=rec.Alg.parXT;parXB=rec.Alg.parXB;
gpuIn=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);if gpuIn;gpuFIn=2;else gpuFIn=0;end 
gpu=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);%YB: set gpu to single(gpuDeviceCount && ~rec.Dyn.BlockGPU)

if size(rec.y,5)>=parXT.maximumDynamics;fprintf('Problem too big to fit in memory (%d dynamics for a limit of %d)\n',size(rec.y,5),parXT.maximumDynamics);rec.Fail=1;return;end
on=cell(1,5);for n=1:5;on{n}=ones(1,n);end
typ2Rec=rec.Dyn.Typ2Rec;

%RECONSTRUCTION PLAN
[resPyr,~,~,resIso,~]=pyramidPlan(voxSiz,parXT.resolMax,parXT.NWend,parXT.accel);% YB: resPy is used for resAni, which in its turn is used in downSampingOperators.m
resPyrMax = [ .25 .5  1  ];
resPyr=rec.Alg.resPyr;%resPyr = resPyrMax(end+1-min(length(rec.Alg.resPyr),length(resPyr)):end) ;
resIso=sqrt(2)*((prod(voxSiz).^(1/3))./resPyr);%resIso = [resIso(end-(length(resPyr)-1):end)]
L = length(resIso);
estT = ones(size(resIso));
estB = ones(size(resIso));

fprintf('Reconstruction plan:\n')
fprintf('   Levels:             %s \n',sprintf(' %d ',1:L))
fprintf('   Effective levels:   %s \n',sprintf(' %d ',max(0, (1:L) - sum(resPyr~=1)) )  )
fprintf('   Resolution:         %s \n',sprintf(' %g ',resPyr))
fprintf('   Motion estimation:  %s \n',sprintf(' %d ',estT))
fprintf('   B0 estimation:      %s \n',sprintf(' %d ',estB))

%ROI COMPUTATION AND EXTRACTION IN THE READOUT DIRECTION
rec.Enc.ROI=computeROI(rec.M,[],[0 0 0],[1 0 0]); %YB: At this point non-logical values
if rec.Dyn.Debug>=2
    fprintf('ROI:\n%s',sprintf(' %d %d %d %d %d %d\n',rec.Enc.ROI'));
    fprintf('Acquisition directions: %s\n',rec.Par.Scan.MPS);
end
[~,MT] = mapNIIGeom([],rec.Par.Mine.APhiRec,'dynInd', {rec.Enc.ROI(1,1):rec.Enc.ROI(1,2),rec.Enc.ROI(2,1):rec.Enc.ROI(2,2),rec.Enc.ROI(3,1):rec.Enc.ROI(3,2)}, size(rec.y));

for n=typ2Rec';datTyp=rec.Plan.Types{n};
    if ~ismember(n,[5 12]);rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,1);end%YB: 1 at the end means 'forrward'. At the end also inverse operation
end

%PERMUTE PHASE ENCODES TO THE FIRST DIMENSIONS
perm=1:rec.Plan.NDims;perm(1:3)=[3 2 1]; %YB: RO-PE1-PE2
for n=typ2Rec'        
    if ~ismember(n,[5 12]);datTyp=rec.Plan.Types{n};%YB: 5 and 12 are noise and reconstruction. rec.x does not exist yet, so don't permute.
        rec.(datTyp)=permute(rec.(datTyp),perm);
    end
end

voxSiz=voxSiz(perm(1:3));parXT.apod=parXT.apod(perm(1:3));
if ~isfield(rec.Par.Mine,'permuteHist'); rec.Par.Mine.permuteHist = [];end
rec.Par.Mine.permuteHist{end+1} = perm(1:4);

MS = voxSiz; 
[~,MT] = mapNIIGeom([],MT,'permute',rec.Par.Mine.permuteHist{end} );%Don't change MS since already permuted in voxSiz
isValidGeom(MS,MT);

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

%APODIZE + REARRANGE DATA + CORRECT FOR INTERLEAVED REPEATS %YB: apodisation means windowing
NY=size(y);NY(end+1:rec.Plan.NDims)=1;
y=bsxfun(@times,y,ifftshift(buildFilter(NY(1:3),'tukey',[10 10 1],gpuIn,parXT.apod))); %Apodize % YB: y in image domain
[y,NY]=resSub(y,5:rec.Plan.NDims);NY(end+1:5)=1;%Different averages along the 5th dimension
y=gather(y);

%ACQUISITION STRUCTURE
NEchos=max(rec.Par.Labels.TFEfactor,rec.Par.Labels.ZReconLength);
NProfs=numel(rec.Assign.z{2});
if NProfs<=1 || NProfs~=numel(rec.Assign.z{3});fprintf('SEPARABLE Y-Z TRAJECTORIES. ALIGNED RECONSTRUCTION IS NOT PERFORMED\n');rec.Fail=1;return;end
if sum(cdfFilt(abs(diff(rec.Assign.z{2}(:))),'med'))<sum(cdfFilt(abs(diff(rec.Assign.z{3}(:))),'med'))
    fprintf('Slow direction is: %s\n',rec.Par.Scan.MPS(4:5));
else; fprintf('Slow direction is: %s\n',rec.Par.Scan.MPS(7:8));
end
%YB: z{2} = PE1 and z{3} = PE2
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
if parXB.useTaylor;typ2RecF=29;else typ2RecF=[];end %29-> Linear coefficients of B0 fields with motion
for n=typ2RecI'
    if ~any(rec.Dyn.Typ2Rec==n);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,n);rec.Dyn.Typ2Wri(n)=1;end
end
typ2Rec=rec.Dyn.Typ2Rec;

%INITIALIZE TRANSFORM, NUMBER OF EXTERN ITERATIONS, TOLERANCES, IDENTIFY 
%STEADY STATE SEQUENCES 
rec.T=zeros([on{4} NSweeps 6],'single');

if isfield(parXT,'corrFact'); E.corrFact = parXT.corrFact;end

subT=ones(1,NStates);%To perform subdivisions of T
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
CX.Ma=rec.M;
if rec.Alg.UseSoftMasking==0 && ~all(ismember(CX.Ma(:),[0 1])) % Need to convert rec.M into discrete mask
    %maxint = max(abs(rec.M(:)));  level = graythresh(abs(rec.M ) / maxint);
    maxint =1;  level = 0.45 * mean(abs(rec.M(:)));
    CX.Ma = single( abs(rec.M)> level*maxint);
    %CX.Ma = imopen(CX.Ma, strel('sphere', 3));
    CX.Ma = imclose(CX.Ma, strel('sphere', 15));
    CX.Ma = imdilate(CX.Ma, strel('sphere', 4));
    rec.M = CX.Ma;
end

if gpuIn;CX.Ma=gpuArray(CX.Ma);end
if rec.Alg.UseSoftMasking==2;CX.Ma(:)=1;end %YB: 2 means not masking
discMa=all(ismember(CX.Ma(:),[0 1])); 
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

yX=mean(dynInd(y,ncx,4),5);
if gpuIn;yX=gpuArray(yX);end
gibbsRing=parXT.UseGiRingi*parXT.GibbsRingi;
if gibbsRing~=0 && gibbsRing<=1;yX=filtering(yX,buildFilter(NY(1:3),'tukeyIso',[],gpuIn,gibbsRing));end
nX=nIt(1);
if parXT.exploreMemory;nX=0;end%To test flow without running the main methods       

if nX~=0
    [xRes,EnX]=CGsolver(yX,E,EH,P,CX,R,x,nX,tolType{1},tol(1));
    rec.x=gather(xRes);
end
if ~isempty(EnX)%Print and store final energy
    EnX=EnX/numel(yX);
    if rec.Dyn.Debug>=2;fprintf('Energy evolution last step:%s\n',sprintf(' %0.6g',sum(EnX,1)));end  
    EnX=[];
end;yX=[];

if gpuIn;y=gpuArray(y);end
Residuals=encode(rec.x,E)-y;
for m=1:2;Residuals=fftGPU(Residuals,m)/sqrt(NY(m));end
Residuals=normm(Residuals,[],3:4);
if rec.Alg.WriteSnapshots
    if strcmp(rec.Par.Labels.FatShiftDir,'F');E.nF=1:floor((1-parXT.redFOV)*NX(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'H');E.nF=floor(1+parXT.redFOV*NX(3)):NX(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.nF=1:NX(3);end
    visSegment(flip(permute(dynInd(rec.x,E.nF,3),[2 1 3]),2),[],2,1,[],[],'Uncorrected',strcat(folderSnapshots,filesep,'Reconstructions'),strcat(fileName,rec.Plan.SuffOu,'_Aq'));
    visResiduals([],[],[],Residuals,kTraj,2,{strcat(folderSnapshots,filesep,'ResidualsStates'),strcat(folderSnapshots,filesep,'ResidualsSpectrum')},strcat(fileName,rec.Plan.SuffOu,sprintf('_l=%d',0)));
end
E.Sf=[];
tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing non-corrected reconstruction: %.3f s\n',tend);end

x=resampling(xRes,NX);

if isfield(rec.Par.Labels ,'H0') && rec.Par.Labels.H0 ==3;MTT=[-1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1];else MTT = eye(4);end
rec.Par.Mine.APhiRec = MTT*rec.Par.Mine.APhiRec;
MT = MTT * MT;

%ITERATING THROUGH THE DIFFERENT LEVELS
for l=1:L    
    %WE RESET MOTION CONVERGENCE INSIDE
    nReset=0;THist=[];
    leff=l-sum(resPyr~=1);

    %UPSAMPLE THE TRANSFORM AND GENERATE MOTION PARAMETERS   
    rec.T=repelem(rec.T,1,1,1,1,subT,1); 
    subTT=repelem(subT,1,subT);
    sampleToSweep=repelem(sampleToSweep,1,subT);
    sweepToSample=cell(1,NSweeps);
    c=0;
    for n=1:length(subT)
        if subT(n)==2
            indState=find(stateSample==(n+c));
            Ni=length(indState);
            indState=indState(floor(Ni/2)+1); 
            stateSample(indState:end)=stateSample(indState:end)+1; 
            sweepToSample{sampleToSweep(n+c)}=[sweepToSample{sampleToSweep(n+c)} n+c:n+c+1];
            c=c+1;
        else
            sweepToSample{sampleToSweep(n+c)}=[sweepToSample{sampleToSweep(n+c)} n+c];
        end
    end 
    timeState=regionprops(stateSample,timeSample,'MeanIntensity');
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
    resAni(indMinRes)=voxSiz(indMinRes)/resPyr(l);
    indNoMinRes=1:3;indNoMinRes(indMinRes)=[];
    resAni(indNoMinRes)=max(resAni(indMinRes),voxSiz(indNoMinRes));  
    resAni=voxSiz./resAni;

    NT=size(rec.T);    
    NStates=NT(5);%YB: NStates can change over resolutions but NSweeps does not
    subT=ones(1,NStates);%To perform subdivisions of T 

    cT=zeros(NT(1:5),'single');%Flag to indicate whether a given motion state can be assumed to have converged     
    w=parXT.winit*ones(NT(1:5),'single');%Weights of LM iteration            
    NSamplesPerState=accumarray(stateSample(:),ones(NProfs,1))';
    voxSizRes = voxSiz .* ( NX ./  ceil(NX(1:3)./(2.^(-log2(resAni)))) );
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
        [xRes,SRes,yRes,timeIndexRes,R.Ti.la,E.rG,E.kG,E.rkG,E.Fof,E.Fob]=downsampleOperators(-log2(resAni),x,S,y,timeIndex,Ti.la,gibbsRing,0,gpuIn); %YB: 0 confirms that y is in image domain
        R.Ti.la=1./(abs(R.Ti.la).^2+0.01);
    else
        [xRes,SRes,yRes,timeIndexRes,CX.Ma,E.rG,E.kG,E.rkG,E.Fof,E.Fob]=downsampleOperators(-log2(resAni),x,S,y,timeIndex,rec.M,gibbsRing,0,gpuIn);%YB:attention: discrete mask provided by CX.Ma not used! 
    end
    
    NXres=size(xRes);NXres(end+1:3)=1;NYres=size(yRes);NYres(end+1:5)=1;
    [MSRes,MTRes] = mapNIIGeom(MS, MT, 'resampling',[], NX, NXres);

    %B0 TAYLOR MODEL INITIALISATION
    if ~parXT.exploreMemory && any(estB) && parXB.useTaylor
        E = updateTaylor (E, xRes, NXres, MSRes, parXB, rec.Par.Labels);
        E.Dl.Geom.APhiRecOrig = rec.Par.Mine.APhiRec;
        E.Dl.Geom.permuteHist = {rec.Par.Mine.permuteHist{end}};
    end
    
    if prod(round(NX(1:3).*resAni))<=rec.Dyn.MaxMem(2) && gpuIn;S=gather(S);rec.M=gather(rec.M);end
    if l==L;S=[];y=[];end
    for n=1:2;yRes=fftGPU(yRes,n)/sqrt(NYres(n));end 
   
    %WEIGHTING FUNCTION FOR MOTION ESTIMATION (ROI)
    if strcmp(rec.Par.Labels.FatShiftDir,'F');E.nF=1:floor((1-parXT.redFOV)*NXres(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'H');E.nF=floor(1+parXT.redFOV*NXres(3)):NXres(3);else fprintf('Readout direction is not FH, performance probably suboptimal\n');E.nF=1:NXres(3);end
    EH.Mbe=zeros([1 1 NYres(3)],'like',real(xRes));EH.Mbe(E.nF)=1;%To use the extracted region for energy computation
    if rec.Dyn.Debug>=2; plotND([], dynInd(resampling(abs(rec.x),NXres), E.nF,3),[],[],0,{[],2},MTRes,{'Extracted region for motion estimation'},[],[],233);end

    %HARMONIC SAMPLING
    E.mSt=timeIndexRes;  
    E.mSt(timeIndexRes~=0)=stateSample(timeIndexRes(timeIndexRes~=0)); % YB: Now mSt a 2D grid with the same states having the value of the state number. Also note that timeIndexRes contains indices of the samples, so biggere than numel(timeIndexRes)
    
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

    cont=0; 
    if isfield(E,'En');E=rmfield(E,'En');end    

    %ITERATIVE SOLVER
    for n=1:rec.Alg.nExtern(l)
        if rec.Dyn.Debug>=2;fprintf('--------------Iteration: %d ---------------\n',n);end
        %RESET CONVERGENCE OF MOTION
        if mod(cont,nReset)==0 || all(E.cT(:)) || parXT.convTransformJoint
            nReset=nReset+1; cont=0;%We increment the number of iterations till next reset
            E.cT(:)=0;%We reset
            E.cT(cT==1)=1;%We fix to 1 the low frequency states in case that for whatever reason cT was set to 1
            if rec.Dyn.Debug>=2;fprintf('Resetting motion convergence\n');end
        end
        if (n==1 || l<L) && any(subTT(:)==2) 
            E.cT(subTT==1)=1;
        end%For within-shot in the first iteration we only estimate those coming from outliers
        if rec.Dyn.Debug>=2;fprintf('Explored motion states: %d of %d\n',NStates-sum(single(E.cT)),NStates);end
        cont=cont+1;

        %TRANSFORM PARAMETERS PER SEGMENT/OUTLIERS/DATA/HARMONICS/SEGMENT INDEXES
        if isfield(parXT,'disableGrouping'); disableGrouping = parXT.disableGrouping;else disableGrouping = 0;end
        [E.Tr,E.Fs,E.nEc]=compressMotion(rec.T,FSt,resIso(end),parXT,[],disableGrouping);
        if rec.Dyn.Debug>=2;fprintf('Number of binned motion states: %d of %d\n',size(E.Tr,5),NT(5));end
        if isfield(E,'Dl');E = transformationConversion(E);end%Make motion parameters for linear model 
        
        %SEGMENTS / MOTION STATES / BLOCK SIZES / TRANSFORM FACTORS 
        if ~isempty(E.Fs);E.NSe=length(E.Fs{1});E.NMs=size(E.Tr,5);else E.NSe=1;E.NMs=1;end%Number of segments (including out of elliptic shutter) / Number of motion states            
        E.NRe=E.NSe-E.NMs; %Number of outern segments
        E.dS(1)=E.NSe;

        %CG SOLVER
        ncx=1:eivaS(1+2*estT(l));
        E.Sf=dynInd(SRes,ncx,4);E.dS(2)=ncx(end);
        if gpuIn;E.Sf=gpuArray(E.Sf);end

        %PRECONDITION
        if ~discMa;P.Se=(normm(E.Sf,[],4)+R.Ti.la).^(-1);%Precondition
        else P.Se=(normm(E.Sf,[],4)+1e-9).^(-1);
        end            

        yX=dynInd(ySt,ncx,4);            
        %LAST ITERATION OUTLIER REJECTION
        EH.We=[]; %YB added
        if ((leff>0 && estT(l)==0)|| parXB.weightedEst) && rec.Alg.AlignedRec>=2 && exist('We','var')%YB: use weihts during all levels + for GDsolverTaylor
            EH.We=ones([size(yX,1) 1],'like',Residuals);
            for s=1:NSweeps;EH.We(mStSweeps==s)=We(s);end %Soft weights
            %visResiduals(We,outlWe,timeSweep,[],[],0, [],[],899);%YB see every iteration    
        end

        if estT(l)
            if l==1 && n==1;nX=nIt(2);
            elseif l==1;nX=max(nX-1,1);
            elseif resPyr(l)~=resPyr(l-1) && n==1;nX=max(nIt(2)-1,1);
            elseif resPyr(l)~=resPyr(l-1);nX=max(nX-1,1);
            elseif resPyr(l)==resPyr(l-1) && n==1;nX=0;
            elseif resPyr(l)==resPyr(l-1) && n==2;nX=max(nIt(2)-2,1);
            elseif resPyr(l)==resPyr(l-1);nX=max(nX-(n-2),1);
            end
        else
            nX=nIt(1);
        end
        if parXT.exploreMemory;nX=0;end%To test flow without running the main methods       
               
        if nX~=0 
            
            if estT(l) && estB(l) && parXB.useTaylor && any(E.Tr(:)~=0) && ( l>1 || n >= parXB.taylorDelay) %%% LINEAR MODEL
                fprintf('** B0 estimation - %d iterations \n', parXB.nGD_Taylor);
                %EnBefore = computeEnergy(yX,xRes,E,R,EH);
                E = GDsolverTaylor(yX,E,EH,P,CX,R,xRes,parXB.nGD_Taylor + 1*single(n==parXB.taylorDelay),tol(estT(l)+1), rec.Dyn.Debug);
                %EnAfter = computeEnergy(yX,xRes,E,R,EH);
                %if rec.Dyn.Debug>=2;fprintf('   Energy before: %0.6g / Energy after: %0.6g\n',sum(EnBefore(:)),sum(EnAfter(:)));end  
                if rec.Dyn.Debug>=2;plotND({abs(xRes),1,0},permute(E.Dl.D,[1:3 6 4 5])*pi/180, [-10 10],[], 1,{1:3,2},MTRes,[],single(abs(xRes)>10000),{3},100,sprintf('Estimated Taylor maps at iteration %d',n));end
            
                %Unwrap LC maps
                if E.Dl.C.unwrapTaylor &&  mod(n,5)==1  ; for ii = 1:2; E.Dl.D = dynInd(E.Dl.D, ii,6, unwrapLinear(dynInd(E.Dl.D,ii,6), abs(xRes),0, E.Dl.TE)) ;end; end
               
            end
            
            fprintf('** Image estimation - %d iterations \n', nX );
            %EnBefore = computeEnergy(yX,xRes,E,R,EH);
            [xRes,EnX]=CGsolver(yX,E,EH,P,CX,R,xRes,nX + 8*single(n==1),tolType{estT(l)+1},tol(estT(l)+1));
            %EnAfter = computeEnergy(yX,xRes,E,R,EH);
            %if rec.Dyn.Debug>=2;fprintf('   Energy before: %0.6g / Energy after: %0.6g\n',sum(EnBefore(:)),sum(EnAfter(:)));end  
            if rec.Dyn.Debug>=2; plotND({abs(xRes),1,0}, angle(xRes),[-pi pi 0 2000000],[],0,{[],2},MTRes,[],[],[],123,sprintf('Estimated image at iteration %d',n));end
        end
    
        if ~isempty(EnX)%Print and store final energy
            EnX=EnX/numel(yX);
            if rec.Dyn.Debug>=2;fprintf('Energy evolution last step:%s\n',sprintf(' %0.6g',sum(EnX,1)));end  
            EnX=[];
        end

        %ENERGY COMPUTATION
        E.Sf=dynInd(SRes,ncx,4);E.dS(2)=ncx(end);
        if gpuIn;E.Sf=gpuArray(E.Sf);end
        yX=dynInd(ySt,ncx,4);
        [Re{l}(:,n),~,fullResid]=computeEnergy(yX,xRes,E,[],EH,[],[],[],[],1);Re{l}(:,n)=Re{l}(:,n)/numel(yX);
        En{l}(:,n)=normm(yX,[],setdiff(1:numDims(yX),1))/numel(yX);
        if rec.Dyn.Debug>=2;fprintf('Residuals: %0.6g\n',sum(Re{l}(:,n)));end        

        %ENERGY PER STATES
        fullResid=gather(dynInd(fullResid,iSt{l},1,fullResid));
        Residuals=Re{l}(:,n);Residuals(iSt{l})=Residuals;
        Residuals=reshape(Residuals,NYres([1:2 5]));
        
        mStSort=mStSweeps;mStSort(iSt{l})=mStSort;
        Sweeps=reshape(mStSort,NYres([1:2 5]));         

        %OUTLIER DETECTION
        if leff>0;percUse=(0.15:0.05:0.85)*100;
        else percUse=(0.3:0.05:0.7)*100;
        end
        WeO=zeros([NSweeps length(percUse)],'like',Residuals);
        for s=1:NSweeps;WeO(s,:)=prctile(log(Residuals(Sweeps==s)),percUse);end
        WeO=permute(WeO,[1 3 2]);           
        Westd=diff(prctile(WeO,parXT.percRobustShot*100,1),1,1)./diff(norminv(parXT.percRobustShot)); 
        Wmea=prctile(WeO,mean(parXT.percRobustShot)*100,1);
        Wmea=Wmea+Westd*sqrt(2)*erfcinv(2*mean(parXT.percRobustShot));%sqrt(2)*erfcinv(2*mean(parXT.percRobustShot))=norminv(mean(parXT.percRobustShot))
        We=bsxfun(@minus, WeO,Wmea);%YB changed for version R2015
        We=bsxfun(@rdivide, We, Westd);%YB changed for version R2015
        We=mean(We,3);         
        We=min((1-normcdf(We))/((1-parXT.enerRobustShot)/NSweeps),1);     
        outlWe=We<1; %YB: outlWe flags which sweep is an outlier

        if rec.Alg.WriteSnapshots && estT(l)==0;visResiduals(We,outlWe,timeSweep,resampling(Residuals,NY(1:2),1),kTraj,2,{strcat(folderSnapshots,filesep,'ResidualsStates'),strcat(folderSnapshots,filesep,'ResidualsSpectrum')},strcat(fileName,rec.Plan.SuffOu,sprintf('_l=%d',l)));end
        %CS-LIKE/ROBUST RECONSTRUCTION
        if leff>0 && estT(l)==0
            subT(:)=1;
            if rec.Alg.AlignedRec<4;outlWeV=find(outlWe);else outlWeV=1:length(outlWe);end%AlignedRec==4 explore all, AlignedRec==3, explore only artifacted
            for o=1:length(outlWeV);subT(sweepToSample{outlWeV(o)})=2;end

            if all(subT(:)==1) || nowithin;L=l;end
        end
        if (leff>0 && estT(l)==0) || (rec.Alg.AlignedRec==2 && leff>=0) 
            EH.We=ones([size(yX,1) 1],'like',Residuals); %THIS COULD BE PROBLEMATIC IN CASE OF SEVERAL REPEATS?
            if ismember(rec.Alg.AlignedRec,[3 5])
                for s=1:NSweeps;EH.We(mStSweeps==s)=1-single(outlWe(s));end 
            else%AlignedRec==4 or AlignedRec==1,AlignedRec ==2: use soft weights 
                for s=1:NSweeps;EH.We(mStSweeps==s)=We(s);end
            end
   
        end
        if leff>0 && (estT(l)==0 || n==rec.Alg.nExtern(l))
            if rec.Alg.WriteSnapshots && max(voxSiz)/min(voxSiz)<1.5;visSegment(flip(permute(dynInd(xRes,E.nF,3),[2 1 3]),2),[],2,1,[],[],sprintf('Corrected ($l=$%d)',l),strcat(folderSnapshots,filesep,'Reconstructions'),strcat(fileName,rec.Plan.SuffOu,sprintf('_Di-l=%d',l)));end
            if parXT.computeCSRecon
                fprintf('** Regularised image estimation - %d iterations \n', nX );
                if nX~=0;xCS=CSsolver(yX,E,EH,P,CX,R,xRes,nX,tolType{estT(l)+1},tol(estT(l)+1),rec.Dyn.Debug);else xCS=xRes;end
                if rec.Alg.WriteSnapshots && max(voxSiz)/min(voxSiz)<1.5;visSegment(flip(permute(dynInd(xCS,E.nF,3),[2 1 3]),2),[],2,1,[],[],sprintf('Corrected regularized ($l=$%d)',l),strcat(folderSnapshots,filesep,'Reconstructions'),strcat(fileName,rec.Plan.SuffOu,sprintf('_Re-l=%d',l)));end
            end
            EH.We=[];%YB: change so you can also have EH.We for CGsolver for robust (not CS) reconstruction
            
        end
        
        % MOTION UPDDATE 
        if estT(l) && ~(n==rec.Alg.nExtern(l))
            fprintf('** Motion estimation\n');
            %TRANSFORM PARAMETERS PER SEGMENT/HARMONICS/BLOCK SIZES
            E.Tr=rec.T;E.Fs=FSt; 
            E = transformationConversion(E);%Make motion parameters for linear model 
            E.NMs=size(E.Tr,5);E.NSe=length(E.Fs{1});
            if parXT.exploreMemory;E.cT(:)=1;end%To test flow without running the main methods       

            %EXTRACT RELEVANT ROI AND COILS AND CALL THE MOTION SOLVER                       
            nct=1:eivaS(2);
            E.Sf=dynInd(SRes,nct,4);
            ySf=dynInd(ySt,nct,4);
            E.dS=[E.NMs eivaS(2)];
            if gpuIn;E.Sf=gpuArray(E.Sf);end
            if parXT.fractionOrder~=0;E.Fd=GRes;end
            [E,xRes]=LMsolver(ySf,E,xRes,CT);
            rec.T=E.Tr;
            if isfield(E,'Dl');E = transformationConversion(E);end%YB: Make sure the motion parameters used in the Taylor model are updated as well.
            THist=cat(1,THist,gather(rec.T));
            visMotion(rec,[],[],0,[],[],[],[],[],456);
            %if isfield(E,'Dl');visMotion(E.Dl.Tr,[],[],0,[],[],[],[],[],457);end
        else
            if l==L && parXT.saveFinal
                recW=rec;recW.Alg.SaveRaw=3;recW.Dyn.Batch=1;
                recW.E=E;recW.EH=EH;recW.P=P;recW.CX=CX;recW.R=R;recW.xRes=xRes;recW.yX=yX;recW.fullResid=fullResid;
                recW.timeState=timeState;recW.We=We;recW.Residuals=Residuals;recW.Sweeps=Sweeps;recW.iSt=iSt{l};
                recW.S=[];recW.y=[];
                writeRaw(recW);recW=[];
            end
            break
        end
        E.Sf=[];

        %CHECK CONVERGENCE
        if all(E.cT(:))%Weak convergence
            %if l==L;estT(l)=0;else break;end%If last level we perform another estimation for x           
            estT(l)=0;            
            if rec.Alg.WriteSnapshots;visMotion(rec,voxSiz,timeState,2,strcat(folderSnapshots,filesep,'Motion'),strcat(fileName,rec.Plan.SuffOu,sprintf('_l=%d',l)));end                
        end
    end

    %RECONSTRUCTION TO BE USED TO INITIALIZE THE NEXT RESOLUTION LEVEL
    x=resampling(xRes,NX);
    if ~parXT.exploreMemory && isfield(E,'Ds') && isfield(E.Ds,'chi') && l~=L ;DC=E.Ds.chi;end %YB: save for next level  
    
    if gpuIn;rec.M=gpuArray(rec.M);end
    if discMa;x=rec.M.*x;end
    
    if leff>0
        rec.d=dynInd(rec.d,leff,4,gather(x));
        if parXT.computeCSRecon
            if discMa;xCS=rec.M.*xCS;end
            rec.r=dynInd(rec.r,leff,4,gather(xCS));
        end
    end
    if ~isempty(En{l})
        EnSort=En{l};EnSort(iSt{l},:)=EnSort;
        ReSort=Re{l};ReSort(iSt{l},:)=ReSort;  
        mStSort=E.mSt;mStSort(iSt{l})=mStSort;
        rec.Par.Mine.DISORDER.Residuals{l}=gather(reshape(ReSort,[NYres([1:2 5]) size(ReSort,2)]));%%%PLAY WITH IST SO THAT THIS CAN BE MADE COMPARABLE AND ALSO LINKED TO MOTIONS
        rec.Par.Mine.DISORDER.Energy{l}=gather(reshape(EnSort,[NYres([1:2 5]) size(EnSort,2)]));%%%PLAY WITH IST SO THAT THIS CAN BE MADE COMPARABLE AND ALSO LINKED TO MOTIONS
        rec.Par.Mine.DISORDER.Motions{l}=THist;
        rec.Par.Mine.DISORDER.States{l}=gather(reshape(mStSort,NYres([1:2 5])));
    end

    if parXT.writeInter && leff>0 || l==L 
        % YB added here before E = [];
        rec.timeState=timeState;
        if parXB.useTaylor && ~parXT.exploreMemory %YB: added  && ~parXT.exploreMemory
                rec.D=gather(E.Dl.D); 
                rec.D = permute(rec.D, perm);
                for tt = 1:2; rec.D =dynInd(rec.D,tt,6,extractROI(dynInd(rec.D,tt,6),rec.Enc.ROI,0,1));end 
                for t=typ2RecF'
                    if ~any(rec.Dyn.Typ2Rec==t);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,t);rec.Dyn.Typ2Wri(t)=1;end
                end
                typ2Rec=rec.Dyn.Typ2Rec;
        end
        rec.E = E;
        E=[];EH=[];

        drec=rec.d;
        if parXT.computeCSRecon;rrec=rec.r;end
        ROIrec=rec.Enc.ROI;          
        for n=typ2RecI';datTyp=rec.Plan.Types{n};
            if gpu;rec.(datTyp)=gpuArray(rec.(datTyp));end
            if n~= 12; rr = 1:leff;  else rr = 1;end %YB added because rec.x only 3D
            rec.(datTyp)=dynInd(rec.(datTyp),rr,4);
            rec.(datTyp)=flip(rec.(datTyp),4);
        end
        if l==L && leff~=1;typ2RecL=setdiff(typ2Rec,[5, 29]);elseif l==L && leff==1;typ2RecL=setdiff(typ2Rec,[5,29]);else typ2RecL=typ2RecI;end%YB changed to account for B0 model
        %FINAL ADJUSTMENTS
        for n=typ2RecL';datTyp=rec.Plan.Types{n};
            rec.(datTyp)=permute(rec.(datTyp),perm);
            rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,0,1);
        end
        for n=typ2RecI';datTyp=rec.Plan.Types{n};
            if rec.Alg.MargosianFilter;rec.(datTyp)=margosianFilter(rec.(datTyp),rec.Enc);end
            rec.(datTyp)=removeOverencoding(rec.(datTyp),rec.Alg.OverDec);%REMOVE OVERDECODING  
        end
        rec.Enc=rmfield(rec.Enc,'ROI');
        rec.Dyn.Typ2Wri(17)=1;%To write the motion estimates to file
        typ2RecI=setdiff(typ2RecI,12);
        writeData(rec);
        if l~=L
            rec.Enc.ROI=ROIrec;            
            rec.d=drec;
            if parXT.computeCSRecon;rec.r=rrec;end
        end
    end
    tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time till level %d: %.3f seconds = %.3fminutes\n',l,tend,tend/60);end
    if l==L
        if rec.Dyn.Debug>=2; diary off;end %YB Stop logging
        break;
    end%Early termination if no outliered states detected
    
end
