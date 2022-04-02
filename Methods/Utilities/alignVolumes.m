
function [xo,T] = alignVolumes(xi, res, orient, ref, exclDim, excludeRef, padding, debug)

%ALIGNVOLUMES  Aligns a set of volumes x to the reference volume x{ref} taking into account a differences in resolution and/or field of view (FOV). 
%   [X,T]=ALIGNVOLUMES(X,{RES},{REF},{EXCLDIM},{EXCLUDEREF},{PADDING},{DEBUG})
%   * X is a cell array with the different volumes to align (w.r.t. x{ref}). Volumes to apply the transformation to (e.g. in fMRI) are stored in the 4th dimension of x.
%   * {RES} is a cell array with the respective resolutions.
%   * {ORIENT} is a cell array with the respective sform matrices.
%   * {REF} the reference volume to align to. Defaults to the first.
%   * {EXCLDIM} the dimension to use for FOV extraction during groupwise volume registration
%   * {EXCLUDEREF} to remove the reference volume from the output to save memory. Defaults to 0.
%   * {PADDING} padding in mm to add to the image to avoid edges of FOV wrapping (as implemented in sincRigidTransform.m).
%   * {DEBUG} a flag to show intermediate results and the region to be excluded during registration (defaults to 0). If 1, show area excluded. If 2 show all intermediate results.
%   ** X are the aligned volumes, arranged in a multi-dimensional array with the same FOV and resolution as x{ref}. Volume x{2,:} are stacked in the 5th dimension
%   ** T cell of motion parameters to tranfsorm from x{vol} to x{ref} in the common image space (see below).
%
%   Yannick Brackenier 2022-01-30

applyToOther = size(xi,1)>1;

if nargin < 2 ; res = []; end 
if nargin < 3 ; orient = []; end 
if nargin < 4 || isempty(ref); ref = 1; end 
if nargin < 5 || isempty(exclDim); exclDim=[];end %Region exclusion will not be used
if nargin < 6 || isempty(excludeRef); excludeRef=0;end
if nargin < 7 || isempty(padding); padding=zeros(1,3);end
if nargin < 8 || isempty(debug); debug=0;end

%%% Make resolutions compatible
if isempty(res); res = cell(1,size(xi,2));end
if length(res)~=size(xi,2); res = dynInd(res, size(res,2)+1:size(xi,2),2,{[]});end
if isempty(res{ref}); res{ref}=ones(1,3);end
idx = cellfun(@isempty, res);
if any(idx); res(idx)= { res{ref} }; end %replace empty resolution with the resolution from the 1st volume

%%% Make sure orient is compatible
if isempty(orient); orient = cell(1,size(xi,2));end
if length(orient)~=size(xi,2); orient = dynInd(orient, size(orient,2)+1:size(xi,2),2,{[]});end

%%% Check whether transformation to be applied to other volumes
if applyToOther
    xToTransform = dynInd(xi,2,1); 
    xi = dynInd(xi,1,1); %take out the volumes from x that are used for registration
    %%% Fill empty cell arrays with the volumes in x
    idx = cellfun(@isempty, xToTransform);
    xToTransform = dynInd( xToTransform, find(idx),2, dynInd(xi,find(idx),2) );  
    %%%Check if number of arrays to transform are same for every pose
    NtoTrans = cellfun(@(x) size(x,4), xToTransform, 'UniformOutput', false);
    NtoTransMax = max(cat(1,NtoTrans{:})); 
    for vol =1:size(xToTransform,2)
        NtoAdd = NtoTransMax - NtoTrans{vol};
        xToTransform{vol} = cat(4,xToTransform{vol},repmat( dynInd( xToTransform{vol}, 1,4) , [1 1 1 NtoAdd]) );
    end
end

%%% Permute so that arrays in correct dimensions
for vol =1:size(xi,2)
    if ~isempty(orient{vol})
        [permMT{vol},flMT{vol}] = T2perm(orient{vol}); permMT{vol}(end+1:12) = length(permMT{vol})+1:12;
        xi{vol} = flipPermute(xi{vol},flMT{vol},permMT{vol}); res{vol} = res{vol}(permMT{vol}(1:3));
        if applyToOther;  xToTransform{vol}= flipPermute(xToTransform{vol},flMT{vol},permMT{vol}); end
    end
end
    
%%% Make same resolution
N = cellfun(@(x) size(x), xi, 'UniformOutput', false);
nVolumes = size(N,2);
%assert( isequal( cellfun(@length, dynInd(N,1,1)), 3*ones(1,nVolumes) ) ,'alignVolumes:: volumes are expected for all arrays in x.')

for vol = 1:nVolumes 
    if ~isequal(res{ref},res{vol})
        Nv=round(N{vol}.*res{vol}./res{ref});
        xi{vol}=resampling(xi{vol},Nv);%x{vol} has same resolution as x{1}
        if applyToOther;  xToTransform{vol} = resampling(xToTransform{vol},[Nv NtoTransMax]);end %assumed that xToTransform{vol} has same resolution as x{vol}
    end     
end
Nnew = cellfun(@size, xi, 'UniformOutput', false);

%%% Bring to common space
Ncom = multDimMax( cat(4,Nnew{:}), 4 );
Ncom = Ncom + round(2*(padding./res{ref}));%factor 2 since at both sides
centCom = ceil((Ncom+1)/2);%centre of common space

for vol = 1:nVolumes
    grid = generateGrid(Nnew{vol}, [],Nnew{vol},ceil((Nnew{vol}+1)/2));
    ind = {centCom(1) + grid{1}, centCom(2) + permute(grid{2},[2 1 3]),centCom(3) + permute(grid{3},[3 1 2]) };
    if vol==ref; indRef = ind;end %save to later extract FOV of first volume

    xi{vol} = dynInd(zeros(Ncom,'like', xi{1}), ind,1:3,xi{vol});%Used to be called xcom (common space) but this way memory is saves
    if applyToOther; xToTransform{vol} = dynInd(zeros([Ncom NtoTransMax ],'like', xi{1}), ind,1:3,xToTransform{vol}); end
end

%%% Registration settings
dimStore = 5;%Fifth dimension so applicable to 4D arrays of xToTransform
pyr=[8 4 2 1];
tolAccel=0.01;
fraction_order = 0.5;
w = cat(dimStore,xi{:});
if applyToOther; wToTransform = cat(dimStore,xToTransform{:});end
    
%%% Region exclusion
if length(exclDim)>1; maskFactor = exclDim(2);else maskFactor = .3;end
NW=size(w); W= ones(NW(1:3),'like',real(w)); %all ones
if ~isempty(exclDim)
    if exclDim(1) < 0 %Inverse direction of the exclDim
        exclDim(1)  = abs(exclDim(1));
        W = dynInd(W, floor(NW(exclDim(1))*(1-maskFactor):NW(exclDim(1))), exclDim(1),0); 
    else%Take absolute but take different direction
        W = dynInd(W, 1:floor(NW(exclDim(1))*maskFactor), exclDim(1),0); 
    end
    if debug>0; plotND([],abs(xi{1}),[],[],0,{[],2},[],[],W,{1}); title('alignVol.m:: Region excluded from rigid registration for ref volume.');end
end            

%%% Register (registration based on magnitude images)
if debug>1; plotND([],abs(w),[],[],0,{[],1},[],[],[],[],200); sgtitle('alignVol.m:: Before alignment');end
[~,Treg,~,w]=groupwiseVolumeRegistration(abs(w),W,[],0,pyr,fraction_order,[],w,tolAccel); % export Treg for motion estimate
W=[];

%%%% Transform back to reference volume x{ref}
[~,kGrid,rkGrid,~,~] = generateTransformGrids(size(dynInd(w,1,dimStore)));
[etRef] = precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(Treg,ref,dimStore),1);     
w = sincRigidTransform(w, etRef, 1, [],[],0); 

if applyToOther
    for vol=1:nVolumes
        di = [0 1];
        perm = [1:3 4:6];
        Tpair = permute(dynInd(Treg,[vol,ref],dimStore), perm); %compositeTransform.m expects states and parameters in 5th and 6th dimension
        Ttemp{vol} = compositeTransform(gather(Tpair) , di);
    end
    Ttemp = cat(5, Ttemp{:}); Ttemp = ipermute(Ttemp, perm);
    [etToTr] = precomputeFactorsSincRigidTransform(kGrid,rkGrid,Ttemp, 1);    
    wToTransform = sincRigidTransform(wToTransform, etToTr, 1); %Implicetely also applies this to other dimensions (4th in this case)
end
if debug>1; plotND([],abs(w),[],[],0,{[],1},[],[],[],[],201); sgtitle('alignVol.m:: After alignment');end

%%% Stack together
if applyToOther; xo = cat(4,w, wToTransform);wToTransform=[]; else xo = w; end;  w=[];%Different poses in 5th dimension

%%% Extract original FOV x{ref}
xo = dynInd(xo,indRef,1:3);

%%% Exclude reference volume
if excludeRef; xo = dynInd(xo, setdiff(1:size(xi,4),ref) ,dimStore); end

%%% Flip/permute first 3 dimensions back to have orientation as the original reference volume
if ~isempty(orient{ref}); xo = flipPermute(xo,flMT{ref},permMT{ref},0); end
        
%%% Output transformation parameters
if nargout>1
    T = Ttemp;
    if excludeRef; T = dynInd(T, setdiff(1:size(T,4),ref) ,4); end
    warning('alignVolumes:: Motion parameters returned might have offset due to difference in FOV.');
end

end
            