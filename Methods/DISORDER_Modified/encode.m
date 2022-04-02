
function xou=encode(x,E)

%ENCODE   Applies a given encoding to the data
%   XOU=ENCODE(X,E)
%   * X is the data before encoding
%   * E is the encoding structure
%   ** XOU is the data after encoding
%

if isempty(E);xou=x;return;end
inBlSz=1e3;

if ~isfield(E,'bS')%Block sizes by default
    E.bS=[size(x,5) size(x,4)];E.dS=[size(x,5) size(x,4)];E.oS=[1 1]; 
    if isfield(E,'Sf');E.bS(2)=size(E.Sf,4);E.dS(2)=size(E.Sf,4);end 
    if isfield(E,'Tr');E.bS(1)=E.NSe;E.dS(1)=E.NSe;end
end
if isfield(E,'Zf') && isempty(E.Zf) && isfield(E,'Sf');E.oS(2)=1;E.bS(2)=size(E.Sf,4);E.dS(2)=size(E.Sf,4);end

if isfield(E,'Af');x=E.Af*x;end%Explicit matrix

if isfield(E,'Ti');xTi=E.Ti.la*x;end%This does not have counterpart in decoding
if isfield(E,'Gf');x=filtering(x,E.Gf,E.mi);end%Filtering
if isfield(E,'Xf');x=bsxfun(@times,x,E.Xf);end%Body coil or mask
if isfield(E,'Ti');x=x+xTi;end%This does not have counterpart in decoding
xou=x;
if isfield(E,'Ff')%Filtering
    xou=[];
    if isfield(E.Ff,'Fd')
        for n=1:E.Ff.ND
            padEl=zeros(1,E.Ff.ND+1);padEl(n)=1;
            xs=padarray(diff(x,1,n)/E.Ff.w(n),padEl,0,'post');
            if n==1;xou=xs;else xou=cat(E.Ff.ND+1,xou,xs);end  
        end
    else    
        ND=numDims(x);    
        N=size(x);
        for n=1:length(E.Ff)
            mirr=zeros(1,ND);mirr(n)=E.ma(n);
            xs=mirroring(x,mirr,1);
            xs=filtering(xs,E.Ff{n},E.mi(n));
            xs=resPop(xs,n,[N(n) 1+E.ma(n)],[n ND+1]);
            if E.Sm;xs=gather(xs);end
            if n==1;xou=xs;else xou=cat(ND+1,xou,xs);end
        end
    end
    x=xou;
end

if isfield(E,'NXAcq');x=resampling(x,E.NXAcq);end    

%SLAB EXTRACTION
if isfield(E,'ZSl') && E.ZSl>0;x=extractSlabs(x,abs(E.ZSl),1,1);end
    
for a=E.oS(1):E.bS(1):E.dS(1);vA=a:min(a+E.bS(1)-1,E.dS(1));%YB: Transformations
    if isfield(E,'vA');vA=vA(ismember(vA,E.vA));end %YB: useful in the LMsolver.m
    xT=x;

    %RIGID TRANSFORM (MOTION STATES)
    if isfield(E,'Tr') && ~isempty(E.Tr)
        if any(vA<=E.NMs)

            if any(E.Tr(:)~=0)              
                Tr=dynInd(E.Tr,vA(vA<=E.NMs),5);
                %DEPHASING - motion dependent
                if isfield(E,'Dl')
                    TrD = dynInd(E.Dl.Tr,vA(vA<=E.NMs),5);
                    xT=bsxfun(@times,xT,dephaseRotationTaylor(TrD,E.Dl.f, E.Dl.D,E.Dl.TE));
                end 
                %TRANSFORM
                if isfield(E,'Tf');Tf=extractFactorsSincRigidTransform(E.Tf,vA(vA<=E.NMs),5);
                else Tf=precomputeFactorsSincRigidTransform(E.kG,E.rkG,Tr,1,0,1,1);
                end
                xT=sincRigidTransform(xT,Tf,1,E.Fof,E.Fob);                       
            end
            if size(xT,5)~=length(vA(vA<=E.NMs)); xT=repmat(xT,[ones(1,4) length(vA(vA<=E.NMs))]);end
            
        end
            
        if any(vA>E.NMs)
            assert(~isempty(xT),'Empty transformed array');
            xT=cat(5,xT,repmat(x,[ones(1,4) sum(vA>E.NMs)]));
        end
    end
            
    %FILTERING (GENERALLY FOR SLICE RECOVERY)
    if isfield(E,'Sp');xT=filtering(xT,E.Sp);end
    
    %SLAB EXTRACTION
    if isfield(E,'ZSl');xT=extractSlabs(xT,abs(E.ZSl),0,1);end
        
    xouT=xT;    
    %COIL PROFILES
    for b=E.oS(2):E.bS(2):E.dS(2);vB=b:min(b+E.bS(2)-1,E.dS(2));%YB: coil sensitivities
        xS=xT;
        if isfield(E,'Sf')
             Saux=dynInd(E.Sf,vB,4);
            if isa(xS,'gpuArray');Saux=gpuArray(Saux);end
            xS=bsxfun(@times,xS,Saux);
            if size(Saux,6)>1;xS=sum(xS,6);end
        end%Sensitivities         
        
        %SENSE
        if isfield(E,'Uf')%Sense folding (first two dimensions)
            for n=1:2
                if ~isempty(E.Uf{n})
                    if E.Uf{n}.NX>E.Uf{n}.NY;xS=fold(xS,n,E.Uf{n}.NX,E.Uf{n}.NY,E.UAf{n});end
                    if E.Uf{n}.NX<E.Uf{n}.NY && (~isfield(E,'Fs') || isempty(E.Fs))
                        NXS=size(xS);NXS(n)=E.Uf{n}.NY;
                        xS=resampling(xS,NXS,2);
                    end     
                end
            end
        end
        
        %MB/NYQUIST
        if isfield(E,'Ef');xS=aplGPU(E.Ef,xS,E.pe);end%Fourier
        if isfield(E,'Bf');xS=bsxfun(@times,xS,dynInd(E.Bf,E.cc,6));end%MB encoding
        if isfield(E,'Gh');xS=bsxfun(@times,xS,E.Gh.Af);end%Nyquist ghosting encoding
        
        if isfield(E,'Uf')%Sense folding (third dimension)    
            for n=3
                if isfield(E,'Bf');xS=resSub(sum(resSub(xS,3,[E.Uf{n}.NY E.Uf{n}.NX/E.Uf{n}.NY]),4),3:4);elseif ~isempty(E.Uf{n});xS=fold(xS,n,E.Uf{n}.NX,E.Uf{n}.NY,E.UAf{n});end  
            end
        end        
        if isfield(E,'Es') && E.Es==1 && isfield(E,'Ef');xS=aplGPU(E.Ef',xS,E.pe);end%Space
                
        xouS=xS;
        %FOURIER DOMAIN (SEGMENTS)
        if isfield(E,'Fs') && ~isempty(E.Fs) 
            for c=1:length(vA)
                xR=dynInd(xS,c,5);xouR=[]; % YB: At this point one motion state extracted
                if ~isempty(E.Fs{1}{vA(c)})
                    Nec=size(E.Fs{1}{vA(c)},1);
                    for d=1:inBlSz:Nec;vD=d:min(d+inBlSz-1,Nec);
                        xQ=xR;        
                        xQ=aplGPU(dynInd(E.Fs{1}{vA(c)},vD,1),xQ,1);%YB: This is Fmatrix for first dimension of 1 motion state, and even subextracted for indiced vD (often inBlsz >> number of samples in state) 
                        if size(E.Fs,1)==2
                            xQ=sum(bsxfun(@times,xQ,dynInd(E.Fs{2}{vA(c)},vD,1)),2);
                        end
                        if d==1;xouR=xQ;else xouR=cat(1,xouR,xQ);end %YB: store all the samples in the first dimension, as ySt is also resized (see HarmonicSampling.m and decode.m)                            
                    end;xQ=[];
                end
                if c==1 || isempty(xouS);xouS=xouR;elseif ~isempty(xouR);xouS=cat(1,xouS,xouR);end
            end;xR=[];xouR=[];
        end     
        
        %FOURIER DOMAIN (MULTISLICE)
        if isfield(E,'Fms') && ~isempty(E.Fms);xouS=matfun(@mtimes,dynInd(E.Fms,vA,5),xouS);end
        
        %SLICE MASK
        if isfield(E,'Bm') && ~isempty(E.Bm);xouS=sum(bsxfun(@times,xouS,E.Bm),6);end
        
        if b==E.oS(2) || isempty(xouT);xouT=xouS;elseif ~isempty(xouS);xouT=cat(4,xouT,xouS);end %YB: Coils are stacked in 4th dimension                    
    end;xS=[];xouS=[];
    if a==E.oS(1) || isempty(xou);xou=xouT;elseif (~isempty(xouT) && isfield(E,'Fs'));xou=cat(1,xou,xouT);elseif ~isempty(xouT);xou=cat(5,xou,xouT);end%YB: Motion states are stacked in 5th dimension
end;xT=[];xouT=[];
