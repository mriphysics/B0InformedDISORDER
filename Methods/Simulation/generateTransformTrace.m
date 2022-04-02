
function [T] = generateTransformTrace(N, type, NinterInfo, tran, rot, meanT, time, MS )

%GENERATETRANSFORMTRACE generates a transformation trace to be used in simulations.
%   [T] = GENERATETRANSFORMTRACE({N},{TYPE},{NINTERINFO},{TRAN},{ROT},{MEANT},TIME,MS) 
%   * {N} the number of motion states
%   * {TYPE} the type of motion
%   * {NINTERINFO} info about the interleaved pattern of motion
%   * {TRAN} the standard deviation of the translation in mm
%   * {ROT} the standard deviation of the rotation in degrees
%   * {MEANT} whether to enforce a motion trace with zero mean
%   * TIME a time array of the acquisition to determine the number of samples per interleave
%   * {MS} the resolution
%   ** T the generated motion trace
%
%   Yannick Brackenier 2022-01-30

if nargin < 1 || isempty(N); N=[]; end 
if nargin < 2 || isempty(type); type  = 'interleaved';end
if nargin < 3 || isempty(NinterInfo) ; NinterInfo = [N 0];end%Second element is the factor to be used for the exponential decay of each segment.
if nargin < 6 || isempty(meanT); meanT = 0;end
if nargin < 8 || isempty(MS); MS = ones(1,3);end
if length(NinterInfo)==1; NinterInfo = [NinterInfo 0];end
time = multDimMea(time,2:ndims(time));
if isempty(N); N = length(time);end
on = oneStruct();

NT=ones(1,6);
NT(5)=N;%motion states in the 5th dimension
NT(6)=6;%3D motion parameters in the 6th dimension
Ninter = max(1,round(time(end) / NinterInfo(1)));
tran = tran./MS; tran = permute(tran(:), [2:6 1]);

if isequal(type,'interleaved') || isequal(type,'interleavedLinear') || isequal(type,'interleavedExponential')
    %%% Make indices for different interleaves
    sameplePerInterval = ceil(N/Ninter);
    idxState = ceil(Ninter*(((1:N)-0.5)/N));

    %%% Make T trace
    T = single(zeros(NT)); %In accordance with compressMotion.m
    for i = 1:Ninter
        idx = find(idxState==i);
        T = dynInd(T, {idx, 1:3}, 5:6,          tran.*( repmat(rand([on{5},3])-0.5, [on{4} ,length(idx)])    ));
        T = dynInd(T, {idx, 4:6}, 5:6,  pi/180* rot*(   repmat(rand([on{5},3])-0.5, [on{4} ,length(idx)])    ));
    end
    %%% Linaer
    if isequal(type,'interleavedLinear')
        Tlin = T;
        for i = 1:Ninter
            idx = find(idxState==i);
            if i~=Ninter; idx = [idx idx(end)+1];end
            if i==1; Tlow=0*dynInd(T,idx(1),5); else; Tlow = dynInd( T,idx(1)-1,5);end
            Tup = dynInd(T,idx(1),5);
            Ttemp = bsxfun(@times, dynInd(T,[ idx(1:end-1) , idx(end-1)],5),permute(linspace(1,2,length(idx))',[2 3 4 5 1]) );
            Ttemp =  rescaleND( Ttemp, cat(2,Tlow,Tup), cat(2,dynInd(Ttemp,1,5), dynInd(Ttemp,length(idx),5) )) ;
            Tlin = dynInd(Tlin,idx,5,Ttemp);
        end
        T=Tlin;
    end
    if isequal(type,'interleavedExponential')
        TinterExp = T;
        for i = 1:Ninter
            idx = find(idxState==i);
            idxExponential = linspace( -1, 1, length(idx) ); %(idx - idx( ceil((length(idx)+1)/2) ) )/ (0.5*length(idx)) ;%ranges from -1 to 1
            factor = NinterInfo(2); 
            b = (log(1-factor) + log(1+factor) )/2;a = log(1-factor)-b;
            modulator = exp(a*idxExponential +b).';
            Ttemp = bsxfun(@times, dynInd(T,idx,5), permute(modulator,[2 3 4 5 1]) );%[ idx(1:end-1) , idx(end-1)] since you don't want the next T to be inclluded in the range
            TinterExp = dynInd(TinterExp,idx,5,Ttemp);
        end
        T=TinterExp;
    end
else
    error('generateTransformTrace:: motion type not implemented.');
end

%%% Remove mean
if meanT
    T = T - multDimMea(T,5);
end
% T = dynInd(T, 4:6, 6 ,dynInd(T,[6 5 4],6));

end