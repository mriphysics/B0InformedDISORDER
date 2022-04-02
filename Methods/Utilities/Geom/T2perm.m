

function [perm, fl] = T2perm(MT)

%T2PERM finds a flip and permution order that (when applied) approximates a given transformation matrix.
%   [PERM,FL]=T2PERM(MT)
%   * MT the transformation matrix.
%   ** PERM the permutation order.
%   ** FL the flip indices for every dimension.
%
%   Yannick Brackenier 2022-01-30


MS=sqrt(sum(MT(1:3,1:3).^2,1));
MT = MT(1:3,1:3);

R = MT/diag(MS);
perm = 1:size(R,1);
fl = zeros(1,size(R,1));

[~,idx] = max(abs(R),[],1);
[~,idxMin] = min((R),[],1);

perm(idx) = perm;
fl(idx==idxMin) = 1;
