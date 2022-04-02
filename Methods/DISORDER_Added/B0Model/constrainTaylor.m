
function [E, Dcon] = constrainTaylor(E,D)

%CONSTRAINTAYLOR  Constrains the Taylor maps of the B0 model. 
%   [E,DCON]=CONSTRAINTAYLOR(E,{D})
%   * E is the encoding structure, which can contain the Taylor maps in E.Dl.D.
%   * {D} are Taylor maps to be contrained. If provided, the maps in E.Dl will not be used.
%   ** E the encoding structure, with optionally the original maps in E.Dl.D
%   ** DCON contrained Taylor maps
%    

if nargin<2; Dcon = E.Dl.D; else; Dcon=D; nd 

% Filtering
if isfield(E.Dl, 'C') && isfield(E.Dl.C, 'H') && ~isempty(E.Dl.C.H)
    fprintf('Filtering Taylor terms.\n');
    for i =1:size(Dcon,6); Dcon = real( dynInd(Dcon,i,6, filtering(dynInd(Dcon,i,6), E.Dl.C.H, 1)));end %Filter in cosine-domain
end
   

end