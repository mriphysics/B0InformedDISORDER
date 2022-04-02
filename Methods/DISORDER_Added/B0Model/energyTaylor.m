

function [EnFi, EnRe, EnReNames] = energyTaylor(y, x, E, EH, W, sep, deb )

%ENERGYTAYLOR  Computes the energy of a solution of the Taylor model of pose-dependent B0 variations. 
%   [ENFI,ENRE,ENRENAMES]=ENERGYTAYLOR(Y,X,E,EH,{W},{SEP},{DEB})
%   * Y is the measured data
%   * X is the reconstructed data
%   * E is the encoding structure
%   * EH is the decoding structure (used for weighted least squares)
%   * {W} is the weights to use for certain regularisers
%   * {SEP} whether to separate the data fidelity loss and the regularisatoin loss
%   * {DEB} whether to display the energies
%   ** ENFI the energy of the data fidelity part
%   ** ENRE the energy of the regularisation part
%   ** ENRENAME the names of the different regularisers
%

if nargin < 5 || isempty(W); W = ones(size(x), 'like',real(x));end
if nargin < 6 || isempty(sep); sep=(nargout>1);end
if nargin < 7 || isempty(deb); deb=0;end

%%%DATA FIDELITY 
EnFi = 0.5 *  multDimSum(computeEnergy(y,x,E,[],EH))/numel(y) ; 

%%%REGULARISATION
EnRe = [];
EnReNames = {};

if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg,'lambda_sparse') && ~isequal(E.Dl.Reg.lambda_sparse,0) %L1 norm 
    EnRe = cat(1,EnRe , 0.5 * E.Dl.Reg.lambda_sparse * multDimSum( abs(bsxfun(@times,ones(size(W)),E.Dl.D))) ) /numel(E.Dl.D);
    EnReNames = cat(1,EnReNames, 'sparse');
end

if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg,'lambda_smooth') && ~isequal(E.Dl.Reg.lambda_smooth,0) %Smoothness 
    EnRe = cat(1,EnRe, 0.5 * E.Dl.Reg.lambda_smooth * ( normm( sqrt(W) .* FiniteDiff(E.Dl.D,1)) + ...  
                                                        normm( sqrt(W) .* FiniteDiff(E.Dl.D,2)) + ...
                                                        normm( sqrt(W) .* FiniteDiff(E.Dl.D,3)) )/numel(E.Dl.D) );
    EnReNames = cat(1,EnReNames, 'smooth');
end

if isfield(E.Dl,'Reg') && isfield(E.Dl.Reg,'lambda_ridge') && ~isequal(E.Dl.Reg.lambda_ridge,0)%Ridge 
    EnRe = cat(1,EnRe, 0.5 * E.Dl.Reg.lambda_ridge * normm(E.Dl.D) /numel(E.Dl.D));
    EnReNames = cat(1,EnReNames, 'ridge');
end

if deb>0
    fprintf('energyTaylor:: Relative to Data Fidelity:\n')  ;   
    names = strcat(EnReNames,{'/'});
    fprintf('%s = %s\n' , strcat(names{:}), sprintf('%d ', (EnRe/EnFi)')  );       
end

%%%GATHER
EnFi = gather(EnFi);
EnRe = gather(EnRe);

if sep==0 && ~isempty(EnRe); EnFi = EnFi + multDimSum(EnRe,1); end

end

