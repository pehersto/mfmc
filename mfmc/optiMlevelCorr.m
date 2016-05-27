%
% See COPYRIGHT for license information 
%
function [ Mlevel, coeff, v ] = optiMlevelCorr( X, w, p )
%OPTIMLEVEL Computes Mlevel (con var) taking correlation into account
% In
%   X       ...     data matrix (each column corresponds to one model)
%   w       ...     runtime measurements
%   p       ...     budget, usually w(FOM)*Minit
% Out
%   Mlevel  ...     number of model evaluations
%   coeff   ...     coefficients
%   v       ...     variance

k = size(X, 2);
sigma = sqrt(var(X));
sigma = sigma(:);
rho = corr(X); % WARNING this is corr NOT cov!
w = w(:);
[Mlevel, coeff, v] = optiMlevelCorrHelper(sigma, rho, w, p);

end

