%
% See COPYRIGHT for license information 
%
function [ Mlevel, coeff, v ] = optiMlevelCorrHelper( sigma, rho, w, p )
%OPTIMLEVEL Computes Mlevel (con var) taking correlation into account
% In
%   sigma   ...     square roots of variance
%   rho     ...     correlation matrix
%   w       ...     costs vector
%   p       ...     budget
% Out
%   Mlevel  ...     number of model evaluations
%   coeff   ...     coefficients
%   v       ...     variance

k = length(sigma);
sigma = sigma(:);
w = w(:);

% sanity check (see paper)
assert(all(rho(2:k - 1, 1).^2 - rho(3:k, 1).^2 > 0));

% coefficients
coeff = sigma(1)*rho(1, :)'./sigma;

% r
r = zeros(k, 1);
r(1) = 1;
r(2:k-1) = sqrt(w(1)*(rho(2:k-1, 1).^2 - rho(3:k, 1).^2)./(w(2:k-1).*(1 - rho(1, 2)^2)));
r(k) = sqrt((w(1)*rho(1, k)^2)/(w(k)*(1 - rho(1, 2)^2)));

m = zeros(k, 1);
m(1) = p/(r'*w);
m(2:end) = r(2:end)*m(1);

v = sigma(1)^2/m(1) + sum((1./m(1:end - 1) - ...
    1./m(2:end)).*(coeff(2:end).^2.*sigma(2:end).^2 - 2*coeff(2:end).*rho(2:end, 1).*sigma(1).*sigma(2:end)));


Mlevel = floor(m);

end
