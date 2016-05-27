%
% See COPYRIGHT for license information 
%
function [ y ] = shortcolROM03( param )
%SHORTCOLROM Surrogate model of short column FOM
% In
%   param       ...     parameter vector
% Out
%   y           ...     output

b = param(:, 1);
h = param(:, 2);
P = param(:, 3);
M = param(:, 4);
Y = param(:, 5);

y = 1 - M./(b.*(h.^2).*Y) - (P./(b.*h.*Y)).^2;
end

