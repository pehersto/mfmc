%
% See COPYRIGHT for license information 
%
function [ y ] = shortcolFOM( param )
%SHORTCOLFOM Summary of this function goes here
%   Detailed explanation goes here

b = param(:, 1);
h = param(:, 2);
P = param(:, 3);
M = param(:, 4);
Y = param(:, 5);

y = 1 - 4*M./(b.*(h.^2).*Y) - (P./(b.*h.*Y)).^2;

end

