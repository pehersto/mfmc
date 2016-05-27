%
% See COPYRIGHT for license information 
%
function [ param ] = drawSamples( M )
%DRAWSAMPLES Draw samples for short column example
% In
%   M       ...     number of samples
% Out
%   param   ...     samples

if(M < 1)
    param = [];
    return;
end

b = rand(M, 1)*(15 - 5) + 5;
h = rand(M, 1)*(25 - 15) + 15;
P = randn(M, 1)*100 + 500;
Mo = randn(M, 1)*400 + 2000;
Y = exp(randn(M, 1)*0.5 + 5);

param = [b h P Mo Y];


end

