%
% See COPYRIGHT for license information 
%
function [ bestSet, bestV, vMC ] = optimalOrderCorr( X, w )
%OPTIMALORDER Algorithm 1 of paper
% In
%   X           ...     data matrix (each column corresponds to one model)
%   w           ...     runtime measurements
% Out
%   bestSet     ...     set of models that lead to lowest variance
%   bestV       ...     estimated variance of MFMC estimator w.r.t. bestSet
%   vMC         ...     variance of vanilla MC estimator

DEBUG = 0;

% define
k = size(X, 2);
sigma = sqrt(var(X));
sigma = sigma(:);
rho = corr(X); % WARNING this is corr NOT cov!
w = w(:);
vMC = sigma(1)^2*w(1)/10; % p = 10

% list of feasible sets of models
bestSet = 1;
bestV = vMC;
for i=2:k
    allSets = combnk(1:k, i);
    for j=1:size(allSets, 1)
        curSet = allSets(j, :);
        % check if FOM is there
        if(isempty(find(curSet == 1, 1)))
            MYDEBUG([num2str(curSet), ' no FOM'], DEBUG);
            continue;
        end
        % order set w.r.t. correlation coefficient
        [~, corrOrder] = sort(rho(1, curSet), 'descend');
        curSet = curSet(corrOrder);
        assert(curSet(1) == 1);
        
        % check if feasability condition is satisified
        if(~isFeasible(rho(1, curSet), w(curSet)))
            MYDEBUG([num2str(curSet), ' not feasible'], DEBUG);
            continue;
        end
        
        % compute variance
        [~, ~, v] = optiMlevelCorr(X(:, curSet), w(curSet), 10);
        myvar = sigma(1)^2/10*(sqrt(w(curSet))'*sqrt(rho(curSet, 1).^2 - [rho(curSet(2:end), 1).^2; 0]))^2;
        assert(abs((v - myvar)/v) < 1e-6); % numerical check (see paper)
        
        MYDEBUG([num2str(curSet), ' v = ', num2str(v)], DEBUG);
        if(bestV > v)
            bestSet = curSet;
            bestV = v;
        end
    end
end
    
end

function y = isFeasible(rho, w)
    rho = rho(:);
    w = w(:);
    rho = [rho; 0];
    y = 1;
    for i=2:length(rho)-1
        y = y*(w(i - 1)/w(i) > (rho(i - 1)^2 - rho(i)^2)/(rho(i)^2 - rho(i + 1)^2));
    end
end


