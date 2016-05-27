%
% See COPYRIGHT for license information 
%

clear variables; close all;

addpath('mfmc/');
addpath('shortcol/');

%% parameters
pList = [1e+4, 1e+5, 1e+6]; % computational budget
% models
modelList = {
    @shortcolFOM, ...
    @shortcolROM02, ...
    @shortcolROM05, ...
    @shortcolROM04, ...
    @shortcolROM03, ...
    };

%% generate samples to estimate correlation and variance
N = 1e+2; % number of samples for corr/var estimation
mu = drawSamples(N);
% eval models and store in matrix X
X = zeros(N, length(modelList));
for i=1:length(modelList)
    X(:, i) = modelList{i}(mu);
end
% assign costs (replace with actual time measurements)
w = [100, 50, 20, 10, 5]';

%% generate optimal ordering of models
[bestSet, bestV, vMC] = optimalOrderCorr(X, w);
modelList = {modelList{bestSet}};
w = w(bestSet);
X = X(:, bestSet);
disp(['Selected models ', num2str(bestSet)]);

%% estimate mean of output
nrRuns = 100; % average of nrRuns runs
yMean = zeros(length(pList), nrRuns);
costs = zeros(length(pList), 1);
for run=1:nrRuns 
    disp(['Run ', num2str(run), ' of ', num2str(nrRuns)]);
    for j=1:length(pList) % compute estimate for different budgets
        % compute optimal number of model evals/coefficients
        [Mlevel, coeff, v] = optiMlevelCorr(X, w, pList(j));
        disp(['Number of model evals ', num2str(Mlevel')]);
        
        p = drawSamples(Mlevel(end));
        % FOM model
        disp(['model 1 ', num2str(Mlevel(1))]);
        y = modelList{1}(p(1:Mlevel(1), :));
        yMeanCur = coeff(1)*mean(y);
        costs(j) = Mlevel(1)*w(1);
    
        % iterate over all other models
        for i=2:length(modelList)
            disp(['model ', num2str(i), ' ', num2str(Mlevel(i))]);
            y = modelList{i}(p(1:Mlevel(i), :));
            yMeanCur = yMeanCur + coeff(i)*(mean(y) - mean(y(1:Mlevel(i - 1))));
            costs(j) = costs(j) + Mlevel(i)*w(i);
        end
        yMean(j, run) = yMeanCur;
    end
end

%% Plot
% estimate with FOM only (modelList{1})
yMeanFOM = zeros(length(pList), nrRuns);
costsFOM = zeros(length(pList), 1);
for run=1:nrRuns
    for j=1:length(pList)
        M = ceil(pList(j)/w(1)); % number of samples
        p = drawSamples(M);
        yMeanFOM(j, run) = mean(modelList{1}(p));
        costsFOM(j) = M*w(1);
    end
end
% estimate reference mean
yMeanRef = mean(modelList{1}(drawSamples(1e+7)));
% estimate MSE
mseMFMC = zeros(length(pList), 1);
mseMC = zeros(length(pList), 1);
for j=1:length(pList)
    mseMFMC(j) = mean((yMeanRef - yMean(j, :)).^2);
    mseMC(j) = mean((yMeanRef - yMeanFOM(j, :)).^2);
end
% plot
loglog(costs, mseMC, '-bo');
hold on;
loglog(costsFOM, mseMFMC, '-xr');
axis([pList(1) pList(end) -Inf Inf]);
xlabel('costs (budget)');
ylabel('estimated MSE');
legend('Monte Carlo', 'Multifidelity Monte Carlo');


        
        
        
        
        
        
