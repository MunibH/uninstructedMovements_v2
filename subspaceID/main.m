
% clear,clc,close all

% add manopt to path
cd_ = pwd; % path to manopt 
addpath(genpath(fullfile(cd_,'manopt')))

%% LOAD SOME DATA

% if working with trial-averaged data, data should be formatted as
%     psth = (bins, neurons, nConditions)

load(fullfile(cd_,'testData.mat')); % (bins,neurons,nConditions)
nBins = size(psth,1);
nNeurons = size(psth,2);
nConditions = size(psth,3);

%% PREPROCESS DATA AS IN ELSAYED 2016

% soft-normalize
% firing rate of a neuron, x of size (time,1), is transformed as:
% x_norm = x / (lambda + max(x) - min(x))
lambda = 5; 
snpsth = psth ./ (lambda + max(psth) - min(psth));

% mean center each neuron (across conditions)
mcpsth = zeros(size(psth));
for iclu = 1:nNeurons
    % find mean at each time point for each condition (time,1)
    mean_across_cond = mean(snpsth(:,iclu,:),3);
    mcpsth(:,iclu,:) = snpsth(:,iclu,:) - repmat(mean_across_cond,[1,1,size(snpsth,3)]);
end

%% DEFINE EPOCHS

% null
nullix = 70:120; % time bins to use to estimate null subspace

% potent
potentix = 125:175; % time bins to use to estimate null subspace

%% CONCATENATE CONDITIONS

% we will estimate covariances from null and potent times using data
% formatted as:
%        (nBins*nConds,nNeurons)

nullpsth = mcpsth(nullix,:,1);
potentpsth = mcpsth(potentix,:,1);
for i = 2:nConditions
    nullpsth = cat(1,nullpsth,mcpsth(nullix,:,i));
    potentpsth = cat(1,potentpsth,mcpsth(potentix,:,i));
end

%% COVARIANCE AND DIMENSIONALITY

% covariances
Cnull = cov(nullpsth);
Cpotent = cov(potentpsth);

var2explain = 90; % dimensionality of subspaces will be equal to how many PCs needed to explain this much variance

[~,~,~,~,explained] = pca(nullpsth);
dNull = find(cumsum(explained)>=var2explain,1,'first');
[~,~,~,~,explained] = pca(potentpsth);
dPotent = find(cumsum(explained)>=var2explain,1,'first');

dMax = max([dNull, dPotent]);

%% MAIN OPTIMIZATION STEP
rng(101) % for reproducibility

alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
[Q,~,P,~,~] = orthogonal_subspaces(Cpotent,dPotent,Cnull,dNull,alpha);


Qpotent = Q*P{1}; % size = (nNeurons,dPotent)
Qnull = Q*P{2};   %        (nNeurons,dNull)

%% PROJECTIONS

proj.potent = tensorprod(mcpsth,Qpotent,2,1); % size = (nBins,nConditions,dPotent)
proj.null = tensorprod(mcpsth,Qnull,2,1);     % size = (nBins,nConditions,dNull)

% plot first dimension (not ordered in any way) of null and potent projs
f = figure;

ax = subplot(2,1,1);
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.TickLength = ax.TickLength .* 2;
hold on;
plot(squeeze(proj.potent(:,1,1)), 'linewidth',2,'color','b')
plot(squeeze(proj.potent(:,2,1)), 'linewidth',2,'color','r')
title('potent')
ylabel('projection (a.u.)')
xlabel('time bins')

ax = subplot(2,1,2);
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.TickLength = ax.TickLength .* 2;
hold on;
plot(squeeze(proj.null(:,1,1)), 'linewidth',2,'color','b')
plot(squeeze(proj.null(:,2,1)), 'linewidth',2,'color','r')
title('null')
ylabel('projection (a.u.)')
xlabel('time bins')

legend('right trials','left trials','edgecolor','none','location','best')
