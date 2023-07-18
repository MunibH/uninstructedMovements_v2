
% clear,clc,close all

% add manopt to path
cd_ = pwd; % path to manopt 
addpath(genpath(fullfile(cd_,'manopt')))

%% LOAD SOME DATA

% if working with single-trial data, data should be formatted as
%    (time bins, trials, neurons)

load(fullfile(cd_,'exampleData.mat')); 

% this data is from an example session:
%   brain region - Left ALM
%   task - delayed response (DR) + water-cued (WC)
%   

% exampleData is a struct containing fields:
%   time - time bins corresponding to first dimension of fields 'seq' and
%          'motionEnergy' (aligned to go cue/water drop)
%   seq  - single trial neural firing rates in 5 ms bins and smoothed with
%          causal gaussian kernel (35 ms s.d.). Shape = (time bins, trials,
%          neurons)
%   motionEnergy - motion energy is used as a measure of amount of movement
%                  at a given frame. The time series has been resampled to 
%                  match neural data. Shape = (time bins, trials)
%   moveMask - logical array of shape (time bins, trials). 0 indicates
%              stationarity, 1 indicates moving. This mask was produced from
%              exampleData.motionEnergy using a manually set threshold
%   rightCorrectTrials - index of correct, right DR trials (corresponds to
%                        dimension 2 of exampleData.seq/motionEnergy/moveMask)
%   leftCorrectTrials  - index of correct, left DR trials

nBins = size(exampleData.seq,1);
nTrials = size(exampleData.seq,2);
nNeurons = size(exampleData.seq,3);

%% PREPROCESS DATA
% choice of normalization/standardization is up to you, here just zscoring


temp = reshape(exampleData.seq,nBins*nTrials,nNeurons);
N.full_cat = zscore(temp); % (time bins * nTrials, nNeurons)
N.full= reshape(N.full_cat,nBins,nTrials,nNeurons); % (time bins, nTrials, nNeurons)


%% DEFINE DATA TO USE FOR NULL AND POTENT SUBSPACES

% for the null subspace, we will use all time points, from all trials, in
% which the animal was stationary.
moveMask = exampleData.moveMask(:); %(time bins * nTrials)
N.null = N.full_cat(~moveMask,:);

% for the potent subspace, we will use all time points, from all trials, in
% which the anima was moving.
N.potent = N.full_cat(moveMask,:);

rez.N = N; % put data into a results struct

%% COVARIANCE AND DIMENSIONALITY

% covariances
rez.cov.null = cov(N.null);
rez.cov.potent = cov(N.potent);

% dimensionality of subspaces
% here I am hard-coding the number of dimensions to be 20. However, one
% could perform further dimensionality reduction or keep more dimensions.
% In our experience, dimensionality >  20 takes an incredibly long time to 
% optimize over.

rez.dNull = 10;
rez.dPotent = 10;

rez.dMax = max([rez.dNull, rez.dPotent]);

%% MAIN OPTIMIZATION STEP
rng(101) % for reproducibility

alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
[Q,~,P,~,~] = orthogonal_subspaces(rez.cov.potent,rez.dPotent,rez.cov.null,rez.dNull,alpha);

% manopt will provide a warning that we haven't provided the Hessian,
% that's ok. Providing the Hessian of the objective function will speed up
% computations, however.

% Q contains the dimensions of both the null and potent subspaces. Here, we
% extract the columns corresponding to null and potent subspaces.
rez.Q.potent = Q*P{1}; % size = (nNeurons,dPotent)
rez.Q.null = Q*P{2};   %        (nNeurons,dNull)

%% PROJECTIONS

% now that we have null and potent subspaces, we can project our neural
% activity onto the subspaces separately.

rez.proj.potent = tensorprod(N.full,rez.Q.potent,3,1); % size = (nBins,nTrials,dPotent)
rez.proj.null = tensorprod(N.full,rez.Q.null,3,1);     % size = (nBins,nTrials,dNull)

%%

% plot the sum squared magnitude across dimensions of activity within the
% null and potent subspaces - only going to plot right and left correct DR
% trials

r = exampleData.rightCorrectTrials;
l = exampleData.leftCorrectTrials;

plt.null.right = sum(rez.proj.null(:,r,:).^2,3); 
plt.null.left = sum(rez.proj.null(:,l,:).^2,3); 
plt.potent.right = sum(rez.proj.potent.^2,3); 
plt.potent.left = sum(rez.proj.potent(:,l,:).^2,3); 


%% PLOT SINGLE TRIALS

f = figure;

xlims = [-2 2];

% plot motion energy
ax = subplot(1,3,1);
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.TickLength = ax.TickLength .* 2;
hold on;
imagesc(exampleData.time,1:(numel(l)+numel(r)),exampleData.motionEnergy(:,[r;l])')
ax.YDir = 'normal';
line([ax.XLim(1) ax.XLim(2)], [numel(r) numel(r)] ,'color','w','linestyle','-') % right trials above line, left trials below
line([0 0], [ax.YLim(1) ax.YLim(2)],'color','w','linestyle','--') % go cue
title('Motion energy')
xlim(xlims);
xlabel('Time from go cue (s)')
ax.FontSize = 12;

% plot potent
ax = subplot(1,3,2);
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.TickLength = ax.TickLength .* 2;
hold on;
temp = cat(2,plt.potent.right,plt.potent.left);
imagesc(exampleData.time,1:(numel(l)+numel(r)),temp')
ax.YDir = 'normal';
line([ax.XLim(1) ax.XLim(2)], [numel(r) numel(r)] ,'color','w','linestyle','-') % right trials above line, left trials below
line([0 0], [ax.YLim(1) ax.YLim(2)],'color','w','linestyle','--') % go cue
title('Potent')
xlim(xlims);
ax.FontSize = 12;

% plot null
ax = subplot(1,3,3);
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.TickLength = ax.TickLength .* 2;
hold on;
temp = cat(2,plt.null.right,plt.null.left);
imagesc(exampleData.time,1:(numel(l)+numel(r)),temp')
ax.YDir = 'normal';
line([ax.XLim(1) ax.XLim(2)], [numel(r) numel(r)] ,'color','w','linestyle','-') % right trials above line, left trials below
line([0 0], [ax.YLim(1) ax.YLim(2)],'color','w','linestyle','--') % go cue
title('Null')
xlim(xlims);
ax.FontSize = 12;


%% PLOT TRIAL AVERAGED DATA

f = figure;

% plot null
ax = subplot(2,1,1);
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.TickLength = ax.TickLength .* 2;
hold on;
% plot right trials
mu = mean(plt.null.right,2); % mean across trials
se = std(plt.null.right,[],2) ./ sqrt(numel(r)); % standard error 
shadedErrorBar(exampleData.time,mu,se,{'Color','b','LineWidth',2},0.2,ax)
% plot left trials
mu = mean(plt.null.left,2); % mean across trials
se = std(plt.null.left,[],2) ./ sqrt(numel(l)); % standard error 
shadedErrorBar(exampleData.time,mu,se,{'Color','r','LineWidth',2},0.2,ax)

line([0 0], [ax.YLim(1) ax.YLim(2)],'color','k','linestyle','--') % go cue
title('Null')
ylabel('projection (a.u.)')
xlabel('Time from go cue (s)')
ax.FontSize = 12;
xlim(xlims)

% plot potent
ax = subplot(2,1,2);
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.TickLength = ax.TickLength .* 2;
hold on;
% plot right trials
mu = mean(plt.potent.right,2); % mean across trials
se = std(plt.potent.right,[],2) ./ sqrt(numel(r)); % standard error 
shadedErrorBar(exampleData.time,mu,se,{'Color','b','LineWidth',2},0.2,ax)
% plot left trials
mu = mean(plt.potent.left,2); % mean across trials
se = std(plt.potent.left,[],2) ./ sqrt(numel(l)); % standard error 
shadedErrorBar(exampleData.time,mu,se,{'Color','r','LineWidth',2},0.2,ax)

line([0 0], [ax.YLim(1) ax.YLim(2)],'color','k','linestyle','--') % go cue
title('potent')
ylabel('projection (a.u.)')
xlabel('Time from go cue (s)')
ax.FontSize = 12;
xlim(xlims)