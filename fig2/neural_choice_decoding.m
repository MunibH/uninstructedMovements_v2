clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
utilspth = '/Users/munib/Economo-Lab/code/uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')))
rmpath(genpath(fullfile(utilspth,'figNP/')))
rmpath(genpath(fullfile(utilspth,'mc_stim/')))


% add paths for figure specific functions
addpath(genpath(pwd))

clc

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off

params.tmin = -3;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


% params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
%     {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance


params.advance_movement = 0.0;


%% SPECIFY DATA TO LOAD


datapth = '/Users/Munib/Documents/Economo-Lab/data/';
datapth = '/users/munib/Economo-Lab/data/';


meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);


% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);


%% choice decoding from dlc features

clearvars -except datapth kin me meta obj params

% params
rez.nFolds = 4; % number of iterations (bootstrap)

rez.binSize = 40; % ms
rez.dt = floor(rez.binSize / (params(1).dt*1000)); % samples
rez.tm = obj(1).time(1:rez.dt:numel(obj(1).time));
rez.numT = numel(rez.tm);

rez.train = 1; % fraction of trials to use for training (1-train for testing)

rez.nShuffles = 2;

% match number of right and left hits, and right and left misses
cond2use = 2:5;
hitcond = [1 3];
misscond = [2 4];



for sessix = 1:numel(obj)
    disp(['Decoding session ' num2str(sessix) ' / ' num2str(numel(obj))])

    % trials

    trials_cond = params(sessix).trialid(cond2use);

    minHitTrials = cellfun(@(x) numel(x),trials_cond(hitcond), 'UniformOutput',false);
    nhits = min(cell2mat(minHitTrials));

    minMissTrials = cellfun(@(x) numel(x),trials_cond(misscond), 'UniformOutput',false);
    nmiss = min(cell2mat(minMissTrials));


    trials_hit = cellfun(@(x) randsample(x,nhits), trials_cond(hitcond), 'UniformOutput', false);
    trialsHit = cell2mat(trials_hit);
    trialsHit = trialsHit(:);

    trials_miss = cellfun(@(x) randsample(x,nmiss), trials_cond(misscond), 'UniformOutput', false);
    trialsMiss = cell2mat(trials_miss);
    trialsMiss = trialsMiss(:);

    %         trials.all = [trialsHit ; trialsMiss];
    trials.all = trialsHit;

    % labels (1 for right choice, 0 for left choice)
    Y = [ones(nhits,1) ; -ones(nhits,1) ; ones(nmiss,1) ; -ones(nmiss,1)]; % right hits, left hits, left miss, right miss

    % input
    X = obj(sessix).trialdat(:,:,trials.all);
    X = permute(X,[1 3 2]); % (time,trials, clu)
    Xmax = max(X(:));
    Xmin = min(X(:));
    Range = Xmax - Xmin;
    Xnrm = ((X - Xmin)/Range - 0.5) * 2; 
    X = Xnrm; clear Xnrm;

    % train/test split
    [trials.train,trials.trainidx] = datasample(trials.all,round(numel(trials.all)*rez.train),'Replace',false);
    trials.testidx = find(~ismember(trials.all,trials.train));
    trials.test = trials.all(trials.testidx);

    in.train.y = Y(trials.trainidx);
    in.test.y  = Y(trials.testidx);
    in.train.X = X(:,trials.trainidx,:);
    in.test.X  = X(:,trials.testidx,:);

    % decoding

    acc(:,sessix) = DLC_ChoiceDecoder(in,rez,trials);


    % shuffle labels for a 'null' distribution


    Y = randsample(Y,numel(Y));

    % train/test split

    in.train.y = Y(trials.trainidx);
    in.test.y  = Y(trials.testidx);

    for ishuf = 1:rez.nShuffles
        acc_shuf(:,sessix,ishuf) = DLC_ChoiceDecoder(in,rez,trials);
    end


end


acc_shuf_ = reshape(acc_shuf,size(acc_shuf,1),size(acc_shuf,2)*size(acc_shuf,3));


%% plot
close all

cols = {'k',[0.6,0.6,0.6]};

alph = 0.15;

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);
trialStart = mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue);

f = figure;
ax = gca;
hold on;

ctrl = mySmooth(acc, 11,'reflect');
shuffed = mySmooth(acc_shuf_, 101,'reflect');


shadedErrorBar(rez.tm(1:end-1),mean(ctrl,2),getCI(ctrl),{'Color',cols{1},'LineWidth',2},alph,ax)
shadedErrorBar(rez.tm(1:end-1),mean(shuffed,2),getCI(shuffed,0),{'Color',cols{2},'LineWidth',2},alph,ax)

% shadedErrorBar(rez.tm(1:end-1),mean(acc,2),std(acc,[],2)./sqrt(numel(obj)),{'Color',cols{1},'LineWidth',2},alph,ax)
% shadedErrorBar(rez.tm(1:end-1),mean(acc_shuf_,2),std(acc_shuf_,[],2)./sqrt(numel(obj)),{'Color',cols{2},'LineWidth',2},alph,ax)
xline(0,'k--','LineWidth',1)
xline(sample,'k--','LineWidth',1)
xline(delay,'k--','LineWidth',1)
xline(trialStart,'k--','LineWidth',1)
xlim([trialStart, params(1).tmax-0.2])
ylim([ax.YLim(1) 1])

xlabel('Time from go cue (s)')
ylabel([num2str(rez.nFolds) '-Fold CV Accuracy'])
title('Neural choice decoding','FontSize',8)

h = zeros(2, 1);
for i = 1:numel(h)
    h(i) = plot(NaN,NaN,'-','Color',cols{i},'LineWidth',2);
end
legString = {'Neural data','Shuffled labels'};

leg = legend(h, legString);
leg.EdgeColor = 'none';
leg.Location = 'best';
leg.Color = 'none';
ax.FontSize = 10;




