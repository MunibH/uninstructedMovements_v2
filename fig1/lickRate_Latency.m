clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig3/')))
rmpath(genpath(fullfile(utilspth,'mc_stim/')))
rmpath(genpath(fullfile(utilspth,'MotionMapper/')))


% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'R&hit&~stim.enable&~autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&hit&~stim.enable&autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % left hits, no stim, aw off

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance


params.advance_movement = 0.0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);


params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

%% LICK LATENCY
gocue = mode(obj(1).bp.ev.goCue);
latency = cell(numel(params(1).condition),1);
for j = 1:numel(params(1).condition)
    latency{j} = [];
    for i = 1:numel(meta)
        trix = params(i).trialid{j};
        lickR = obj(i).bp.ev.lickR(trix);
        lickL = obj(i).bp.ev.lickL(trix);
        licktimes = cellfun( @(x,y) sort([x,y]), lickR, lickL, 'UniformOutput', false );
        latency{j} = [latency{j} ; cell2mat(cellfun(@(x) x(1),licktimes, 'UniformOutput',false )) - gocue];
    end
    latency{j}(latency{j}<0) = [];
end


figure; hold on;
rng(pi) % just to reproduce the random data I used
div = 1.3;

clrs = getColors();
fns = fieldnames(clrs);
for i = 1:numel(fns)
    cols{i} = clrs.(fns{i});
end

xs = [1 2 4 5];
for i = 1:numel(latency)
    b(i) = bar(xs(i),median(latency{i}));
    b(i).FaceColor = cols{i};
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.8;
%     vs(i) = scatter(xs(i)*ones(size(latency{i})),latency{i},60,'MarkerFaceColor',cols{i}./div,...
%         'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25);
    errorbar(b(i).XEndPoints,median(latency{i}),std(latency{i}),'LineStyle','none','Color','k','LineWidth',1)


end


xticklabels([" " "Right 2AFC" "Left 2AFC" " " "Right AW" "Left AW"])
ylabel("First lick latency (s)")
ylim([0,1])
ax = gca;
ax.FontSize = 12;



%% LICK RATE
gocue = mode(obj(1).bp.ev.goCue);
lickrate = cell(numel(params(1).condition),1);
for j = 1:numel(params(1).condition)
    lickrate{j} = [];
    for i = 1:numel(meta)
        trix = params(i).trialid{j};
        lickR = obj(i).bp.ev.lickR(trix);
        lickL = obj(i).bp.ev.lickL(trix);
        licktimes = cellfun( @(x,y) sort([x,y])-gocue, lickR, lickL, 'UniformOutput', false );
        licktimes = cellfun(@(x) x(x>0),licktimes, 'UniformOutput',false);
        lr = 1./cell2mat(cellfun(@(x) mean(diff(x)), licktimes, 'UniformOutput',false)); % hz
        lickrate{j} = [lickrate{j} ; lr];
    end
end


figure; hold on;
rng(pi) % just to reproduce the random data I used
div = 1.3;

clrs = getColors();
fns = fieldnames(clrs);
for i = 1:numel(fns)
    cols{i} = clrs.(fns{i});
end

xs = [1 2 4 5];
for i = 1:numel(lickrate)
    b(i) = bar(xs(i),nanmean(lickrate{i}));
    b(i).FaceColor = cols{i};
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.8;
%     vs(i) = scatter(xs(i)*ones(size(lickrate{i})),lickrate{i},60,'MarkerFaceColor',cols{i}./div,...
%         'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25);
    errorbar(b(i).XEndPoints,nanmean(lickrate{i}),nanstd(lickrate{i}),'LineStyle','none','Color','k','LineWidth',1)


end


xticklabels([" " "Right 2AFC" "Left 2AFC" " " "Right AW" "Left AW"])
ylabel("Lick rate (Hz)")
% ylim([0,1])
ax = gca;
ax.FontSize = 12;

























