clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
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
params.condition(1)     = {'R&~stim.enable&~autowater&((1:Ntrials)>20)'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&~stim.enable&~autowater&((1:Ntrials)>20)'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&~stim.enable&autowater&((1:Ntrials)>20)'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&~stim.enable&autowater&((1:Ntrials)>20)'};             % left hits, no stim, aw off

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

params.behav_only = 1;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB15_ALMVideo(meta,datapth);


params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params, params.behav_only);


%% LICK RASTER

close all

sessix = 2; % jeb7 4-29


clrs = getColors();
fns = fieldnames(clrs);
for i = 1:numel(fns)
    cols{i} = clrs.(fns{i});
end


cond2use = 1:4;

lw = 1.5;
lwx = 0.5;
ms = 3;

goCue = 0;
sample = mode(obj(sessix).bp.ev.sample) - mode(obj(sessix).bp.ev.goCue);
delay = mode(obj(sessix).bp.ev.delay) - mode(obj(sessix).bp.ev.goCue);



nTrials{1} = numel(params(sessix).trialid{3});
nTrials{2} = numel(params(sessix).trialid{1});



f = figure; hold on;
f.Position = [520   614   189   365];
trialOffset = 1;


% right trials (2afc + aw)
cond2use = [3 1];
for cix = 1:numel(cond2use)
    if cix > 1; trialOffset = trialOffset + 10; end

    for trix = 1:nTrials{cix}
        check1 = 0;
        check2 = 0;

        trial = params(sessix).trialid{cond2use(cix)}(trix);
        lickL =  obj(sessix).bp.ev.lickL{trial} - obj(sessix).bp.ev.goCue(trial);
        lickL(lickL > 2) = [];
        lickR =  obj(sessix).bp.ev.lickR{trial} - obj(sessix).bp.ev.goCue(trial);
        lickR(lickR > 2) = [];

        plot([sample sample], trialOffset+[-0.5 0.5], 'k--', 'LineWidth', lwx);
        plot([delay delay], trialOffset+[-0.5 0.5], 'k--', 'LineWidth', lwx);
        plot([goCue goCue], trialOffset+[-0.5 0.5], 'k--', 'LineWidth', lwx);

        if ~isempty(lickL)
            plot(lickL, trialOffset*ones(size(lickL)), '.', 'Color', cols{cond2use(cix)+1}, 'MarkerSize',ms);
            check1 = 1;
        end

        if ~isempty(lickR)
            plot(lickR, trialOffset*ones(size(lickR)), '.', 'Color', cols{cond2use(cix)}, 'MarkerSize',ms);
            check2 = 1;
        end

        if obj(sessix).bp.hit(trial)
            fill([2.4 2.65 2.65 2.4], [trialOffset-0.5 trialOffset-0.5 trialOffset+0.5 trialOffset+0.5], [150 150 150]./255,'EdgeColor','none')
        elseif obj(sessix).bp.miss(trial)
            fill([2.4 2.65 2.65 2.4], [trialOffset-0.5 trialOffset-0.5 trialOffset+0.5 trialOffset+0.5], [0 0 0]./255,'EdgeColor','none')
        elseif obj(sessix).bp.no(trial)
            %         fill([2.4 2.65 2.65 2.4], [trialOffset trialOffset trialOffset+1 trialOffset+1], [1 1 1]./255,'EdgeColor','k')
        end

        if check1 || check2
            trialOffset = trialOffset + 1;
        end

    end

end
xlim([-2.5 2.7]);
ylim([0 160])
xlabel('Time (s) from go cue')
ylabel('Trials')
ax = gca;
ax.FontSize = 9;

%

nTrials{1} = numel(params(sessix).trialid{4});
nTrials{2} = numel(params(sessix).trialid{2});


f = figure; hold on;
f.Position = [520   614   189   365];
trialOffset = 1;

% left trials (2afc + aw)
cond2use = [4 2];
for cix = 1:numel(cond2use)
    if cix > 1; trialOffset = trialOffset + 10; end

    for trix = 1:nTrials{cix}
        check1 = 0;
        check2 = 0;

        trial = params(sessix).trialid{cond2use(cix)}(trix);
        lickL =  obj(sessix).bp.ev.lickL{trial} - obj(sessix).bp.ev.goCue(trial);
        lickL(lickL > 2) = [];
        lickR =  obj(sessix).bp.ev.lickR{trial} - obj(sessix).bp.ev.goCue(trial);
        lickR(lickR > 2) = [];

        plot([sample sample], trialOffset+[-0.5 0.5], 'k--', 'LineWidth', lwx);
        plot([delay delay], trialOffset+[-0.5 0.5], 'k--', 'LineWidth', lwx);
        plot([goCue goCue], trialOffset+[-0.5 0.5], 'k--', 'LineWidth', lwx);

        if ~isempty(lickL)
            plot(lickL, trialOffset*ones(size(lickL)), '.', 'Color', cols{cond2use(cix)}, 'MarkerSize',ms);
            check1 = 1;
        end

        if ~isempty(lickR)
            plot(lickR, trialOffset*ones(size(lickR)), '.', 'Color', cols{cond2use(cix)-1}, 'MarkerSize',ms);
            check2 = 1;
        end

        if obj(sessix).bp.hit(trial)
            fill([2.4 2.65 2.65 2.4], [trialOffset-0.5 trialOffset-0.5 trialOffset+0.5 trialOffset+0.5], [150 150 150]./255,'EdgeColor','none')
        elseif obj(sessix).bp.miss(trial)
            fill([2.4 2.65 2.65 2.4], [trialOffset-0.5 trialOffset-0.5 trialOffset+0.5 trialOffset+0.5], [0 0 0]./255,'EdgeColor','none')
        elseif obj(sessix).bp.no(trial)
            %         fill([2.4 2.65 2.65 2.4], [trialOffset trialOffset trialOffset+1 trialOffset+1], [1 1 1]./255,'EdgeColor','k')
        end

        if check1 || check2
            trialOffset = trialOffset + 1;
        end

    end

end
xlim([-2.5 2.7]);
ylim([0 110])
xlabel('Time (s) from go cue')
ylabel('Trials')
ax = gca;
ax.FontSize = 9;












