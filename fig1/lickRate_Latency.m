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
% meta = loadJEB15_ALMVideo(meta,datapth);


params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

%% LICK LATENCY

close all

latency = cell(numel(meta),1);
cond2use = 1:4;
for sessix = 1:numel(meta)
    for cix = 1:numel(cond2use)
        latency{sessix}{cix} = nan((numel(params(sessix).trialid{cond2use(cix)})),1);
        for trix = 1:numel(params(sessix).trialid{cond2use(cix)})
            trial = params(sessix).trialid{cond2use(cix)}(trix);

            lickR = obj(sessix).bp.ev.lickR{trial};
            lickL = obj(sessix).bp.ev.lickL{trial};
            licktimes = sort([lickR lickL]) - obj(sessix).bp.ev.goCue(trial);
            licktimes(licktimes <= 0) = [];
            latency{sessix}{cix}(trix) = licktimes(1);
        end
    end
end

% concatenate latency for each mouse
anms = unique({meta.anm});
latency_by_mouse = cell(numel(unique(anms)),1);
for i = 1:numel(meta) % preallocating in this loop
    anmix = find(ismember(anms,meta(i).anm));
    for cix = 1:numel(cond2use)
        latency_by_mouse{anmix}{cix} = [];
    end
end
for i = 1:numel(meta)
    anmix = find(ismember(anms,meta(i).anm));
    for cix = 1:numel(cond2use)
        latency_by_mouse{anmix}{cix} = [latency_by_mouse{anmix}{cix} ; latency{i}{cix}];
    end
end
% average across trials/sessions for each mouse
for i = 1:numel(anms)
    latencies{i} = cellfun(@(x) mean(x), latency_by_mouse{i}, 'UniformOutput',false);                                                                                                                      
end

% concatenate data 
lates = zeros(numel(anms),numel(cond2use));
for sessix = 1:numel(anms)
    for cix = 1:numel(cond2use)
        lates(sessix,cix) = latencies{sessix}{cix};
    end
end

clrs = getColors();
fns = fieldnames(clrs);
for i = 1:numel(fns)
    cols{i} = clrs.(fns{i});
end
div = 1.3;

f = figure; hold on;
f.Position = [680   755   244   223];
xs = [1 2 4 5];
for i = 1:numel(xs)
    b(i) = bar(xs(i),mean(lates(:,i)));
    b(i).FaceColor = cols{i};
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.8;
    vs(i) = scatter(xs(i)*ones(size(lates(:,i))),lates(:,i),60,'MarkerFaceColor',cols{i}./div,...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25);
    errorbar(b(i).XEndPoints,mean(lates(:,i)),std(lates(:,i)),'LineStyle','none','Color','k','LineWidth',1)


end

% xticklabels(["Right 2AFC" "Left 2AFC" " " "Right AW" "Left AW"])
ylabel("First lick latency (s)")
% ylim([0,1])
ax = gca;
ax.FontSize = 12;
ax.XTick = [];


%% LICK RATE

lickrate = cell(numel(meta),1);
cond2use = 1:4;
for sessix = 1:numel(meta)
    for cix = 1:numel(cond2use)
        lickrate{sessix}{cix} = nan((numel(params(sessix).trialid{cond2use(cix)})),1);
        for trix = 1:numel(params(sessix).trialid{cond2use(cix)})
            trial = params(sessix).trialid{cond2use(cix)}(trix);

            lickR = obj(sessix).bp.ev.lickR{trial};
            lickL = obj(sessix).bp.ev.lickL{trial};
            licktimes = sort([lickR lickL]) - obj(sessix).bp.ev.goCue(trial);
            licktimes(licktimes <= 0) = [];
            lickrate{sessix}{cix}(trix) = 1 ./ mean(diff(licktimes)); % hz
        end
    end
end

% concatenate rates for each mouse
anms = unique({meta.anm});
rates_by_mouse = cell(numel(unique(anms)),1);
for i = 1:numel(meta) % preallocating in this loop
    anmix = find(ismember(anms,meta(i).anm));
    for cix = 1:numel(cond2use)
        rates_by_mouse{anmix}{cix} = [];
    end
end
for i = 1:numel(meta)
    anmix = find(ismember(anms,meta(i).anm));
    for cix = 1:numel(cond2use)
        rates_by_mouse{anmix}{cix} = [rates_by_mouse{anmix}{cix} ; lickrate{i}{cix}];
    end
end
% average across trials/sessions for each mouse
for i = 1:numel(anms)
    rates{i} = cellfun(@(x) mean(x), rates_by_mouse{i}, 'UniformOutput',false);                                                                                                                      
end

% concatenate data 
rates_ = zeros(numel(anms),numel(cond2use));
for sessix = 1:numel(anms)
    for cix = 1:numel(cond2use)
        rates_(sessix,cix) = rates{sessix}{cix};
    end
end

clrs = getColors();
fns = fieldnames(clrs);
for i = 1:numel(fns)
    cols{i} = clrs.(fns{i});
end
div = 1.3;

f = figure; hold on;
f.Position = [680   755   244   223];

xs = [1 2 4 5];
for i = 1:numel(xs)
    b(i) = bar(xs(i),nanmean(rates_(:,i)));
    b(i).FaceColor = cols{i};
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.8;
    vs(i) = scatter(xs(i)*ones(size(rates_(:,i))),rates_(:,i),60,'MarkerFaceColor',cols{i}./div,...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25);
    errorbar(b(i).XEndPoints,nanmean(rates_(:,i)),nanstd(rates_(:,i)),'LineStyle','none','Color','k','LineWidth',1)


end

% xticklabels(["Right 2AFC" "Left 2AFC" " " "Right AW" "Left AW"])
ylabel("Lick rate (Hz)")
% ylim([0,1])
ax = gca;
ax.FontSize = 12;
ax.XTick = [];


























