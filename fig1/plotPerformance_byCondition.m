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

% performance by mouse on AFC and AW trials

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'~stim.enable&~autowater&((1:Ntrials)>20)'};             % afc
params.condition(end+1) = {'~stim.enable&autowater&((1:Ntrials)>20)'};              % aw

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

% --- ALM --- 
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);
% meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params, params.behav_only);

%% PERFORMANCE
rez = getPerformanceByCondition(meta,obj,params);

% concatenate perf for each mouse
anms = unique({meta.anm});
perf = cell(numel(unique(anms)),1);
for i = 1:numel(meta)
    anmix = find(ismember(anms,meta(i).anm));
    perf{anmix} = [perf{anmix} ; rez(i).perf];
end
% average across sessions for each mouse
perf = cellfun(@(x) mean(x,1)*100, perf, 'UniformOutput',false);

perf = cell2mat(perf); % (animals, cond)


%% PLOT

close all

figure;
ax = gca;
hold on;

clrs = getColors();
cols{1} = clrs.afc;
cols{2} = clrs.aw;

xs = [0 1];
for i = 1:size(perf,2)
    b(i) = bar(xs(i),mean(perf(:,i)));
    b(i).FaceColor = cols{i};
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.5;
    vs(i) = scatter(randn(numel(anms),1) * 0.1 + xs(i)*ones(size(perf(:,i))),perf(:,i),30,'MarkerFaceColor',cols{i},...
        'MarkerEdgeColor','k','LineWidth',1);%,'XJitter','randn','XJitterWidth',0.25);
    errorbar(b(i).XEndPoints,mean(perf(:,i)),std(perf(:,i)),'LineStyle','none','Color','k','LineWidth',1)

    xs_(:,i) = vs(i).XData';
    ys_(:,i) = perf(:,i);
end

for i = 1:size(xs_,1)
    patchline(xs_(i,:),ys_(i,:),'EdgeAlpha',0.4)
end


ax.XTick = xs;
xticklabels(["2AFC" "AW"])
ylabel("Performance (%)")
ylim([0,100])
ax = gca;
ax.FontSize = 12;


