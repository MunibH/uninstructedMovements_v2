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
params.condition(end+1) = {'R&hit&~stim.enable&autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % right hits, no stim, aw on
params.condition(end+1) = {'L&hit&~stim.enable&autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % left hits, no stim, aw on

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
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

%% REACTION TIME (FIRST TONGUE APPEARANCE AFTER GO CUE)

close all

% rt = firstTongueRT(obj);
rt = firstJawRT(obj);

% group reaction time by mouse, condition, only for trials we want
anms = unique({meta.anm});
cond2use = 1:4;
for ianm = 1:numel(anms) % for each animal
    [~,mask] = patternMatchCellArray({meta(:).anm}',{anms{ianm}},'all');
    sessix = find(mask);
    anm_rt_cell{ianm} = rt(sessix);

    anm_rt{ianm} = cell(1,numel(cond2use));
    for isess = 1:numel(sessix)
        for icond = 1:numel(cond2use)
            trix = params(sessix(isess)).trialid{cond2use(icond)};
            anm_rt{ianm}{icond} = [anm_rt{ianm}{icond} anm_rt_cell{ianm}{isess}(trix)];
        end
    end

end


% get means by condition for each mouse
for anmix = 1:numel(anm_rt)
    temp = anm_rt{anmix};
    for icond = 1:numel(temp)
        mean_rt_by_mouse_by_cond{anmix}(icond) = nanmean(temp{icond});
    end
end

% concatenate conditions for all mice
rt_all_mice_by_cond = cell(numel(cond2use),1);
for icond = 1:numel(cond2use)
    rt_all_mice_by_cond{icond} = [];
    for anmix = 1:numel(anm_rt)
        temp = anm_rt{anmix}{icond};
        rt_all_mice_by_cond{icond} = [rt_all_mice_by_cond{icond} temp];
    end
end

%%
close all

clrs = getColors();
fns = fieldnames(clrs);
for i = 1:numel(fns)
    cols{i} = clrs.(fns{i});
end
div = 1.3;

f = figure; hold on;
f.Position = [680   739   270   239];
xs = [1 2 4 5];
for i = 1:numel(xs)
    b(i) = bar(xs(i),nanmean(rt_all_mice_by_cond{i}));
    b(i).FaceColor = cols{i};
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.8;

    for a = 1:numel(mean_rt_by_mouse_by_cond)
        vs(i) = scatter(randn(1) * 0.1 + xs(i)*ones(size(mean_rt_by_mouse_by_cond{a}(i))),mean_rt_by_mouse_by_cond{a}(i),20,'MarkerFaceColor',cols{i}./div,...
            'MarkerEdgeColor','k','LineWidth',1);
    end


    errorbar(b(i).XEndPoints,nanmean(rt_all_mice_by_cond{i}),nanstd(rt_all_mice_by_cond{i}),'LineStyle','none','Color','k','LineWidth',1)
end

ylabel("Reaction time (s)")
ylim([0,0.85])
ax = gca;
ax.FontSize = 12;
ax.XTick = [1 2 4 5];
xticklabels([ "Right 2AFC" "Left 2AFC" "Right AW" "Left AW"])



%% distribution of rts

allrt = cell2mat(rt');
figure; 
histogram(allrt,200,'EdgeColor','none')
xline(nanmedian(allrt),'k--')
ax = gca;
ax.FontSize = 13;
xlabel('Reaction time (s)')
ylabel('# of trials')


