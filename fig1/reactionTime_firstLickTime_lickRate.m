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

%% REACTION TIME 

close all

% rt = firstTongueRT(obj);
rt = firstJawRT(obj);

% group  by mouse, condition, only for trials we want
anms = unique({meta.anm});
cond2use = 1:2;
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
        mean_rt_by_mouse_by_cond(anmix,icond) = nanmean(temp{icond}); % (anm,cond)
    end
end

%% FIRST LICK TIME (FIRST LICKPORT CONTACT)

firstLT = firstLickTime(obj);

% group  by mouse, condition, only for trials we want
anms = unique({meta.anm});
cond2use = 1:2;
for ianm = 1:numel(anms) % for each animal
    [~,mask] = patternMatchCellArray({meta(:).anm}',{anms{ianm}},'all');
    sessix = find(mask);
    anm_LT_cell{ianm} = firstLT(sessix);

    anm_LT{ianm} = cell(1,numel(cond2use));
    for isess = 1:numel(sessix)
        for icond = 1:numel(cond2use)
            trix = params(sessix(isess)).trialid{cond2use(icond)};
            anm_LT{ianm}{icond} = [anm_LT{ianm}{icond} anm_LT_cell{ianm}{isess}(trix)];
        end
    end

end

% get means by condition for each mouse
for anmix = 1:numel(anm_LT)
    temp = anm_LT{anmix};
    for icond = 1:numel(temp)
        mean_LT_by_mouse_by_cond(anmix,icond) = nanmean(temp{icond}); % (anm,cond)
    end
end



%%
close all

toplot = cat(2,mean_rt_by_mouse_by_cond,mean_LT_by_mouse_by_cond);

clrs = getColors();
cols{1} = clrs.afc;
cols{2} = clrs.aw;
cols{3} = clrs.afc;
cols{4} = clrs.aw;

f = figure; hold on;
f.Position = [680   739   270   239];

xs = 1:size(toplot,2);
for i = 1:numel(xs)
    b(i) = bar(xs(i),nanmean(toplot(:,i)));
    b(i).FaceColor = cols{i};
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.5;

    vs(i) = scatter(randn(numel(anms),1) * 0.1 + xs(i)*ones(size(toplot(:,i))),toplot(:,i),30,'MarkerFaceColor',cols{i},...
        'MarkerEdgeColor','k','LineWidth',1);

    errorbar(b(i).XEndPoints,nanmean(toplot(:,i)),nanstd(toplot(:,i)),'LineStyle','none','Color','k','LineWidth',1)

    xs_(:,i) = vs(i).XData';
    ys_(:,i) = toplot(:,i);
end

for i = 1:size(xs_,1)
    patchline(xs_(i,1:2),ys_(i,1:2),'EdgeAlpha',0.4)
    patchline(xs_(i,3:4),ys_(i,3:4),'EdgeAlpha',0.4)
end

ylabel("Time (s)")
ylim([0,1])
ax = gca;
ax.FontSize = 10;
ax.XTick = xs;
xticklabels([ "RT, 2AFC" "RT, AW" "L1, 2AFC" "L1, AW "])



%% distribution of rts

allrt = cell2mat(rt');
figure; 
histogram(allrt,200,'EdgeColor','none')
xline(nanmedian(allrt),'k--')
ax = gca;
ax.FontSize = 13;
xlabel('Reaction time (s)')
ylabel('# of trials')



%% LICK RATE

lickrate = lickRate(obj);

% group  by mouse, condition, only for trials we want
anms = unique({meta.anm});
cond2use = 1:2;
for ianm = 1:numel(anms) % for each animal
    [~,mask] = patternMatchCellArray({meta(:).anm}',{anms{ianm}},'all');
    sessix = find(mask);
    anm_rate_cell{ianm} = lickrate(sessix);

    anm_rate{ianm} = cell(1,numel(cond2use));
    for isess = 1:numel(sessix)
        for icond = 1:numel(cond2use)
            trix = params(sessix(isess)).trialid{cond2use(icond)};
            anm_rate{ianm}{icond} = [anm_rate{ianm}{icond} anm_rate_cell{ianm}{isess}(trix)];
        end
    end

end

% get means by condition for each mouse
for anmix = 1:numel(anm_rate)
    temp = anm_rate{anmix};
    for icond = 1:numel(temp)
        mean_rate_by_mouse_by_cond(anmix,icond) = nanmean(temp{icond}); % (anm,cond)
    end
end

%%
close all

clear xs_ ys_

toplot = (mean_rate_by_mouse_by_cond);

clrs = getColors();
cols{1} = clrs.afc;
cols{2} = clrs.aw;
cols{3} = clrs.afc;
cols{4} = clrs.aw;

f = figure;
ax = gca;
hold on;
f.Position = [680   739   270   239];

xs = 1:size(toplot,2);
for i = 1:numel(xs)
    b(i) = bar(xs(i),nanmean(toplot(:,i)));
    b(i).FaceColor = cols{i};
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.5;

    vs(i) = scatter(randn(numel(anms),1) * 0.1 + xs(i)*ones(size(toplot(:,i))),toplot(:,i),30,'MarkerFaceColor',cols{i},...
        'MarkerEdgeColor','k','LineWidth',1);

    errorbar(b(i).XEndPoints,nanmean(toplot(:,i)),nanstd(toplot(:,i)),'LineStyle','none','Color','k','LineWidth',1)

    xs_(:,i) = vs(i).XData';
    ys_(:,i) = toplot(:,i);
end

for i = 1:size(xs_,1)
    patchline(xs_(i,:),ys_(i,:),'EdgeAlpha',0.4)
end

ylabel("Lick Rate (Hz)")
% ylim([0,1])
ax = gca;
ax.FontSize = 10;
ax.XTick = xs;
xticklabels([ "2AFC" "AW" ])



