clear, clc, close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'figNP/')))
rmpath(genpath(fullfile(utilspth,'fig1/')))
rmpath(genpath(fullfile(utilspth,'musall2019/')))


% add paths for figure specific functions
addpath(genpath(pwd))

clc

%% TODO
% - created new objs from pipeline and nidq data - only use those
% - need to handle frameTimes appropriately now (subtract 0.5 secs)
% - for choice decoding, just use avg during stim period

%% default params (same for each session)

dfparams = [];

% -- time alignment params --
% dfparams.alignEv = 'goCue';
% dfparams.times = [-2.5 2.5]; % relative to goCue
dfparams.alignEv = 'delay';
dfparams.times = [-1.0 2.5]; % relative to delay

dfparams.dt_vid = 0.0025;
dfparams.time = dfparams.times(1):dfparams.dt_vid:dfparams.times(2);

dfparams.warp = 0; % 0 means no warping, 1 means warp delay period to 


% -- trial type params --
dfparams.cond(1) = {'(hit|miss|no)&~stim.enable&~autowater&~autolearn&~early'}; % all trials, no stim, no autowater, no autolearn
dfparams.cond(end+1) = {'(hit|miss|no)&stim.enable&~autowater&~autolearn&~early'};  % all trials trials, stim, no autowater, no autolearn
dfparams.cond(end+1) = {'R&~stim.enable&~autowater&~autolearn&~early'}; % right trials, no stim, no autowater
dfparams.cond(end+1) = {'R&stim.enable&~autowater&~autolearn&~early'};  % right trials, stim, no autowater
dfparams.cond(end+1) = {'L&~stim.enable&~autowater&~autolearn&~early'}; % left trials, no stim, no autowater
dfparams.cond(end+1) = {'L&stim.enable&~autowater&~autolearn&~early'};  % left trials, stim, no autowater
dfparams.cond(end+1) = {'R&hit&~stim.enable&~autowater&~autolearn&~early'}; % right hit trials, no stim, no autowater
dfparams.cond(end+1) = {'R&hit&stim.enable&~autowater&~autolearn&~early'};  % right hit trials, stim, no autowater
dfparams.cond(end+1) = {'L&hit&~stim.enable&~autowater&~autolearn&~early'}; % left hit trials, no stim, no autowater
dfparams.cond(end+1) = {'L&hit&stim.enable&~autowater&~autolearn&~early'};  % left hit trials, stim, no autowater
% dfparams.cond(end+1) = {'R&miss&~stim.enable&~autowater&~autolearn'}; % left hit trials, no stim, no autowater
% dfparams.cond(end+1) = {'L&miss&stim.enable&~autowater&~autolearn'};  % left hit trials, stim, no autowater

% -- stim types --
dfparams.stim.types = {'Bi_MC','Right_MC','Left_MC','Bi_ALM','Bi_M1TJ','Right_ALM','Right_M1TJ','Left_ALM','Left_M1TJ'}; 
% dfparams.stim.num   = logical([1 1 1 1 1 1 1 1 1]);   % ALL
% dfparams.stim.num   = logical([0 0 0 0 0 1 1 1 1]);   % Right_ALM / Left_ALM / Right_M1TJ / Left_M1TJ
% dfparams.stim.num   = logical([0 0 0 1 0 0 0 0 0]);   % Bi_ALM
% dfparams.stim.num   = logical([0 0 0 0 1 0 0 0 0]);   % Bi_M1TJ
dfparams.stim.num   = logical([1 0 0 0 0 0 0 0 0]);   % Bi_MC
% dfparams.stim.num   = logical([0 0 0 1 1 0 0 0 0]);   % Bi_M1TJ Bi_ALM
% dfparams.stim.num   = logical([0 1 1 0 0 0 0 0 0]);   % Right_MC
% dfparams.stim.num   = logical([0 0 1 0 0 0 0 0 0]);   % Left_MC
% dfparams.stim.num   = logical([0 0 0 0 0 1 0 0 0]);   % Right_ALM
% dfparams.stim.num   = logical([0 0 0 0 0 0 0 1 0]);   % Left_ALM
% dfparams.stim.num   = logical([0 0 0 0 0 0 1 0 0]);   % Right_M1TJ
% dfparams.stim.num   = logical([0 0 0 0 0 0 0 0 1]);   % Left_M1TJ

dfparams.stim.pow.types = [1.08, 3.14, 8]; % mW, 3.14 for all unilateral sessions and most bilateral sessions
dfparams.stim.pow.num = logical([1 1 1]);


% -- plotting params --
dfparams.plt.color{1}     = [10, 10, 10];
dfparams.plt.color{end+1} = [120, 117, 117];
dfparams.plt.color{end+1} = [31, 23, 252];
dfparams.plt.color{end+1} = [22, 172, 247];
dfparams.plt.color{end+1} = [252, 23, 23];
dfparams.plt.color{end+1} = [252, 23, 130];
dfparams.plt.color = cellfun(@(x) x./255, dfparams.plt.color, 'UniformOutput', false);
dfparams.plt.ms = {'.','.','x','x','o','o'};



%% load data objects

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];
meta = loadMAH13_MCStim(meta,datapth);
meta = loadMAH14_MCStim(meta,datapth);
% meta = loadMAH17_MCStim(meta,datapth);

% subset based on stim types
stim2use = dfparams.stim.types(dfparams.stim.num);
use = false(size(meta));
for sessix = 1:numel(use)
    [~,mask] = patternMatchCellArray({meta(sessix).stimLoc}, stim2use,'any');
    if mask
        use(sessix) = true;
    end
end
meta = meta(use);

% subset based on stim pow
pow2use = dfparams.stim.pow.types(dfparams.stim.pow.num);
pows = [meta(:).stimPow];
use = ismember(pows,pow2use);
meta = meta(use);

disp(' ')
disp(['Loading ' num2str(numel(meta)) ' sessions with the following parameters:'])
disp(['Animals: ' unique({meta.anm})])
disp(['Stim Loc: ' stim2use])
disp(['    Stim Power: ' num2str(pow2use)])
disp('')


obj = loadObjs(meta);


%% find trials for each condition

for i = 1:numel(meta)
    params(i).trialid = findTrials(obj(i), dfparams.cond);
end

% %%
% if dfparams.warp
%     obj = warpData(obj,params);
% end

%% behavioral performance
close all

rez = getPerformance(meta,obj,params); % rez is struct array, each entyr is an animal. perf is (sessions,conditions)

% plots
cond2use = 1:6;
connectConds = 0;
plotPerformanceAllMice(meta,obj,rez,dfparams,params,cond2use,connectConds)
% plotPerformanceEachMouse(meta,obj,rez,dfparams,params) % TODO


%% number of early licks per trial

% % TODO: classify early licks by kinematics, not bpod data
% 
% cond2use = 3:6;
% early_per_trial = nan(numel(obj),numel(cond2use)); % (sessions,conds)
% for sessix = 1:numel(obj)
%     for condix = 1:numel(cond2use)
%         trials2use = params(sessix).trialid{cond2use(condix)};
%         nTrials = numel(trials2use);
% 
%         nEarly = sum(obj(sessix).bp.early(trials2use));
%         early_per_trial(sessix,condix) = nEarly / nTrials;
%     end
% end
% 
% f = figure; 
% f.Position = [316          73        1232         905];
% ax = axes(f);
% nCond = numel(cond2use);
% violincols = reshape(cell2mat(dfparams.plt.color(cond2use)),3,nCond)';
% vs = violinplot(early_per_trial,dfparams.cond(cond2use),...
%     'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1}, 'ViolinColor', violincols);
% ylabel('fraction of early licks per trial')
% ylim([0,1])
% title('early licks all sessions, all mice')
% ax.FontSize = 20;

%% kinematics

for sessix = 1:numel(obj)
    if ~isstruct(obj(sessix).me)
        temp = obj(sessix).me;
        obj(sessix).me = [];
        obj(sessix).me.data = temp;
        clear temp
    end
end
% kin.dat (don't use)
% kin.featLeg corresponds to 3rd dimension of kinfeats/kinfeats_norm
[kin,kinfeats,kinfeats_norm] = getKin(meta,obj,dfparams,params);

%% plot kinematics
close all
% feats2plot = {'tongue_ydisp_view1',...
%               'tongue_yvel_view1',...
%               'jaw_ydisp_view1',...
%               'jaw_yvel_view1'};

% feats2plot = {'tongue_ydisp_view1',...
%     'jaw_ydisp_view2',...
%     'jaw_yvel_view2'};
% feats2plot = {'tongue_ydisp_view1','jaw_ydisp_view1','motion_energy'};
feats2plot = {'motion_energy'};
cond2plot = 1:2;
% cond2plot = 3:6;
% cond2plot = 7:10;
sav = 0;

plotKinfeats(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)



%% avg jaw velocity during stim
close all
if strcmpi(dfparams.alignEv,'delay') % function depends on data aligned to delay period, since we use the stim period to measure avgjawvel
%     feats2plot = {'jaw_ydisp_view1',...
%         'jaw_yvel_view1',...
%         'motion_energy'};
%     feats2plot = {'jaw_yvel_view1',...
%         'tongue_ydisp_view1'};
    feats2plot = {'motion_energy'};
    cond2plot = 3:6;
    sav = 0;

    plotAvgFeatValDuringStim_v2(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)

%     plotAvgFeatValDuringStim(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)
%     plotAvgFeatValDuringStim_singleTrials(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)
end








