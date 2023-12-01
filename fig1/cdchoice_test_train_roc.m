clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig3/')))
rmpath(genpath(fullfile(utilspth,'fig2/')))
rmpath(genpath(fullfile(utilspth,'figx/')))
rmpath(genpath(fullfile(utilspth,'mc_stim/')))
rmpath(genpath(fullfile(utilspth,'MotionMapper/')))
rmpath(genpath(fullfile(utilspth,'musall2919/')))

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
params.condition(1)     = {'(hit|miss|no)'};                                              % all trials       (1)
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};      % right hits, 2afc (2)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};      % left hit, 2afc   (3)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};     % right miss, 2afc (4)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};     % left miss, 2afc  (5)
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};        % 2afc hits        (6)
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};         % aw hits          (7)
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};       % right hits, aw   (8)
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};       % left hits, aw    (9)
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};        % ramping          (10)

% params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % left hits, aw    (9)


params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;
params.bctype = 'reflect'; % reflect, zeropad, none

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance


params.advance_movement = 0.0;


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
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);


params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

%% PARAMETERS

p.train = 0.3; % trial percentage
p.test = 1 - p.train;
p.cond2use = [2 3]; % left hit , right hit


p.epoch = 'goCue';
ix = findTimeIX(obj(1).time,[-0.4 0]);
p.ix = ix(1):ix(2);

%% CODING DIMENSIONS

for isess = 1:numel(meta)

    % sample equal number of trials across cond2use
    trials = params(isess).trialid(p.cond2use);
    nTrials = min(cell2mat(cellfun(@(x) numel(x), trials, 'uni', 0)));
    trials = cellfun(@(x) randsample(x,nTrials,false), trials, 'uni', 0);

    % test train split
    cv = cvpartition(nTrials,'HoldOut',p.test);
    idx = cv.test;
    t.train = cellfun(@(x) x(~idx), trials, 'uni', 0);
    t.test = cellfun(@(x) x(idx), trials, 'uni', 0);

    % predict choice (resp)
    clear y
    for i = 1:numel(trials)
        y.train{i} = i * ones(size(t.train{i}));
        y.test{i} = i * ones(size(t.test{i}));
    end
    y.train = cell2mat(y.train'); y.test = cell2mat(y.test');
    y.train(y.train==1) = 0; y.test(y.test==1) = 0;
    y.train(y.train==2) = 1; y.test(y.test==2) = 1;


    % get data
    for i = 1:numel(trials)
        dat.train{i} = permute(obj(isess).trialdat(:,:,t.train{i}),[1 3 2]); % (time,trials,neurons)
        datmu.train{i} = squeeze(nanmean(dat.train{i},2)); % (time,neurons)

        dat.test{i} = permute(obj(isess).trialdat(:,:,t.test{i}),[1 3 2]); % (time,trials,neurons)
        datmu.test{i} = squeeze(nanmean(dat.train{i},2)); % (time,neurons)
    end


    % get coding dimension (choice)
    % find time points to use
    psth = (cat(3,datmu.train{1},datmu.train{2})); % (time,neurons,cond)
    mode = calcCD(psth,p.ix,[1 2]); % (neurons,1)

    % proj single trials
    trialdat = cat(2,dat.test{1},dat.test{2}); % (time,trials,neurons)
    proj = tensorprod(trialdat,mode,3,1); % (time,trials)
    % figure; imagesc(proj')


    % classifier
    X = mean(proj(p.ix,:),1); % (trials,1)
    mdl = fitglm(X,y.test,'Distribution','binomial','Link','logit');

    scores = mdl.Fitted.Probability;
    [Xperf{isess},Yperf{isess},~,AUC(isess)] = perfcurve(y.test,scores,1);
end

%% plot
close all

load('session_sortix_roc_cdchoice.mat') % comes from codingDirections.m -> sortix of late delay selectivity

Xplot = Xperf(sortix);
Yplot = Yperf(sortix);
AUCplot = AUC(sortix);

c1 = [0 0 0]./255;
c2 = [0.8 0.8 0.8];
cmap_ = flip(createcolormap(numel(meta), c1, c2));

f = figure;
f.Position = [680   582   337   296];
f.Renderer = 'painters';
ax = prettifyPlot(gca);
hold on;
for isess = 1:numel(meta)
    plot(Xplot{isess},Yplot{isess},'Color',cmap_(isess,:),'LineWidth',1.5)
end
xlabel('False positive rate')
ylabel('True positive rate')
ln = line([0 1],[0 1]);
ln.LineWidth = 1;
ln.Color = 'k';

vio = getXCoordsForViolin(AUC, []);
xx = (vio.jitter.*vio.jitterStrength);
ax = prettifyPlot(axes('Position',[.6 .2 .2 .3]));
hold on;
b = bar(0,mean(AUCplot));
b.EdgeColor = 'none';
b.FaceColor = [76 0 153] ./ 255;
b.BarWidth = 0.5;
scatter(xx,AUCplot,15,cmap_,'filled','markeredgecolor','w')
xlabel(ax,'');
xticks(ax,[])
ylabel(ax,'AUC')
ax.FontSize = 11;

