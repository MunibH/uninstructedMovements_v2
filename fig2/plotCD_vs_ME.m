clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
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
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off

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

% --- ALM --- 
% meta = loadJEB6_ALMVideo(meta,datapth);
% meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth);
% meta = loadEKH3_ALMVideo(meta,datapth);
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

meta = meta(3);

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

% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end

% ------------------------------------------
% -- Kinematics --
% kin (struct array) - one entry per session
% TODO - save kin to obj and load them
% ------------------------------------------
disp('Loading Kinematics')
for sessix = 1:numel(meta)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end

%% CODING DIMENSIONS

clearvars -except obj meta params sel_corr_mat me

cond2use = [2 3]; % left hit, right hit
cond2proj = [2 3 4 5];
proj_trialdat = 1;
rez = getCodingDimensions(obj,params,cond2use,cond2proj,proj_trialdat);

% allrez = concatRezAcrossSessions(rez);

%% PLOTS
close all

sav = 0; % 1=save, 0=no_save

% plot CD delay and avg motion energy for each session
cond2use = [2 3];
% plotSingleTrial_CDdelay_ME(obj,params,rez,me,meta,cond2use);
% plotCorrelation_Cddelay_ME(obj,params,rez,me,meta,cond2use);
plotCorrelation_Cddelay_ME_v2(obj,params,rez,me,meta,cond2use);

% % these plots are really for figure 1, works here as well
% plotmiss = 1;
% plotCDProj(allrez,rez,sav,plotmiss)
% plotCDVarExp(allrez,sav)
% plotSelectivity(allrez,rez,sav)
% plotSelectivityExplained(allrez,rez,sav)

%% cd vs me (tavg)

close all

cols = getColors();
clrs{1} = cols.lhit;
clrs{2} = cols.rhit;

cond2use = [2 3];
xlims = [-1 0.05];

sample = mode(obj.bp.ev.sample) - 2.5;
delay = mode(obj.bp.ev.delay) - 2.5;

figure; 
ax = gca; hold on;
for i = 1:numel(cond2use)
    trix = params.trialid{cond2use(i)};
    plotme = me.data(:,trix);
    shadedErrorBar(obj.time,mean(plotme,2),std(plotme,[],2)./sqrt(numel(trix)/2), {'Color',clrs{i},'LineWidth',2},0.5,ax);
end
ys = ax.YLim;
p = patch([-0.9 0 0 -0.9],[ys(1) ys(1) ys(2) ys(2)],'k');
p.EdgeColor = 'none';
% p.FaceColor = [120 120 120]./255;
p.FaceColor = [62, 168, 105]./255;
p.FaceAlpha = 0.1;
xline(sample,'k:','LineWidth',2)
xline(delay,'k:','LineWidth',2)
xline(0,'k:','LineWidth',2)
xlabel('Time (s) from go cue')
ylabel('Motion Energy')
xlim(xlims)
ylim([ax.YLim(1)+5 ax.YLim(2)-30])

% cd
figure; 
ax = gca; hold on;
for i = 1:numel(cond2use)
    trix = params.trialid{cond2use(i)};
    plotcd = rez.cd_proj_trialdat(:,trix);
    shadedErrorBar(obj.time,mean(plotcd,2),std(plotcd,[],2)./sqrt(numel(trix)/2), {'Color',clrs{i},'LineWidth',2},0.5,ax);
end
ys = ax.YLim;
p = patch([-0.9 0 0 -0.9],[ys(1) ys(1) ys(2) ys(2)],'k');
p.EdgeColor = 'none';
% p.FaceColor = [120 120 120]./255;
p.FaceColor = [62, 168, 105]./255;
p.FaceAlpha = 0.1;
xline(sample,'k:','LineWidth',2)
xline(delay,'k:','LineWidth',2)
xline(0,'k:','LineWidth',2)
xlabel('Time (s) from go cue')
ylabel('CD Delay')
xlim(xlims)
ylim([ax.YLim(1)+10 ax.YLim(2)-30])


