clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')));
rmpath(genpath(fullfile(utilspth,'fig2/')));
rmpath(genpath(fullfile(utilspth,'mc_stim/')));

% add paths for figure specific functions
addpath(genpath(pwd))

clc

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off
params.condition(end+1) = {'R&no&~stim.enable&~autowater'};            % no right, no stim, aw off
params.condition(end+1) = {'L&no&~stim.enable&~autowater'};            % no left, no stim, aw off

params.tmin = -3.0;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
% params.quality = {'Excellent','Great','Good','Fair','Multi'};

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;


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


% --- M1TJ ---
% meta = loadJEB13_M1TJVideo(meta,datapth);

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
    %     kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end



%% Null and Potent Space

clearvars -except obj meta params me sav datapth kin

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -- Calculate coding directions from null and potent spaces --
% -----------------------------------------------------------------------

for sessix = 1:numel(meta)

    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat,obj(sessix));

    % -- null and potent spaces
    cond2use = [2 3 4 5]; % right hit, left hit, right miss, left miss
    cond2proj = [2:7];
    nullalltime = 0; % use all time points to estimate null space if 1
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj, nullalltime);

    % -- coding dimensions
    cond2use = [1 2]; % right hits, left hits (corresponding to null/potent psths in rez)
    cond2proj = [1:6]; % right hits, left hits, right miss, left miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
    cd_null(sessix) = getCodingDimensions(rez(sessix).N_null_psth,rez(sessix).N_null,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);
    cd_potent(sessix) = getCodingDimensions(rez(sessix).N_potent_psth,rez(sessix).N_potent,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);

end

cd_null_all = concatRezAcrossSessions(cd_null);
cd_potent_all = concatRezAcrossSessions(cd_potent);

%% plot coding dimensions from activity N/P spaces
close all

sav = 0;

plotmiss = 0;
plotno = 0;

titlestring = 'Null';
axnull = plotCDProj(cd_null_all,cd_null,sav,titlestring,plotmiss,plotno);
% plotCDVarExp(cd_null_all,sav,titlestring)
% plotSelectivity(cd_null_all,cd_null,sav,titlestring)
% plotSelectivityExplained(cd_null_all,cd_null,sav,titlestring)

titlestring = 'Potent';
axpotent = plotCDProj(cd_potent_all,cd_potent,sav,titlestring,plotmiss,plotno);
% plotCDVarExp(cd_potent_all,sav,titlestring)
% plotSelectivity(cd_potent_all,cd_potent,sav,titlestring)
% plotSelectivityExplained(cd_potent_all,cd_potent,sav,titlestring)

% match y axis limits across null and potent space CDs
for i = 1:numel(axnull)
    ys = [ min(axnull(i).YLim(1),axpotent(i).YLim(1)) max(axnull(i).YLim(2),axpotent(i).YLim(2)) ];
    axnull(i).YLim = ys;
    axpotent(i).YLim = ys;
end

titlestring = 'Null | Potent CDs';
% plotCDProj_NP(cd_potent_all,cd_null_all,cd_potent,cd_null,sav,titlestring,plotmiss)
% plotSelectivityExplained_NP(cd_potent_all,cd_null_all,cd_potent,cd_null,sav,titlestring)

%% variance explained

plotVarianceExplained_NP(rez);

plotVarianceExplained_DelayResponse(rez);

%% single trial magnitude activity in N/P
close all

cond2use = 2:3;
% plotME_NPMagnitude_singleTrials(meta,obj,params,me,rez,cond2use)

corr_ME_NPMagnitude(meta,obj,params,me,rez);

%% trial-averaged magnitude of activity in N/P
close all

cond2use = [2 3];
plotSubpsaceMagnitude_TrialAvg(obj,meta,params,rez,cond2use)

% plotSubpsaceMagnitude_TrialAvg_byRT(obj,meta,params,rez,cond2use)

%% session-averaged magnitude of activity in N/P
close all

cond2use = [2 3];
plotSubpsaceMagnitude_SessionAvg(obj,meta,params,rez,me,cond2use)

% plotSubpsaceMagnitude_TrialAvg_byRT(obj,meta,params,rez,cond2use)

%% magnitude of activity in N/P aligned to movement bouts


stitch_dist = 0.05; % in seconds, stitch together movement bouts shorter than this % 0.025
purge_dist = 0.1; % in seconds, remove move bouts shorter than this value, after stitching complete % 0.1
tbout = 0.3; % move/non-move bout/transition required to be at least this long in seconds % 0.3
[dat.mdat,dat.mdat_leg,dat.qdat,dat.qdat_leg,newme] = nonGC_moveTransitions(obj,me,params,stitch_dist,purge_dist,tbout);

% plot

close all

ndims = 10;
sav = 0;
% plot_nonGC_moveTransitions_trialAvg_v6(dat,obj,newme,rez,params,meta,ndims) % sumsqmag across dimensions, move to quiet
plot_nonGC_moveTransitions_trialAvg_v7(dat,obj,newme,rez,params,meta,ndims) % sumsqmag across dimensions, quiet to move


% plotSubspaceActivityAlignedToMoveBoutEnd(dat,obj,me,rez,params,meta)



%% plot NP projections (trial-averaged and single-trials)

close all

sav = 0;
% -----------------------------------------------------------------------
% -- Null and Potent Space Single Trial Projections --
% -----------------------------------------------------------------------

% % - projections showing move / quiet somehow
cond2use = [2 3]; % right hits, left hits
ndims = 4; % top ndims variance explaining dimensions
% plotSingleTrialNPHeatmaps(rez,params,me,ndims,cond2use,meta);

% -----------------------------------------------------------------------
% -- Null and Potent Space Trial-averaged Projections --
% -----------------------------------------------------------------------

ndims = 5; % how many n/p dimensions to plot, in order of variance explained
cond2plot = [1 2]; % right hit, left hit
% plot_NP_PSTH(rez,obj,params,ndims,cond2plot,meta)

%% PCA on trial-avg N/P activity

close all

cond2use = 1:2; % right/left hits, corresponding to rez.N_potent/null_psth
pcaNP(meta,obj,params,rez,cond2use)

% pause, press any key to cycle through sessions


%% selectivity in N/P activity
close all

cond2use_ta = [1 2]; % right and left hits, corresponding to trial-avg projs onto n/p
cond2use_st = [2 3]; % right and left hits, corresponding to single-trial projs onto n/p
plotSelectivityNP(meta,obj,params,rez,cond2use_ta,cond2use_st)


















