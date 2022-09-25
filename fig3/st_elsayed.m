clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')))

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
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off

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

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];
%
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);

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



%% Null and Potent Space

clearvars -except obj meta params me sav

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
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat);

    % -- null and potent spaces
    cond2use = [2 3 4 5]; % right hit, left hit, right miss, left miss
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use);

    % -- coding dimensions
    cond2use = [1 2]; % right hits, left hits (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
    cd_null(sessix) = getCodingDimensions(rez(sessix).N_null_psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat);
    cd_potent(sessix) = getCodingDimensions(rez(sessix).N_potent_psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat);
end



cd_null_all = concatRezAcrossSessions(cd_null);
cd_potent_all = concatRezAcrossSessions(cd_potent);


%% plots

close all

sav = 0;


% -----------------------------------------------------------------------
% -- Null and Potent Space Single Trial Projections --
% -----------------------------------------------------------------------

% % - projections showing move / quiet somehow
cond2use = [2 3]; % right hits, left hits
ndims = 4; % top ndims variance explaining dimensions
% plotSingleTrialNPHeatmaps(rez,params,mendims,cond2use,meta);


% % - how much variance in move and non-move time points
cond2use = [2 3]; % right hits, left hits
% plotVarianceInEpochs(rez,me,params,cond2use);

% % - ve
% plotVarianceExplained_NP(rez);

% % - ve over time (TODO)
% % % plotVarianceExplained_NP_overTime(rez);


% -----------------------------------------------------------------------
% -- Null and Potent Space Trial-averaged Projections --
% -----------------------------------------------------------------------

ndims = 5; % how many n/p dimensions to plot, in order of variance explained
cond2plot = [1 2]; % right hit, left hit
% plot_NP_PSTH(rez,obj,params,ndims,cond2plot,meta)


% -----------------------------------------------------------------------
% -- Coding Dimensions --
% -----------------------------------------------------------------------

titlestring = 'Null';
% plotCDProj(cd_null_all,cd_null,sav,titlestring)
% plotCDVarExp(cd_null_all,sav,titlestring)
% plotSelectivity(cd_null_all,cd_null,sav,titlestring)
% plotSelectivityExplained(cd_null_all,cd_null,sav,titlestring)

titlestring = 'Potent';
% plotCDProj(cd_potent_all,cd_potent,sav,titlestring)
% plotCDVarExp(cd_potent_all,sav,titlestring)
% plotSelectivity(cd_potent_all,cd_potent,sav,titlestring)
% plotSelectivityExplained(cd_potent_all,cd_potent,sav,titlestring)



%% t=0 is the go cue, but only on trials where the animals were not moving PRIOR to the go cue
% same plots as plotSelectivityExplained

%% t=0 is transitions between non-movement and movement that do not coincide with the go cue
% same plots as plotSelectivityExplained

stitch_dist = 0.1; % in seconds, stitch together movement bouts shorter than this
purge_dist = 0.5; % in seconds, remove move bouts shorter than this value, after stitching complete
tbout = 0.2; % move/non-move bout/transition required to be at least this long in seconds
[dat.mdat,dat.mdat_leg,dat.qdat,dat.qdat_leg,newme] = nonGC_moveTransitions(obj,me,params,stitch_dist,purge_dist,tbout);

%%
close all

dim = 1; % dim to plot (most to least variance explained)
% plot_nonGC_moveTransitions_singleTrials(dat,obj,newme,rez,params,dim,meta)

ndims = 10;
% plot_nonGC_moveTransitions_trialAvg(dat,obj,me,rez,params,meta,ndims) % separate tiles for each dimension
% plot_nonGC_moveTransitions_trialAvg_v2(dat,obj,me,rez,params,meta,ndims) % all dimensions plotted on same axis, 
plot_nonGC_moveTransitions_trialAvg_v3(dat,obj,me,rez,params,meta,ndims) % mean,stderr across dimensions


%% PROBABLY DON'T NEED TO DO THIS (will come back if needed)
% from movement and non-movement transitions, shrink each epoch to tbout
% length 
% e.g. moveboutstart = 0.2; moveboutend = 0.6; quietboutstart = 0.6;
% quietboutend = 0.9; 
% trim this to: moveboutstart = 0.4; moveboutend = 0.6; quietboutstart =
% 0.6; quietboutend = 0.8;

% for sessix = 1:numel(m2q)
%     for trix = 1:numel(m2q(sessix).moveStart)
%         % trim move-to-quiet transitions to size (each epoch tbout long)
%         mdat = m2q(sessix);
%         for i = 1:numel(mdat.moveStart{trix})
%             % for m2q, we will trim moveStart and quietEnd (if needed)
%             mstart = mdat.moveStart{trix}(i);
%             mend = mdat.moveEnd{trix}(i);
%             qstart = mdat.quietStart{trix}(i);
%             qend = mdat.quietEnd{trix}(i);
%             ml = mend - mstart;
%         end
%         
% 
%         % trim quiet-to-move transitions to size (each epoch tbout long)
%         qdat = q2m(sessix);
%     end
% end


%% plot motion energy 
% close all
% rng(pi)
% for sesh = 1:numel(meta)
%     temp = newme(sesh).data;
%     temp2 = newme(sesh).move;
%     z = temp;
%     z(~temp2) = nan;
%     y = temp;
%     y(temp2) = nan;
%     f = figure; f.Position = [589    65   663   927];
%     hold on;
%     offset = 0; dy = 55;
%     k = 30;
%     nTrials = size(z,2);
%     trix = randsample(nTrials,k,false);
%     for i = 1:k
%         plot(y(:,trix(i)) + offset,'k')
%         plot(z(:,trix(i)) + offset,'r')
%         offset = offset + dy;
%     end
% end
% 


% 
% figure; plot(newme(sessix).data(:,trix))
% hold on
% plot(newme(sessix).move(:,trix)*100)



























