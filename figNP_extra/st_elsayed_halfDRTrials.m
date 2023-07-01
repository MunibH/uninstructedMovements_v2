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

% use half DR trials to estimate potent space
%   reconstruct all neural activity
% project other half DR on potent space
%   reconstruct all neural activity
% project WC trials on potent space
%   reconstruct all neural activity
%
% compare reconstructions

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warp params
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials (1)
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off (2)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off (3)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off (4)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off (5)
params.condition(end+1) = {'R&no&~stim.enable&~autowater&~early'};            % no right, no stim, aw off (6)
params.condition(end+1) = {'L&no&~stim.enable&~autowater&~early'};            % no left, no stim, aw off (7)

% for projections
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off (8)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off (9)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off (10)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off (11)

params.condition(end+1) = {'R&~stim.enable&~autowater&~early'}; % (12)
params.condition(end+1) = {'L&~stim.enable&~autowater&~early'}; % (13)

% for ramping
params.condition(end+1) = {'hit&~stim.enable&~autowater'};               % all hits, no stim, aw off (14)

% autowater
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};             % right hits, no stim, aw off (15)
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};             % left hits, no stim, aw off  (16)
params.condition(end+1) = {'R&miss&~stim.enable&autowater&~early'};             % right hits, no stim, aw off (17)
params.condition(end+1) = {'L&miss&~stim.enable&autowater&~early'};             % left hits, no stim, aw off  (18)

params.tmin = -2.4;
params.tmax = 2.5;
params.dt = 1/50;

% smooth with causal gaussian kernel
params.smooth = 15;
params.bctype = 'reflect';

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
% params.quality = {'Excellent','Great','Good','Fair','Multi'};
% params.quality = {'Excellent','Great','Good','Fair'};

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;


%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\munib\Documents\Economo-Lab\data';

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth); % selectivity in ME
% % % % meta = loadEKH1_ALMVideo(meta,datapth); % selectivity in ME
meta = loadEKH3_ALMVideo(meta,datapth); % selectivity in ME
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth); % selectivity in ME % go cue is at 2.3 instead of 2.5 like all other sessions??
meta = loadJEB15_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);


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
    disp(['Loading kinematics - Session ' num2str(sessix) '/' num2str(numel(meta))])
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
    % kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end



%% Null and Potent Space

clearvars -except obj meta params me sav datapth kin rt

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -- Calculate coding directions from null and potent spaces --
% -----------------------------------------------------------------------
disp('finding null and potent spaces')
for sessix = 1:numel(meta)
    disp([meta(sessix).anm ' ' meta(sessix).date])
    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix));

    % -- null and potent spaces
    cond2use = [2:5]; % right hit, left hit, right miss, left miss, 2afc and aw
    cond2proj = [8:11 14];
    nullalltime = 0; % use all time points to estimate null space if 1
    onlyAW = 0; % only use AW trials
    delayOnly = 0; % only use delay period
    responseOnly = 0; % only use response period
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj, nullalltime, onlyAW, delayOnly, responseOnly);
end
disp('DONE')



%% compare reconstructions from dr1,dr2, and wc trials

fns = fieldnames(rez(1).trials);
np = {'null','potent'};
for i = [1:7 19:numel(meta)]%1:numel(meta) % aw trials only
    for j = 1:numel(fns)
        trix.(fns{j}) = rez(i).trials.(fns{j});
    end
    trix.wc = cell2mat(params(i).trialid(15:18)');
    
    % reconstruct null
    tfns = {'dr1','dr2','wc'};
    for k = 1:numel(tfns)
        t = trix.(tfns{k});

        % reconstruct neural activity from null and potent spaces
        for npix = 1:numel(np)
            Q = rez(i).(['Q' np{npix}]);
            proj.(np{npix}).(tfns{k}) = rez(i).(['N_' np{npix}])(:,t,:); % (time,trials,dims)
            recon.(np{npix}).(tfns{k}) = tensorprod(proj.(np{npix}).(tfns{k}),Q,3,2); % (time, trials, neurons)
        end        
    end

    % compare reconstructions
    dr1 = squeeze(mean(recon.potent.dr1,2));
    dr2 = squeeze(mean(recon.potent.dr2,2));
    wc  = squeeze(mean(recon.potent.wc,2));

    % -- dr1 to dr2
    

    % -- dr1 to wc


    break
end






