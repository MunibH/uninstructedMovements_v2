% decode reaction time from delay period N/P CDChoice


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

% time warp params
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials (1)
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off (2)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off (3)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off (4)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off (5)
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
    cond2use = [2:5 15:18]; % right hit, left hit, right miss, left miss, 2afc and aw
    cond2proj = [8:11 14];
    nullalltime = 0; % use all time points to estimate null space if 1
    onlyAW = 0; % only use AW trials
    delayOnly = 0; % only use delay period
    responseOnly = 0; % only use response period
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj, nullalltime, onlyAW, delayOnly, responseOnly);

    % -- coding dimensions
    cond2use = [1 2]; % right hits, left hits (corresponding to null/potent psths in rez)
    cond2proj = [1:4]; % right hits, left hits, right miss, left miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
    rampcond = 5; % corresponding to cond2proj in null/potent analysis
    %     cd_null(sessix) = getCodingDimensions(rez(sessix).N_null_psth,rez(sessix).N_null,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,rampcond);
    %     cd_potent(sessix) = getCodingDimensions(rez(sessix).N_potent_psth,rez(sessix).N_potent,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,rampcond);

    cd_null(sessix) = getCodingDimensions(rez(sessix).recon_psth.null,trialdat_zscored,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj, rampcond);
    cd_potent(sessix) = getCodingDimensions(rez(sessix).recon_psth.potent,trialdat_zscored,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj, rampcond);
end
disp('DONE')


cd_null_all = concatRezAcrossSessions(cd_null);
cd_potent_all = concatRezAcrossSessions(cd_potent);

%% single trial projs onto cd choice

sm = 21;

for sessix = 1:numel(meta)
    Wchoice = cd_null(sessix).cd_mode_orth(:,1);
    trialdat.null{sessix} = tensorprod(rez(sessix).recon.null,Wchoice,3,1); % (time,trials)
    trialdat.null{sessix} = mySmooth(trialdat.null{sessix},sm,'reflect');

    Wchoice = cd_potent(sessix).cd_mode_orth(:,1);
    trialdat.potent{sessix} = tensorprod(rez(sessix).recon.potent,Wchoice,3,1); % (time,trials)
    trialdat.potent{sessix} = mySmooth(trialdat.potent{sessix},sm,'reflect');
end



%% DECODING PARAMETERS

% input data = neural data (time*trials,neurons)
% output data = kin data   (time*trials,kin feats)

par.pre=6; % time bins prior to output used for decoding
par.post=0; % time bins after output used for decoding
par.dt = params(1).dt; % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using. Set params.dt and pre/post appropriately for you analysis
par.post_s = par.post .* params(1).dt;

% these parameters above are important for decoding accuracy. for actual
% analyses (like a final analysis to be included in a publication),
% you should vary these parameters only if you have a validation
% set that you won't test on until these parameters are set. Otherwise,
% there's a high risk of overfitting

% cross val folds
par.nFolds = 4;

% data sets
par.train = 1; % fraction of trials (just using cross-val here, matlab's kFoldPredict uses held out data for testing)
par.test = 1 - par.train;

% trials
par.cond2use = [6 7];

par.regularize = 1; % if 0, linear regression. if 1, ridge regression

par.prctile = [20 50];

% time2use = [-0.3 0];
time2use = [2 2.48];
ix = findTimeIX(obj(1).time,time2use);
par.time2use = ix(1):ix(2);


%% DECODING

% rt = firstTongueRT(obj);
% rt = firstJawRT(obj);
% prevchoice = getPrevChoice(obj);
outcome = getOutcome(obj);

% ideas
% avg delay period decode rt
% when does rt become decodable (moving bins)

close all

for sessix = 1:numel(meta)
    disp([num2str(sessix) ' / ' num2str(numel(meta))])

    %%
    % reaction time
    clear thisrt x mask rtcat Y X xx mask rtcat nanmask
    thisrt = rt{sessix};
    x = prctile(thisrt,par.prctile);
    mask(:,1) = thisrt<x(1);
    mask(:,2) = thisrt>=x(1)&thisrt<x(2);
    mask(:,3) = thisrt>=x(2);
    thisrt(mask(:,1)) = -1;
    thisrt(mask(:,2)) = 0;
    thisrt(mask(:,3)) = 1;
    rtcat = cell(size(thisrt));
    rtcat(thisrt==-1) = {'fast'};
    rtcat(thisrt==0) = {'mid'};
    rtcat(thisrt==1) = {'slow'};

    thisoutcome = outcome{sessix};
    rtcat = cell(size(thisrt));
    rtcat(thisoutcome==1) = {'hit'};
    rtcat(thisoutcome==0) = {'miss'};
    rtcat(isnan(thisoutcome)) = {'none'};

    % % predict RT from single trial neural data in null and potent spaces
    % X = trialdat.null{sessix};
    % xx(:,1) = mean(X(:,mask(:,1)),2);
    % xx(:,2) = mean(X(:,mask(:,2)),2);
    % xx(:,3) = mean(X(:,mask(:,3)),2);
    % f = figure;
    % hold on;
    % plot(obj(1).time,xx(:,1)); plot(obj(1).time,xx(:,2)); plot(obj(1).time,xx(:,3));
    % X = trialdat.potent{sessix};
    % xx(:,1) = mean(X(:,mask(:,1)),2);
    % xx(:,2) = mean(X(:,mask(:,2)),2);
    % xx(:,3) = mean(X(:,mask(:,3)),2);
    % f = figure;
    % hold on;
    % plot(obj(1).time,xx(:,1)); plot(obj(1).time,xx(:,2)); plot(obj(1).time,xx(:,3));

    %% Null
    % nanmask = isnan(thisrt);
    nanmask = false(size(thisoutcome));

    X.train = squeeze(mean(rez(sessix).N_null(par.time2use,~nanmask,:),1)); % (trials,dims)
    X.size = size(X.train);
    Y.train = rtcat(~nanmask)';
    t = templateSVM('Standardize',true);
    mdl = fitcecoc(X.train,Y.train,'Learners',t,...
        'ClassNames',unique(rtcat));
    cvmdl = crossval(mdl);
    genError.null(sessix) = kfoldLoss(cvmdl);

    X.train = squeeze(mean(rez(sessix).N_potent(par.time2use,~nanmask,:),1)); % (trials,dims)
    X.size = size(X.train);
    t = templateSVM('Standardize',true);
    mdl = fitcecoc(X.train,Y.train,'Learners',t,...
        'ClassNames',unique(rtcat));
    cvmdl = crossval(mdl);
    genError.potent(sessix) = kfoldLoss(cvmdl);

    % X.train = trialdat.null{sessix}(par.time2use,~nanmask); % trials    
    % X.size = size(X.train);
    % Y.train = rtcat';
    % Y.train = Y.train(~nanmask);
    % newY = repmat(Y.train,X.size(1),1);
    % X.train = reshape(X.train, size(X.train,1)*size(X.train,2), 1);
    % X.train = reshapePredictors(X.train,par);
    % Y.train = reshape(newY,size(newY,1)*size(newY,2),1);
    % t = templateSVM('Standardize',true);
    % mdl = fitcecoc(X.train,Y.train,'Learners',t,...
    %     'ClassNames',{'fast','mid','slow'});
    % cvmdl = crossval(mdl);
    % genError.null(sessix) = kfoldLoss(cvmdl);

    % % Potent
    % X.train = trialdat.potent{sessix}(par.time2use,~nanmask); % trials    
    % X.size = size(X.train);
    %  X.train = reshape(X.train, size(X.train,1)*size(X.train,2), 1);
    % X.train = reshapePredictors(X.train,par);
    % t = templateSVM('Standardize',true);
    % mdl = fitcecoc(X.train,Y.train,'Learners',t,...
    %     'ClassNames',{'fast','mid','slow'});
    % cvmdl = crossval(mdl);
    % genError.potent(sessix) = kfoldLoss(cvmdl);

  

end

%%
f = figure;
f.Position = [748   394   193   272];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

cols = getColors;
col{1} = cols.null;
col{2} = cols.potent;

xs = [1 2];
fns = fieldnames(genError);
for i = 1:numel(xs)
    dat = (1 - genError.(fns{i})) * 100;
    b(i) = bar(xs(i),mean(dat,'omitmissing'));
    b(i).FaceColor = col{i};
    b(i).EdgeColor = 'none';
    % b(i).FaceAlpha = 0.5;
    b(i).BarWidth = 0.7;
    vs(i) = scatter(xs(i)*ones(numel(dat),1),dat,20,'MarkerFaceColor','k',...
        'MarkerEdgeColor','w','XJitter','randn','XJitterWidth',0.25);
    errorbar(b(i).XEndPoints,nanmean(dat),nanstd(dat),'LineStyle','none','Color','k','LineWidth',1)
end
ax.XTick = [];
ylabel("Accuracy")
ax = gca;
% ax.FontSize = 12;
yline(33.3,'k-')
title('pred. outcome from NP single trials')







