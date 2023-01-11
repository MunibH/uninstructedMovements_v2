clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')));
rmpath(genpath(fullfile(utilspth,'mc_stim/')));

% add paths for figure specific functions
addpath(genpath(pwd))

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

params.tmin = -2.5;
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
end



%% Null and Potent Space

clearvars -except obj meta params me sav datapth


for sessix = 1:numel(meta)

    rez(sessix).psth = obj(sessix).psth(:,:,[2,3]); % right, left hits, 2afc
    rez(sessix).time = obj(sessix).time;

    % soft-normalize
    rez(sessix).softnorm_lambda = 0.01;
    % firing rate of a neuron, x of size (time,1), is transformed as:
    % x_norm = x / (lambda + max(x) - min(x))
    snpsth = rez(sessix).psth ./ (rez(sessix).softnorm_lambda + max(rez(sessix).psth) - min(rez(sessix).psth));

    % mean center across conditions
    for clu = 1:size(snpsth,2)
        % find mean at each time point for each condition (time,1)
        mean_across_cond = mean(snpsth(:,clu,:),3);
        snpsth(:,clu,:) = snpsth(:,clu,:) - repmat(mean_across_cond,[1,1,size(snpsth,3)]);
    end
    rez(sessix).psth_processed = snpsth; % soft normed and mean centered

    % epochs
    % prep
    edges = [-1.2 -0.2];
    [~,e1] = min(abs(rez(sessix).time - edges(1)));
    [~,e2] = min(abs(rez(sessix).time - edges(2)));
    rez(sessix).prepix = e1:e2;

    % move
    edges = [0.02 1.2];
    [~,e1] = min(abs(rez(sessix).time - edges(1)));
    [~,e2] = min(abs(rez(sessix).time - edges(2)));
    rez(sessix).moveix = e1:e2;

    % concatenates psth to (ct x n)
    psthprep = rez(sessix).psth_processed(rez(sessix).prepix,:,1);
    psthmove = rez(sessix).psth_processed(rez(sessix).moveix,:,1);
    for i = 2:size(rez(sessix).psth_processed,3)
        psthprep = [psthprep ; rez(sessix).psth_processed(rez(sessix).prepix,:,i)];
        psthmove = [psthmove ; rez(sessix).psth_processed(rez(sessix).moveix,:,i)];
    end

    % computes covariance of each epoch
    rez(sessix).Cprep = cov(psthprep);
    rez(sessix).Cmove = cov(psthmove);

    [~,~,explained] = myPCA(psthprep);
    rez(sessix).dPrep = numComponentsToExplainVariance(explained, 90);
    if rez(sessix).dPrep == 1 % always at least 2 dimensions
        rez(sessix).dPrep = 2;
    end


    [~,~,explained] = myPCA(psthmove);
    rez(sessix).dMove = numComponentsToExplainVariance(explained, 90);
    if rez(sessix).dMove == 1
        rez(sessix).dMove = 2;
    end


    rez(sessix).dMax = max(rez(sessix).dMove,rez(sessix).dPrep);

    % main optimization step
    rez(sessix).alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
    [Q,~,P,~,~] = orthogonal_subspaces(rez(sessix).Cmove,rez(sessix).dMove, ...
        rez(sessix).Cprep,rez(sessix).dPrep,rez(sessix).alpha);


    rez(sessix).Qpotent = Q*P{1};
    rez(sessix).Qnull = Q*P{2};

    % variance explained

    prepeigs = sort(eig(rez(sessix).Cprep),'descend');
    moveeigs = sort(eig(rez(sessix).Cmove),'descend');

    prepproj = rez(sessix).Qnull'*rez(sessix).Cprep*rez(sessix).Qnull;
    moveproj = rez(sessix).Qpotent'*rez(sessix).Cmove*rez(sessix).Qpotent;

    crossprepproj = rez(sessix).Qpotent'*rez(sessix).Cprep*rez(sessix).Qpotent;
    crossmoveproj = rez(sessix).Qnull'*rez(sessix).Cmove*rez(sessix).Qnull;

    rez(sessix).Qprep_null_ve = trace(prepproj) / sum(prepeigs(1:rez(sessix).dPrep)) * 100;
    rez(sessix).Qmove_potent_ve = trace(moveproj) / sum(moveeigs(1:rez(sessix).dMove)) * 100;

    rez(sessix).Qprep_potent_ve = trace(crossprepproj) / sum(prepeigs(1:rez(sessix).dPrep)) * 100;
    rez(sessix).Qmove_null_ve = trace(crossmoveproj) / sum(moveeigs(1:rez(sessix).dMove)) * 100;

    % projections

    Q = rez(sessix).Qnull;

    rez(sessix).proj_null = zeros(size(rez(sessix).psth,1),size(Q,2),size(rez(sessix).psth,3)); % (time,dims,conditions)
    for i = 1:2 % condition
        rez(sessix).proj_null(:,:,i) = rez(sessix).psth_processed(:,:,i) * Q;
    end

    Q = rez(sessix).Qpotent;

    rez(sessix).proj_potent = zeros(size(rez(sessix).psth,1),size(Q,2),size(rez(sessix).psth,3)); % (time,dims,conditions)
    for i = 1:2 % condition
        rez(sessix).proj_potent(:,:,i) = rez(sessix).psth_processed(:,:,i) * Q;
    end


    %single trial projections

    dat{1} = obj(sessix).trialdat(:,:,params(sessix).trialid{2});
    dat{2} = obj(sessix).trialdat(:,:,params(sessix).trialid{3});

    projdat = cellfun(@(x) permute(x, [1 3 2]),dat,'UniformOutput',false);
    projdat_reshape = cellfun(@(x) reshape(x,size(x,1)*size(x,2),size(x,3)), projdat, 'UniformOutput',false);


    for i = 1:numel(dat)
        projdat_null{i} = projdat_reshape{i} * rez(sessix).Qnull;
        projdat_potent{i} = projdat_reshape{i} * rez(sessix).Qpotent;
    end

    for i = 1:2
        projdat_null{i} = reshape(projdat_null{i}, numel(rez(sessix).time), numel(params(sessix).trialid{i+1}), rez(sessix).dPrep);
        projdat_potent{i} = reshape(projdat_potent{i}, numel(rez(sessix).time), numel(params(sessix).trialid{i+1}), rez(sessix).dMove);
    end

    null_trialdat = cat(2,projdat_null{1}, projdat_null{2});
    potent_trialdat = cat(2,projdat_potent{1}, projdat_potent{2});

    rez(sessix).null_ssm = sum(null_trialdat.^2,3);
    rez(sessix).potent_ssm = sum(potent_trialdat.^2,3);

end


%% plots


for sessix = 1:numel(rez)
    ai(1,sessix) = rez(sessix).Qprep_null_ve ./ 100;
    ai(2,sessix) = rez(sessix).Qprep_potent_ve ./ 100;
    ai(3,sessix) = rez(sessix).Qmove_null_ve ./ 100;
    ai(4,sessix) = rez(sessix).Qmove_potent_ve ./ 100;
end


defcols = getColors();
cols(1,:) = defcols.null * 255;
cols(2,:) = defcols.potent * 255;
cols(3,:) = defcols.null * 255;
cols(4,:) = defcols.potent * 255;
cols = cols ./ 255;

figure; ax = gca; hold on;
xs = [1 2 4 5];
for i = 1:size(ai,1)
    temp = ai(i,:);
    h(i) = bar(xs(i),mean(temp));
    h(i).FaceColor = cols(i,:);
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.5;
    scatter(xs(i)*ones(size(temp)),temp,60,'MarkerFaceColor',cols(i,:), ...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
        'MarkerFaceAlpha',0.7)
    errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);
end
% ylim([-0.001 ax.YLim(2)])
ax.XTick = xs;
xlabels  = {'null-delay', 'null-response', 'potent-delay', 'potent-response'};
xticklabels(xlabels);

ylabel('Frac of nVE / Alignment Index')
ax.FontSize = 15;































