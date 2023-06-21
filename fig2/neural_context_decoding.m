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
params.condition(1)     = {'R&hit&~stim.enable&~autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&hit&~stim.enable&autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % left hits, no stim, aw off

params.tmin = -3;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


% params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
%     {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance


params.advance_movement = 0.0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);

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

%% context decoding from dlc features

clearvars -except datapth kin me meta obj params acc acc_shuf

% params
rez.nFolds = 4; % number of iterations (bootstrap)

rez.binSize = 75; % ms
rez.dt = floor(rez.binSize / (params(1).dt*1000)); % samples
rez.tm = obj(1).time(1:rez.dt:numel(obj(1).time));
rez.numT = numel(rez.tm);

rez.train = 1; % fraction of trials to use for training (1-train for testing)

rez.nShuffles = 4;

% match number of right and left hits, and right and left misses
cond2use = 1:4;
afccond = [1 2];
awcond = [3 4];

% featGroups = {{'tongue'},...
%     {'jaw','trident'},...
%     {'nose','nostril'},...
%     {'paw'},...
%     {'motion_energy'}};

featGroups = {'all'};

for ifeat = 1:numel(featGroups)
    disp(['Feature Group ' num2str(ifeat) '/' num2str(numel(featGroups))])
    %     rez.feats2use = kin(1).featLeg;
    % rez.feats2use = {'jaw_ydisp_view1'};
    % rez.feats2use = {'motion_energy'};
    % rez.feats2use = {'view2'};


    rez.feats2use = featGroups{ifeat};

    if strcmpi(rez.feats2use,'all')
        rez.featix = 1:numel(kin(1).featLeg);
    else

        if size(rez.feats2use,1) == 1
            rez.feats2use = rez.feats2use';
        end

        [~,mask] = patternMatchCellArray(kin(1).featLeg,rez.feats2use,'any');
        if sum(mask) == 0
            [~,mask] = patternMatchCellArray(kin(1).featLeg,rez.feats2use,'all');
        end
        if sum(mask) == 0
            use = zeros(size(kin(1).featLeg));
            for i = 1:numel(kin(1).featLeg)
                for j = 1:numel(rez.feats2use)
                    if contains(kin(1).featLeg{i},rez.feats2use{j})
                        use(i) = 1;
                    end
                end
            end
            mask = logical(use);
        end
        if sum(mask) == 0
            error('didnt find features to use')
        end
        rez.featix = find(mask);

    end



    for sessix = 1:numel(obj)
        disp(['Decoding session ' num2str(sessix) ' / ' num2str(numel(obj))])

        % trials

        trials_cond = params(sessix).trialid(cond2use);

        minAFCTrials = cellfun(@(x) numel(x),trials_cond(afccond), 'UniformOutput',false);
        nafc = min(cell2mat(minAFCTrials));

        minAWTrials = cellfun(@(x) numel(x),trials_cond(awcond), 'UniformOutput',false);
        naw = min(cell2mat(minAWTrials));

        nTrials = min(nafc,naw);

        trials_afc = cellfun(@(x) randsample(x,nTrials), trials_cond(afccond), 'UniformOutput', false);
        trialsAFC = cell2mat(trials_afc);
        trialsAFC = trialsAFC(:);

        trials_aw = cellfun(@(x) randsample(x,nTrials), trials_cond(awcond), 'UniformOutput', false);
        trialsAW = cell2mat(trials_aw);
        trialsAW = trialsAW(:);

        trials.all = [trialsAFC ; trialsAW];

        % labels (1 for right choice, 0 for left choice)
        Y = [ones(nTrials,1) ; ones(nTrials,1) ; -ones(nTrials,1) ; -ones(nTrials,1)]; % right hits, left hits, left miss, right miss
        % Y = [ones(nhits,1) ; -ones(nhits,1) ; ones(nmiss,1) ; -ones(nmiss,1)]; % right hits, left hits, left miss, right miss

        % input
        X = obj(sessix).trialdat(:,:,trials.all);
        X = permute(X,[1 3 2]); % (time,trials, clu)
        Xmax = max(X(:));
        Xmin = min(X(:));
        Range = Xmax - Xmin;
        Xnrm = ((X - Xmin)/Range - 0.5) * 2;
        X = Xnrm; clear Xnrm;

        % train/test split
        [trials.train,trials.trainidx] = datasample(trials.all,round(numel(trials.all)*rez.train),'Replace',false);
        trials.testidx = find(~ismember(trials.all,trials.train));
        trials.test = trials.all(trials.testidx);

        in.train.y = Y(trials.trainidx);
        in.test.y  = Y(trials.testidx);
        in.train.X = X(:,trials.trainidx,:);
        in.test.X  = X(:,trials.testidx,:);

        % decoding

        acc(:,sessix,ifeat) = DLC_ChoiceDecoder(in,rez,trials);


        % shuffle labels for a 'null' distribution


        Y = randsample(Y,numel(Y));

        % train/test split

        in.train.y = Y(trials.trainidx);
        in.test.y  = Y(trials.testidx);

        for ishuf = 1:rez.nShuffles
            acc_shuf(:,sessix,ishuf,ifeat) = DLC_ChoiceDecoder(in,rez,trials);
        end


    end

end

acc_shuf_ = reshape(acc_shuf,size(acc_shuf,1),size(acc_shuf,2)*size(acc_shuf,3));


%% t-test

for itime = 1:size(acc,1)
    [h(itime),p(itime)] = ttest2(acc(itime,:),acc_shuf_(itime,:));
end
figure;
hold on
plot(rez.tm(1:end-1),p,'LineWidth',2)


%% plot

close all

cols = {'k',[0.6,0.6,0.6]};

alph = 0.2;

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);
trialStart = mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue);

f = figure;
f.Position = [644   517   335   258];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

tempacc = cat(1,mean(acc(1,:),2)*ones(size(acc)),acc);
ctrl = mySmooth(tempacc, 11,'reflect');
ctrl = ctrl(size(tempacc)/2+1:end,:);
acc_shuff_ = cat(2,acc_shuf_,acc_shuf_);
shuffed = mySmooth(acc_shuff_, 15,'reflect');

shadedErrorBar(rez.tm(1:end-1),mean(ctrl,2),getCI(ctrl),{'Color',cols{1},'LineWidth',2},alph,ax)
shadedErrorBar(rez.tm(1:end-1),mean(shuffed,2),getCI(shuffed,0)/2,{'Color',cols{2},'LineWidth',2},alph,ax)

% shadedErrorBar(rez.tm(1:end-1),mean(acc,2),std(acc,[],2)./sqrt(numel(obj)),{'Color',cols{1},'LineWidth',2},alph,ax)
% shadedErrorBar(rez.tm(1:end-1),mean(acc_shuf_,2),std(acc_shuf_,[],2)./sqrt(numel(obj)),{'Color',cols{2},'LineWidth',2},alph,ax)
xx=xline(0,'k--'); uistack(xx,'bottom');
xx=xline(sample,'k--'); uistack(xx,'bottom');
xx=xline(delay,'k--'); uistack(xx,'bottom');
xx=xline(trialStart,'k--'); uistack(xx,'bottom');
xlim([trialStart+0.1, params(1).tmax-0.2])
ylim([ax.YLim(1) 1])

xlabel('Time from go cue (s)')
ylabel([num2str(rez.nFolds) '-Fold CV Accuracy'])
title('Neural context decoding','FontSize',8)

% h = zeros(2, 1);
% for i = 1:numel(h)
%     h(i) = plot(NaN,NaN,'-','Color',cols{i},'LineWidth',2);
% end
% legString = {'Kinematics','Shuffled labels'};
%
% leg = legend(h, legString);
% leg.EdgeColor = 'none';
% leg.Location = 'best';
% leg.Color = 'none';
% ax.FontSize = 10;


