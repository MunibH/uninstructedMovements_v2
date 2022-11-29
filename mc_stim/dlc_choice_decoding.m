clear, clc, close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig3/')))
rmpath(genpath(fullfile(utilspth,'fig1/')))
rmpath(genpath(fullfile(utilspth,'musall2019/')))

% add paths for figure specific functions
addpath(genpath(pwd))
rmpath(genpath(fullfile(utilspth,'mc_stim/mc_stim_funcs/oldGetKin/')))


clc


%% default params (same for each session)

dfparams = [];

% -- time alignment params --
dfparams.alignEv = 'goCue';
dfparams.times = [-3 2.5]; % relative to goCue
% dfparams.alignEv = 'delay';
% dfparams.times = [-1.0 2.5]; % relative to delay

dfparams.dt_vid = 0.0025;
dfparams.time = dfparams.times(1):dfparams.dt_vid:dfparams.times(2);

dfparams.warp = 0; % 0 means no warping, 1 means warp delay period to

% dfparams.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
%     {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

dfparams.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

dfparams.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

dfparams.advance_movement = 0.0;


% -- trial type params --
dfparams.cond(1)     = {'R&hit&~stim.enable&~autowater&~autolearn&~early'};       % R hits, no stim, no autowater, no autolearn
dfparams.cond(end+1) = {'L&hit&~stim.enable&~autowater&~autolearn&~early'};       % L hits, no stim, no autowater, no autolearn

dfparams.cond(end+1) = {'R&miss&~stim.enable&~autowater&~autolearn&~early'};      % R misses, no stim, no autowater, no autolearn
dfparams.cond(end+1) = {'L&miss&~stim.enable&~autowater&~autolearn&~early'};      % L misses, no stim, no autowater, no autolearn

dfparams.cond(end+1) = {'R&hit&stim.enable&~autowater&~autolearn&~early'};        % R hits, stim, no autowater, no autolearn
dfparams.cond(end+1) = {'L&hit&stim.enable&~autowater&~autolearn&~early'};        % L hits, stim, no autowater, no autolearn

dfparams.cond(end+1) = {'R&miss&stim.enable&~autowater&~autolearn&~early'};       % R misses, stim, no autowater, no autolearn
dfparams.cond(end+1) = {'L&miss&stim.enable&~autowater&~autolearn&~early'};       % L misses, stim, no autowater, no autolearn

% -- stim types --
dfparams.stim.types = {'Bi_MC','Right_MC','Left_MC','Bi_ALM','Bi_M1TJ','Right_ALM','Right_M1TJ','Left_ALM','Left_M1TJ'};
% dfparams.stim.num   = logical([1 1 1 1 1 1 1 1 1]);   % ALL
% dfparams.stim.num   = logical([0 0 0 0 0 1 1 1 1]);   % Right_ALM / Left_ALM / Right_M1TJ / Left_M1TJ
% dfparams.stim.num   = logical([0 0 0 1 0 0 0 0 0]);   % Bi_ALM
dfparams.stim.num   = logical([0 0 0 0 1 0 0 0 0]);   % Bi_M1TJ
% dfparams.stim.num   = logical([1 0 0 0 0 0 0 0 0]);   % Bi_MC
% dfparams.stim.num   = logical([0 0 0 1 1 0 0 0 0]);   % Bi_M1TJ Bi_ALM
% dfparams.stim.num   = logical([0 1 1 0 0 0 0 0 0]);   % Right_MC
% dfparams.stim.num   = logical([0 0 1 0 0 0 0 0 0]);   % Left_MC
% dfparams.stim.num   = logical([0 0 0 0 0 1 0 0 0]);   % Right_ALM
% dfparams.stim.num   = logical([0 0 0 0 0 0 0 1 0]);   % Left_ALM
% dfparams.stim.num   = logical([0 0 0 0 0 0 1 0 0]);   % Right_M1TJ
% dfparams.stim.num   = logical([0 0 0 0 0 0 0 0 1]);   % Left_M1TJ

dfparams.stim.pow.types = [1.08, 3.14, 8]; % mW, 3.14 for all unilateral sessions and most bilateral sessions
dfparams.stim.pow.num = logical([0 1 0]);


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

% meta = meta(1:5);

% LOAD DATA OBJECTS
obj = loadObjs(meta);

% TRIALS | MOTION ENERGY | KINEMATICS
for sessix = 1:numel(meta)
    % trials
    params(sessix).trialid = findTrials(obj(sessix), dfparams.cond);

    % motion energy
    if isstruct(obj(sessix).me)
        me(sessix).data = obj(sessix).me.data;
        me(sessix).moveThresh = obj(sessix).me.moveThresh;
    else
        me(sessix).data = obj(sessix).me;
        me(sessix).moveThresh = nan;
    end
    me(sessix) = processME(obj(sessix),me(sessix),dfparams);

    % kinematics
    obj(sessix).time = dfparams.time;
    params(sessix).traj_features = dfparams.traj_features;
    params(sessix).advance_movement = dfparams.advance_movement;
    params(sessix).feat_varToExplain = dfparams.feat_varToExplain;
    params(sessix).alignEvent = dfparams.alignEv;
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end

%% choice decoding from dlc features (train on ctrl, test on stim)

clearvars -except datapth kin kinfeats kinfeats_norm me meta obj params dfparams

% params
rez.nFolds = 10; % number of iterations (bootstrap)

rez.binSize = 50; % ms
rez.dt = floor(rez.binSize / (dfparams.dt_vid*1000)); % samples
rez.tm = dfparams.time(1:rez.dt:numel(dfparams.time));
rez.numT = numel(rez.tm);

rez.train = 1; % fraction of trials to use for training ( (1-train) becomes test fraction)

rez.nShuffles = 20;

% match number of right and left hits, and right and left misses
rez.cond.ctrl.hit = 1:2;
rez.cond.ctrl.miss = 3:4;
rez.cond.stim.hit = 5:6;
rez.cond.stim.miss = 7:8;



% featGroups = {{'tongue'},...
%     {'jaw'},...
%     {'paw'},...
%     {'nose','nostril'},...
%     {'motion_energy'}};


featGroups = {'all'};
% featGroups = {'tongue'};
% featGroups = {'jaw'}; % [x]
% featGroups = {'nose','nostril'}; % [x]
% featGroups = {'paw'};
% featGroups = {'motion_energy'};

% featGroups = {'tongue'};

for ifeat = 1:numel(featGroups)
    disp(['Feature Group ' num2str(ifeat) '/' num2str(numel(featGroups))])

    rez.feats2use = featGroups{ifeat};

    if strcmpi(rez.feats2use,'all')
        rez.featix = (1:size(kin(1).featLeg))';
    else
        if numel(rez.feats2use) == 1
            rez.feats2use = rez.feats2use';
        end
        if ischar(rez.feats2use)
            rez.feats2use = {rez.feats2use};
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

        % trials, % labels (1 for right choice, -1 for left choice)

        trialid = params(sessix).trialid;
        % match number of hits, misses (separately for ctrl and stim)

        % ctrl hit trials
        nRight = numel(trialid{1});
        nLeft = numel(trialid{2});
        nTrials = min(nRight,nLeft);
        trialid{1} = randsample(trialid{1},nTrials,false);
        trialid{2} = randsample(trialid{2},nTrials,false);

        % ctrl miss trials
        nRight = numel(trialid{3});
        nLeft = numel(trialid{4});
        nTrials = min(nRight,nLeft);
        trialid{3} = randsample(trialid{3},nTrials,false);
        trialid{4} = randsample(trialid{4},nTrials,false);

        % stim hit trials
        nRight = numel(trialid{5});
        nLeft = numel(trialid{6});
        nTrials = min(nRight,nLeft);
        trialid{5} = randsample(trialid{5},nTrials,false);
        trialid{6} = randsample(trialid{6},nTrials,false);

        % stim miss trials
        nRight = numel(trialid{7});
        nLeft = numel(trialid{8});
        nTrials = min(nRight,nLeft);
        trialid{7} = randsample(trialid{7},nTrials,false);
        trialid{8} = randsample(trialid{8},nTrials,false);



        trials.ctrl.all = [trialid{1} ; ... % right hits
            trialid{2} ; ... % left hits
            trialid{3} ; ... % right miss
            trialid{4} ; ... % left miss
            ];

        Y.ctrl = [ones(numel(trialid{1}),1) ; ...
            zeros(numel(trialid{2}),1) ; ...
            zeros(numel(trialid{3}),1) ; ...
            ones(numel(trialid{4}),1) ; ...
            ];

        trials.stim.all = [trialid{5} ; ... % right hits
            trialid{6} ; ... % left hits
            trialid{7} ; ... % right miss
            trialid{8} ; ... % left miss
            ];

        % labels (1 for right choice, 0 for left choice)
        Y.stim = [ones(numel(trialid{5}),1) ; ...
            zeros(numel(trialid{6}),1) ; ...
            zeros(numel(trialid{7}),1) ; ...
            ones(numel(trialid{8}),1) ; ...
            ];

        % input
        % use all features
        X.ctrl = kin(sessix).dat_std(:,trials.ctrl.all,rez.featix);
        X.stim = kin(sessix).dat_std(:,trials.stim.all,rez.featix);
        % fill missing values
        for featix = 1:size(X.ctrl,3)
            X.ctrl(:,:,featix) = fillmissing(X.ctrl(:,:,featix),"constant",0);
            X.stim(:,:,featix) = fillmissing(X.stim(:,:,featix),"constant",0);
        end

        % train/test split
        trials.train = trials.ctrl.all;
        trials.test  = trials.stim.all;

        in.train.y = Y.ctrl;
        in.test.y = Y.stim;
        in.train.X = X.ctrl;
        in.test.X = X.stim;


        %         [trials.train,trials.trainidx] = datasample(trials.ctrl.all,round(numel(trials.ctrl.all)*rez.train),'Replace',false);
        %         trials.testidx = find(~ismember(trials.all,trials.train));
        %         trials.test = trials.all(trials.testidx);

        %         in.train.y = Y(trials.trainidx);
        %         in.test.y  = Y(trials.testidx);
        %         in.train.X = X(:,trials.trainidx,:);
        %         in.test.X  = X(:,trials.testidx,:);

        % decoding

        [acc_ctrl(:,sessix,ifeat), acc_stim(:,sessix,ifeat)] = ...
            DLC_ChoiceDecoder_TrainOnNonStim_TestOnStim(in,rez);



        %     % shuffle labels for a 'null' distribution
        %
        %
        Y = randsample(Y,numel(Y));

        % train/test split

        in.train.y = randsample(Y.ctrl,numel(Y.ctrl));
        in.test.y  = randsample(Y.stim,numel(Y.stim));


        for ishuf = 1:rez.nShuffles
            [acc_shuf(:,sessix,ifeat),~] = ...
                DLC_ChoiceDecoder_TrainOnNonStim_TestOnStim(in,rez);
        end

    end

end

acc_shuf_ = reshape(acc_shuf,size(acc_shuf,1),size(acc_shuf,2)*size(acc_shuf,3));


%% plot

cols = {'k','c',[0.6,0.6,0.6]};

alph = 0.15;

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);
trialStart = mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue);

figure;
ax = gca;
hold on;

shadedErrorBar(rez.tm(1:end-1),mean(acc_ctrl,2),std(acc_ctrl,[],2)./sqrt(numel(obj)),{'Color',cols{1},'LineWidth',2},alph,ax)
shadedErrorBar(rez.tm(1:end-1),mean(acc_stim,2),std(acc_stim,[],2)./sqrt(numel(obj)),{'Color',cols{2},'LineWidth',2},alph,ax)
shadedErrorBar(rez.tm(1:end-1),mean(acc_shuf_,2),std(acc_shuf_,[],2)./sqrt(numel(obj)),{'Color',cols{3},'LineWidth',2},alph,ax)
xline(0,'k:','LineWidth',2)
xline(sample,'k:','LineWidth',2)
xline(delay,'k:','LineWidth',2)
xline(trialStart,'k:','LineWidth',2)
xlim([trialStart, dfparams.times(2)-0.2])
ylim([ax.YLim(1) 1])

xlabel('Time (s) from go cue')
ylabel([num2str(rez.nFolds) '-Fold CV Accuracy'])
title('Choice Decoding from DLC Features')

h = zeros(3, 1);
for i = 1:numel(h)
    h(i) = plot(NaN,NaN,'-','Color',cols{i},'LineWidth',2);
end
legString = {'ctrl','stim','shuf'};

leg = legend(h, legString);
leg.EdgeColor = 'none';
leg.Location = 'best';

ax.FontSize = 15;


%% plot

cols = linspecer(size(acc_ctrl,3));

alph = 0.5;

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);

figure;
ax = gca;
hold on;
for ifeat = 1:size(acc_ctrl,3)
    temp = acc_ctrl(:,:,ifeat);
    shadedErrorBar(rez.tm(1:end-1),mean(temp,2),std(temp,[],2)./sqrt(numel(obj)),{'Color',cols(ifeat,:),'LineWidth',2},alph,ax)
end
xline(0,'k:','LineWidth',2)
xline(sample,'k:','LineWidth',2)
xline(delay,'k:','LineWidth',2)

ylim([ax.YLim(1) 1])

xlabel('Time (s) from go cue')
ylabel([num2str(rez.nFolds) '-Fold CV Accuracy'])
title('Choice Decoding from DLC Features')

h = zeros(numel(featGroups), 1);
for i = 1:numel(h)
    h(i) = plot(NaN,NaN,'-','Color',cols(i,:),'LineWidth',2);
end
legString = cellfun(@(x) strrep(x{1},'_',' '), featGroups,'UniformOutput',false);

leg = legend(h, legString);
leg.EdgeColor = 'none';
leg.Location = 'best';

% same plot but for each session

cols = linspecer(size(acc_ctrl,3));

alph = 0.5;

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);

figure;
t = tiledlayout('flow');
for sessix = 1:size(acc_ctrl,2)
    ax = nexttile;
    hold on
    temp = mySmooth(squeeze(acc_ctrl(:,sessix,:)),21);
    for ifeat = 1:size(acc_ctrl,3)
        temp2 = temp(:,ifeat);
        shadedErrorBar(rez.tm(1:end-1),mean(temp2,2),std(temp2,[],2)./sqrt(numel(obj)),{'Color',cols(ifeat,:),'LineWidth',2},alph,ax)
    end
    title([meta(sessix).anm ' ' meta(sessix).date])
    xline(0,'k:','LineWidth',2)
    xline(sample,'k:','LineWidth',2)
    xline(delay,'k:','LineWidth',2)

    ylim([0.4 1])

    %     if sessix == 1
    %         h = zeros(numel(featGroups), 1);
    %         for i = 1:numel(h)
    %             h(i) = plot(NaN,NaN,'-','Color',cols(i,:),'LineWidth',2);
    %         end
    %         legString = cellfun(@(x) strrep(x{1},'_',' '), featGroups,'UniformOutput',false);
    %
    %         leg = legend(h, legString);
    %         leg.EdgeColor = 'none';
    %         leg.Location = 'best';
    %     end
end


xlabel(t,'Time (s) from go cue')
ylabel(t,[num2str(rez.nFolds) '-Fold CV Accuracy'])


%% choice decoding from dlc features (stim trials)

clearvars -except datapth kin kinfeats kinfeats_norm me meta obj params dfparams acc_ctrl acc_stim

% params
rez.nFolds = 4; % number of iterations (bootstrap)

rez.binSize = 30; % ms
rez.dt = floor(rez.binSize / (dfparams.dt_vid*1000)); % samples
rez.tm = dfparams.time(1:rez.dt:numel(dfparams.time));
rez.numT = numel(rez.tm);

rez.train = 1; % fraction of trials to use for training (1-train for testing)

% match number of right and left hits, and right and left misses
cond2use = 5:8;
hitcond = [1 3];
misscond = [2 4];

featGroups = {{'tongue'},...
    {'jaw','trident'},...
    {'nose','nostril'},...
    {'paw'},...
    {'motion_energy'}};

for ifeat = 1:numel(featGroups)
    disp(['Feature Group ' num2str(ifeat) '/' num2str(numel(featGroups))])
    %     rez.feats2use = kin(1).featLeg;
    % rez.feats2use = {'jaw_ydisp_view1'};
    % rez.feats2use = {'motion_energy'};
    % rez.feats2use = {'view2'};
    rez.feats2use = featGroups{ifeat};
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



    for sessix = 1:numel(obj)
        disp(['Decoding session ' num2str(sessix) ' / ' num2str(numel(obj))])

        % trials

        trials_cond = params(sessix).trialid(cond2use);

        minHitTrials = cellfun(@(x) numel(x),trials_cond(hitcond), 'UniformOutput',false);
        nhits = min(cell2mat(minHitTrials));

        minMissTrials = cellfun(@(x) numel(x),trials_cond(misscond), 'UniformOutput',false);
        nmiss = min(cell2mat(minMissTrials));


        trials_hit = cellfun(@(x) randsample(x,nhits), trials_cond(hitcond), 'UniformOutput', false);
        trialsHit = cell2mat(trials_hit);
        trialsHit = trialsHit(:);

        trials_miss = cellfun(@(x) randsample(x,nmiss), trials_cond(misscond), 'UniformOutput', false);
        trialsMiss = cell2mat(trials_miss);
        trialsMiss = trialsMiss(:);

        trials.all = [trialsHit ; trialsMiss];

        % labels (1 for right choice, 0 for left choice)
        Y = [ones(nhits,1) ; -ones(nhits,1) ; ones(nmiss,1) ; -ones(nmiss,1)]; % right hits, left hits, left miss, right miss

        % input
        % use all features
        X = kinfeats{sessix}(:,trials.all,rez.featix);
        % fill missing values
        for featix = 1:size(X,3)
            X(:,:,featix) = fillmissing(X(:,:,featix),"constant",0);
        end

        % train/test split

        [trials.train,trials.trainidx] = datasample(trials.all,round(numel(trials.all)*rez.train),'Replace',false);
        trials.testidx = find(~ismember(trials.all,trials.train));
        trials.test = trials.all(trials.testidx);

        in.train.y = Y(trials.trainidx);
        in.test.y  = Y(trials.testidx);
        in.train.X = X(:,trials.trainidx,:);
        in.test.X  = X(:,trials.testidx,:);

        % decoding

        acc_stim(:,sessix,ifeat) = DLC_ChoiceDecoder(in,rez,trials);


        % %     % shuffle labels for a 'null' distribution
        % %
        % %
        % %     Y = randsample(Y,numel(Y));
        % %
        % %     % train/test split
        % %
        % %     in.train.y = Y(trials.trainidx);
        % %     in.test.y  = Y(trials.testidx);
        % %
        % %     acc_shuffled(:,sessix) = DLC_ChoiceDecoder(in,rez,trials);


    end

end

% plot

cols = linspecer(size(acc_stim,3));

alph = 0.5;

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);

figure;
ax = gca;
hold on;
for ifeat = 1:size(acc_stim,3)
    temp = acc_stim(:,:,ifeat);
    shadedErrorBar(rez.tm(1:end-1),nanmean(temp,2),nanstd(temp,[],2)./sqrt(numel(obj)),{'Color',cols(ifeat,:),'LineWidth',2},alph,ax)
end
xline(0,'k:','LineWidth',2)
xline(sample,'k:','LineWidth',2)
xline(delay,'k:','LineWidth',2)

ylim([ax.YLim(1) 1])

xlabel('Time (s) from go cue')
ylabel([num2str(rez.nFolds) '-Fold CV Accuracy'])
title('Choice Decoding from DLC Features')

h = zeros(numel(featGroups), 1);
for i = 1:numel(h)
    h(i) = plot(NaN,NaN,'-','Color',cols(i,:),'LineWidth',2);
end
legString = cellfun(@(x) strrep(x{1},'_',' '), featGroups,'UniformOutput',false);

leg = legend(h, legString);
leg.EdgeColor = 'none';
leg.Location = 'best';

% same plot but for each session

cols = linspecer(size(acc_stim,3));

alph = 0.5;

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);

figure;
t = tiledlayout('flow');
for sessix = 1:size(acc_stim,2)
    ax = nexttile;
    hold on
    temp = mySmooth(squeeze(acc_stim(:,sessix,:)),21);
    for ifeat = 1:size(acc_stim,3)
        temp2 = temp(:,ifeat);
        shadedErrorBar(rez.tm(1:end-1),nanmean(temp2,2),nanstd(temp2,[],2)./sqrt(numel(obj)),{'Color',cols(ifeat,:),'LineWidth',2},alph,ax)
    end
    title([meta(sessix).anm ' ' meta(sessix).date])
    xline(0,'k:','LineWidth',2)
    xline(sample,'k:','LineWidth',2)
    xline(delay,'k:','LineWidth',2)

    ylim([0.4 1])

    %     if sessix == 1
    %         h = zeros(numel(featGroups), 1);
    %         for i = 1:numel(h)
    %             h(i) = plot(NaN,NaN,'-','Color',cols(i,:),'LineWidth',2);
    %         end
    %         legString = cellfun(@(x) strrep(x{1},'_',' '), featGroups,'UniformOutput',false);
    %
    %         leg = legend(h, legString);
    %         leg.EdgeColor = 'none';
    %         leg.Location = 'best';
    %     end
end


xlabel(t,'Time (s) from go cue')
ylabel(t,[num2str(rez.nFolds) '-Fold CV Accuracy'])

%% plot ctrl - stim

% close all

% plot

delta_acc = acc_ctrl - acc_stim;

cols = linspecer(size(delta_acc,3));

alph = 0.1;

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);

figure;
ax = gca;
hold on;
for ifeat = 1:size(delta_acc,3)
    temp = delta_acc(:,:,ifeat);
    shadedErrorBar(rez.tm(1:end-1),nanmean(temp,2),nanstd(temp,[],2)./sqrt(numel(obj)),{'Color',cols(ifeat,:),'LineWidth',2},alph,ax)
end
xline(0,'k:','LineWidth',2)
xline(sample,'k:','LineWidth',2)
xline(delay,'k:','LineWidth',2)

% ylim([ax.YLim(1) 1])

xlabel('Time (s) from go cue')
ylabel([num2str(rez.nFolds) '-Fold CV Accuracy (Ctrl - Stim)'])
title('Choice Decoding from DLC Features')

h = zeros(numel(featGroups), 1);
for i = 1:numel(h)
    h(i) = plot(NaN,NaN,'-','Color',cols(i,:),'LineWidth',2);
end
legString = cellfun(@(x) strrep(x{1},'_',' '), featGroups,'UniformOutput',false);

leg = legend(h, legString);
leg.EdgeColor = 'none';
leg.Location = 'best';

% same plot but for each session

cols = linspecer(size(delta_acc,3));

alph = 0.1;

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);

figure;
t = tiledlayout('flow');
for sessix = 1:size(delta_acc,2)
    ax = nexttile;
    hold on
    temp = mySmooth(squeeze(delta_acc(:,sessix,:)),21);
    for ifeat = 1:size(delta_acc,3)
        temp2 = temp(:,ifeat);
        shadedErrorBar(rez.tm(1:end-1),nanmean(temp2,2),nanstd(temp2,[],2)./sqrt(numel(obj)),{'Color',cols(ifeat,:),'LineWidth',2},alph,ax)
    end
    title([meta(sessix).anm ' ' meta(sessix).date])
    xline(0,'k:','LineWidth',2)
    xline(sample,'k:','LineWidth',2)
    xline(delay,'k:','LineWidth',2)

    ylim([-0.1 0.5])

    %     if sessix == 1
    %         h = zeros(numel(featGroups), 1);
    %         for i = 1:numel(h)
    %             h(i) = plot(NaN,NaN,'-','Color',cols(i,:),'LineWidth',2);
    %         end
    %         legString = cellfun(@(x) strrep(x{1},'_',' '), featGroups,'UniformOutput',false);
    %
    %         leg = legend(h, legString);
    %         leg.EdgeColor = 'none';
    %         leg.Location = 'best';
    %     end
end


xlabel(t,'Time (s) from go cue')
ylabel(t,[num2str(rez.nFolds) '-Fold CV Accuracy (Ctrl - Stim)'])






