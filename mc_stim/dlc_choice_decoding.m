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


%% TODO
% - created new objs from pipeline and nidq data - only use those
% - need to handle frameTimes appropriately now (subtract 0.5 secs)
% - for choice decoding, just use avg during stim period

%% default params (same for each session)

dfparams = [];

% -- time alignment params --
% dfparams.alignEv = 'goCue';
% dfparams.times = [-2.5 2.5]; % relative to goCue
dfparams.alignEv = 'goCue';
dfparams.times = [-2.5 2.5]; % relative to delay

dfparams.dt_vid = 0.0025;
dfparams.time = dfparams.times(1):dfparams.dt_vid:dfparams.times(2);

dfparams.warp = 0; % 0 means no warping, 1 means warp delay period to 


% -- trial type params --
% control
dfparams.cond(1) =     {'R&hit&~stim.enable&~autowater&~autolearn&~early'};  % right hits, no stim, aw off
dfparams.cond(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % left miss, no stim, aw off
dfparams.cond(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
dfparams.cond(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % right miss, no stim, aw off
% stim
dfparams.cond(end+1) = {'R&hit&stim.enable&~autowater&~autolearn&~early'};  % right hits, stim, aw off
dfparams.cond(end+1) = {'L&miss&stim.enable&~autowater&~early'};            % left miss, stim, aw off
dfparams.cond(end+1) = {'L&hit&stim.enable&~autowater&~early'};             % left hits, stim, aw off
dfparams.cond(end+1) = {'R&miss&stim.enable&~autowater&~early'};            % right miss, stim, aw off

% -- stim types --
dfparams.stim.types = {'Bi_MC','Right_MC','Left_MC','Bi_ALM','Bi_M1TJ','Right_ALM','Right_M1TJ','Left_ALM','Left_M1TJ'}; % ALM_Bi is MC_Bi
% dfparams.stim.num   = logical([0 0 0 1 0 0 0 0 0]); % Bi_ALM
dfparams.stim.num   = logical([0 0 0 0 1 0 0 0 0]); % Bi_M1TJ
% dfparams.stim.num   = logical([0 0 0 1 1 0 0 0 0]); % Bi_ALM & Bi_M1TJ

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


obj = loadObjs(meta);


%% find trials for each condition

for i = 1:numel(meta)
    params(i).trialid = findTrials(obj(i), dfparams.cond);
end

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

%% choice decoding from dlc features (ctrl trials)

clearvars -except datapth kin kinfeats kinfeats_norm me meta obj params dfparams acc_ctrl acc_stim

% params
rez.nFolds = 4; % number of iterations (bootstrap)

rez.binSize = 30; % ms
rez.dt = floor(rez.binSize / (dfparams.dt_vid*1000)); % samples
rez.tm = dfparams.time(1:rez.dt:numel(dfparams.time));
rez.numT = numel(rez.tm);

rez.train = 1; % fraction of trials to use for training (1-train for testing)

% match number of right and left hits, and right and left misses
cond2use = 1:4;
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

        acc_ctrl(:,sessix,ifeat) = DLC_ChoiceDecoder(in,rez,trials);


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






