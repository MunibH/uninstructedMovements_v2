clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')));
rmpath(genpath(fullfile(utilspth,'fig2/')));
rmpath(genpath(fullfile(utilspth,'figx/')));
rmpath(genpath(fullfile(utilspth,'MotionMapper/')));
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
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off
params.condition(end+1) = {'R&no&~stim.enable&~autowater&~early'};            % no right, no stim, aw off
params.condition(end+1) = {'L&no&~stim.enable&~autowater&~early'};            % no left, no stim, aw off

params.tmin = -2.4;
params.tmax = 2.5;
params.dt = 1/50;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality (except 'garbage')
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
% meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB13_ALMVideo(meta,datapth);
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

%     % mean center across conditions
%     for clu = 1:size(snpsth,2)
%         % find mean at each time point for each condition (time,1)
%         mean_across_cond = mean(snpsth(:,clu,:),3);
%         snpsth(:,clu,:) = snpsth(:,clu,:) - repmat(mean_across_cond,[1,1,size(snpsth,3)]);
%     end
    rez(sessix).psth_processed = snpsth; % soft normed and mean centered
    rez(sessix).psth_processed = rez(sessix).psth;

    % epochs
    % prep
    edges = [-1 -0.2];
    [~,e1] = min(abs(rez(sessix).time - edges(1)));
    [~,e2] = min(abs(rez(sessix).time - edges(2)));
    rez(sessix).prepix = e1:e2;

    % move
    edges = [0.1 1];
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
    
    trialdat = zscore_singleTrialNeuralData(obj(sessix)); % (time,trials,clu)
%     trialdat = obj(sessix).trialdat; % (time,clu,trials)
    rez(sessix).null_trialdat = tensorprod(trialdat,rez(sessix).Qnull,3,1);
    rez(sessix).potent_trialdat = tensorprod(trialdat,rez(sessix).Qpotent,3,1);


    rez(sessix).null_ssm = mean(rez(sessix).null_trialdat.^2,3);
    rez(sessix).potent_ssm = mean(rez(sessix).potent_trialdat.^2,3);

end


%% variance explained plots
close all

for sessix = 1:numel(rez)
    ai(1,sessix) = rez(sessix).Qprep_null_ve ./ 100;
    ai(2,sessix) = rez(sessix).Qprep_potent_ve ./ 100;
    ai(3,sessix) = rez(sessix).Qmove_null_ve ./ 100;
    ai(4,sessix) = rez(sessix).Qmove_potent_ve ./ 100;
end

[objix,uAnm] = groupSessionsByAnimal(meta);
nAnm = numel(uAnm);

for ianm = 1:nAnm
    ix = find(objix{ianm});
    anmVE.null_prep(ianm) = mean(ai(1,ix));
    anmVE.null_move(ianm) = mean(ai(3,ix));
    anmVE.potent_prep(ianm) = mean(ai(2,ix));
    anmVE.potent_move(ianm) = mean(ai(4,ix));
end


defcols = getColors();
clear cols
cols(1,:) = defcols.null * 255;
cols(2,:) = defcols.potent * 255;
cols(3,:) = defcols.null * 255;
cols(4,:) = defcols.potent * 255;
cols = cols ./ 255;

fns = fieldnames(anmVE);

f = figure;
f.Position = [688   501   291   319];
ax = gca; hold on;
xs = [1 2 4 5];
for i = 1:numel(fns)
    temp = anmVE.(fns{i});
    h(i) = bar(xs(i),mean(temp)); % mean across dims and sessions
    cix = i;
    h(i).FaceColor = cols(cix,:);
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.5;
    %     scatter(xs(i)*ones(size(temp)),temp,10,'MarkerFaceColor','k', ...
    %         'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
    %         'MarkerFaceAlpha',1)
    e = errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);
    e.LineWidth = 0.5;
    e.CapSize = 2;

    vs(i) = scatter(randn(nAnm,1) * 0.1 + xs(i)*ones(nAnm,1),temp,10,'MarkerFaceColor','k',...
        'MarkerEdgeColor','k','LineWidth',1);



    xs_(:,i) = vs(i).XData';
    ys_(:,i) = temp;
end

for i = 1:size(xs_,1)
    patchline(xs_(i,1:2),ys_(i,1:2),'EdgeAlpha',0.4,'LineWidth',0.1)
    patchline(xs_(i,3:4),ys_(i,3:4),'EdgeAlpha',0.4,'LineWidth',0.1)
end

% ylim([-0.001 ax.YLim(2)])
ax.XTick = xs;
xlabels  = strrep(fns,'_','-');
xticklabels(xlabels);

ylabel('Fraction of normVE')
ax.FontSize = 12;

title('Covariances')


%% single trial r2 b/w (me,potent) and (me,null)
close all



% correlate magnitude of activity in N/Ps and ME
% plot by mouse and session

ms = {'o','<','^','v','>','square','diamond','o','<'};
sz = 50;

[objix,uAnm]  = groupSessionsByAnimal(meta);
nAnm = numel(uAnm);

cond2use = [2 3];

f = figure;
ax = gca;
hold on;

for ianm = 1:nAnm
    % get sessions for current animal
    ix = find(objix{ianm});

    % for each session, calculate avg feat value over trials
    for isess = 1:numel(ix)
        sessix = ix(isess);

        trix = cell2mat(params(sessix).trialid(cond2use)');

        thisme = me(sessix).data(:,trix);
        null = rez(sessix).null_ssm;
        potent = rez(sessix).potent_ssm;

        r2.null{ianm}(isess) = corr(thisme(:),null(:)).^2;
        r2.potent{ianm}(isess) = corr(thisme(:),potent(:)).^2;
    end

    mu.potent = nanmean(r2.potent{ianm});
    mu.null = nanmean(r2.null{ianm});
    if ianm > 7
        s = scatter(mu.null,mu.potent,sz, ...
            'MarkerEdgeColor','k','MarkerFaceColor','none','Marker',ms{ianm});
    else
        s = scatter(mu.null,mu.potent,sz, ...
            'MarkerEdgeColor','w','MarkerFaceColor','k','Marker',ms{ianm});
    end
end

xlabel('R2(ME,Null)', 'Interpreter','none')
ylabel('R2(ME,Potent)', 'Interpreter','none')
ax.FontSize = 9;
axis(ax,'equal')

mins = min([ax.XLim(1) ax.YLim(1)]);
maxs = max([ax.XLim(2) ax.YLim(2)]);

xlim([0 maxs])
ylim([0 maxs])
ax = gca;
plot(ax.XLim,ax.YLim,'k--','LineWidth',2)
% view([90 -90])


h = zeros(nAnm, 1);
for ianm = 1:numel(h)
    if ianm > 7
        h(ianm) = scatter(NaN,NaN,sz,'MarkerEdgeColor','k','MarkerFaceColor','none','Marker',ms{ianm});
    else
        h(ianm) = scatter(NaN,NaN,sz,'MarkerEdgeColor','none','MarkerFaceColor','k','Marker',ms{ianm});
    end

end
legString = uAnm;

leg = legend(h, legString);
leg.EdgeColor = 'none';
leg.Color = 'none';
leg.Location = 'best';


%% single trial cross-corr b/w (me,potent) and (me,null)
close all
clear r

dt = params(1).dt;

% maxlag = int64( ./ dt);
maxlag = 2 / dt;

tt = [-2 0]; % only use time points before go cue
for i = 1:numel(tt)
    [~,tix(i)] = min(abs(obj(1).time - tt(i)));
end
tix = tix(1):tix(2);

cond2use = [2 3];
for sessix = 1:numel(meta)

    trix = cell2mat(params(sessix).trialid(cond2use)');

%     thisme = zscore(me(sessix).data(:,trix));
    thisme = normalize(me(sessix).data(tix,trix));
    null = normalize(rez(sessix).null_ssm(tix,:));
    potent = normalize(rez(sessix).potent_ssm(tix,:));

    for t = 1:numel(trix)

        [r.null{sessix}(:,t),lagtm] = xcorr(thisme(:,t),null(:,t),maxlag,'normalized');
        r.potent{sessix}(:,t) = xcorr(thisme(:,t),potent(:,t),maxlag,'normalized');
    end


end

mu.null = cellfun(@(x) squeeze(nanmean(x,2)),r.null,'UniformOutput',false); % mean across trials for each session
mu.potent = cellfun(@(x) squeeze(nanmean(x,2)),r.potent,'UniformOutput',false);
nullcc = cell2mat(mu.null);
potentcc = cell2mat(mu.potent);

%%
col = getColors();
lw = 2;
alph = 0.2;
f = figure;
f.Position = [680   694   383   284];
ax = gca;
hold on;
shadedErrorBar(lagtm*dt,mean(nullcc,2),getCI(nullcc),{'Color',col.null,'LineWidth',lw},alph,ax);
shadedErrorBar(lagtm*dt,mean(potentcc,2),getCI(potentcc),{'Color',col.potent,'LineWidth',lw},alph,ax);

xlabel('Time lag, movement and subspace activity (s)')
ylabel('Correlation')
title('ta-elsayed')




%% psth magnitudes in NP

close all
clear col cols
cols = getColors;
col.null{1} = [3, 3, 173] ./ 255;
col.null{2} = [173, 3, 3] ./ 255;

col.potent{1} = [115, 169, 245] ./ 255;
col.potent{2} = [240, 110, 138] ./ 255;

alph = 0.2;
lw = 2;

sm = 11;

xlims = [-2.3,2];

sample = (mode(obj(1).bp.ev.sample) - 2.5);
delay = (mode(obj(1).bp.ev.delay) - 2.5);


cond2use = [2 3];
for sessix = 1:numel(meta)
    f = figure;
    f.Position = [698   436   343   230];
    ax = gca;
    hold on;
    for c = 1:numel(cond2use)
        null = rez(sessix).null_trialdat(:,params(sessix).trialid{cond2use(c)},:);
        null = mean(null.^2,3);
        mu = nanmean(null,2);
        sig = nanstd(null,[],2) ./ sqrt(rez(sessix).dPrep);
        shadedErrorBar(obj(1).time,mu,getCI(null),{'Color',col.null{c},'LineWidth',lw},alph,ax)

        potent = rez(sessix).potent_trialdat(:,params(sessix).trialid{cond2use(c)},:);
        potent = mean(potent.^2,3);
        mu = nanmean(potent,2);
        sig = nanstd(potent,[],2) ./ sqrt(rez(sessix).dMove);
        shadedErrorBar(obj(1).time,mu,getCI(potent),{'Color',col.potent{c},'LineWidth',lw},alph,ax)
    end
    title([meta(sessix).anm ' ' meta(sessix).date])
    ax.FontSize = 12;
    xline(sample,'k--')
    xline(delay,'k--')
    xline(0,'k--')
    xlim(xlims)
    xlabel('Time from go cue (s)')
    ylabel('Activity')

end



%% plot me np single trial heatmaps

close all

edges = [-0.4 -0.02]; % relative to to go cue
for i = 1:numel(edges)
    [~,tix(i)] = min(abs(obj(1).time - edges(i)));
end

sm = 7;

for sessix = 15%:numel(rez)
    temprez = rez(sessix);
    trix = cell2mat(params(sessix).trialid(2:3)');

    tempme = me(sessix).data(:,trix);
    
    null = rez(sessix).null_ssm;
    potent = rez(sessix).potent_ssm;
   
    [~,ix] = sort(mean(tempme(tix(1):tix(2),:),1),'descend');
    trix = ix;

    f = figure;
    f.Position = [680   591   321   387];

%     ax1 = subplot(1,3,1);
    plotme = mySmooth(tempme(:,trix),2);
    imagesc(obj(sessix).time,1:numel(trix),plotme');
    colormap(linspecer)
    lims = clim;
%     clim([lims(1) lims(2) / 1])
    colorbar;
%     ylabel('Trials')

%     ax2 = subplot(1,3,2);
    f = figure;
    f.Position = [680   591   321   387];
    potent = mySmooth(normalize(potent(:,trix)),sm);
    imagesc(obj(sessix).time,1:numel(trix),potent');
    colormap(linspecer);
    lims = clim;
%     clim([lims(1) lims(2) / 1.5])
    colorbar;
%     xlabel('Time from go cue (s)')

%     ax3 = subplot(1,3,3);
    f = figure;
    f.Position = [680   591   321   387];
    null = mySmooth(normalize(null(:,trix)),sm);
    imagesc(obj(sessix).time,1:numel(trix), null');
    colormap(linspecer)
    lims = clim;
%     clim([lims(1) lims(2) / 1.5])
    colorbar;

%     sgtitle([meta(sessix).anm ' ' meta(sessix).date])
%     pause
%     close all

end










