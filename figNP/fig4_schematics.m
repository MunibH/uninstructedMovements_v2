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

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;


%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\munib\Documents\Economo-Lab\data';

meta = [];

% --- ALM ---
% meta = loadJEB6_ALMVideo(meta,datapth);
% meta = loadJEB7_ALMVideo(meta,datapth); % selectivity in ME
% % meta = loadEKH1_ALMVideo(meta,datapth); % selectivity in ME
% meta = loadEKH3_ALMVideo(meta,datapth); % selectivity in ME
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB13_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth); % selectivity in ME % go cue is at 2.3 instead of 2.5 like all other sessions??
meta = loadJEB15_ALMVideo(meta,datapth);


% --- M1TJ ---
% meta = loadJEB13_M1TJVideo(meta,datapth);

meta = meta(2);

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


%%
% close all

% f = figure;
% f.Position = [437   139   498   659];
% ax = gca;
% hold on;

ts = mode(obj.bp.ev.bitStart) - 2.5;
sample = mode(obj.bp.ev.sample) - 2.5;
delay = mode(obj.bp.ev.delay) - 2.5;
gc = 0;

nsm = 5;
msm = 7;
nclu = 5;

cols = getColors;
lw = 2;
alph = 0.23;
alph2 = 0.13;
cscale = 0.8;

% while true
for trix = 58
    if obj.bp.early(trix)
        continue
    end
    clu = randsample(size(obj.psth,2),nclu);
    %     trix = randsample(obj.bp.Ntrials,1);
%     trix = 139; [58 164 170]
    clu(1) = 36;
    clu(2) = 23;
    clu(3) = 22;
    clu(4) = 51;
    clu(5) = 47;
    
    me = stitchAndPurgeMoveBouts(me,params,0.05,0.01);

    m = mySmooth(normalize(me.data(:,trix),'range',[0 2]),msm);
    move = me.move(:,trix);
    mnull = m;
    mpotent = m;
    mnull(~move,:) = nan;
    mpotent(move,:) = nan;


    for i = 1:nclu
        n{i} = mySmooth(normalize(squeeze(obj.trialdat(:,clu(i),trix)),'range',[0 1]),nsm);
        nnull{i} = n{i};
        npotent{i} = n{i};
        nnull{i}(~move,:) = nan;
        npotent{i}(move,:) = nan;
    end

    f = figure;
    f.Position = [437   139   498   659];
    ax = gca;
    hold on;

    cc = cols.potent * cscale;
    plot(obj.time,m,'Color',cc,'LineWidth',lw)
    cc = cols.null * cscale;
    plot(obj.time,mpotent,'Color',cc,'LineWidth',lw)
%     for i = 1:nclu
%         plot(obj.time,n{i}-(1.5*i),'Color',cols.potent,'LineWidth',lw)
%         plot(obj.time,npotent{i}-(1.5*i),'Color',cols.null,'LineWidth',lw)
%     end
%     plot(obj.time,m,'Color','k','LineWidth',lw)
%     plot(obj.time,mpotent,'Color','k','LineWidth',lw)
    for i = 1:nclu
        plot(obj.time,n{i}-(1.5*i),'Color','k','LineWidth',lw)
        plot(obj.time,npotent{i}-(1.5*i),'Color','k','LineWidth',lw)
    end
    title(['Trial ' num2str(trix) ' - Clu ' num2str(clu')])
    xline(ts,'k--'); xline(delay,'k--'); xline(sample,'k--'); xline(gc,'k--')
    xlim([-2.35 2])

    y = ax.YLim;
    Y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];

    [nstart, nend, nlen] = ZeroOnesCount(~move);
    nend = nend + 1;
    if nend(end) > numel(obj.time)
        nend(end) = nend(end) - 1;
    end
    ns = obj.time(nstart);
    ne = obj.time(nend);
    for ii = 1:numel(nstart)
        ff(ii) = fill([ns(ii) ne(ii) ne(ii) ns(ii)], Y, cols.null);
        ff(ii).FaceAlpha = alph;
        ff(ii).EdgeColor = 'none';
        uistack(ff(ii),'bottom');
        ff2(ii) = fill([ns(ii) ne(ii) ne(ii) ns(ii)], [-0.3 -0.3 ax.YLim(2) ax.YLim(2)], cols.null);
        ff2(ii).FaceAlpha = alph2;
        ff2(ii).EdgeColor = 'none'; 
    end
    [mstart, mend, mlen] = ZeroOnesCount(move);
    mend = mend + 1;
    if mend(end) > numel(obj.time)
        mend(end) = mend(end) - 1;
    end
    ms = obj.time(mstart);
    me_ = obj.time(mend);
    for ii = 1:numel(mstart)
        ff(ii) = fill([ms(ii) me_(ii) me_(ii) ms(ii)], Y, cols.potent);
        ff(ii).FaceAlpha = alph;
        ff(ii).EdgeColor = 'none';
        uistack(ff(ii),'bottom');
        ff2(ii) = fill([ms(ii) me_(ii) me_(ii) ms(ii)], [-0.3 -0.3 ax.YLim(2) ax.YLim(2)], cols.potent);
        ff2(ii).FaceAlpha = alph2;
        ff2(ii).EdgeColor = 'none'; 
%         uistack(ff2(ii),'bottom');
    end
%     uistack(ff,'bottom')
    ax.Visible = 'off';
%     break
%     pause
%     cla(ax);
end
