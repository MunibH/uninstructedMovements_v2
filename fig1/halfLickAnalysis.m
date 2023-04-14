clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig3/')))
rmpath(genpath(fullfile(utilspth,'mc_stim/')))
rmpath(genpath(fullfile(utilspth,'MotionMapper/')))

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
params.condition(1)     = {'hit&~stim.enable&~autowater&~early'};            % not early
params.condition(end+1) = {'hit&~stim.enable&~autowater&early'};             % early
params.condition(end+1) = {'hit&~stim.enable&~autowater'};             % early

params.tmin = -2.4;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance


params.advance_movement = 0.0;

params.behav_only = 0;


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
meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params, params.behav_only);

% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%%





%% Grand average PSTH
close all

cond2use = [1 2];

c = linspecer(numel(cond2use),'qualitative');
ls = '-';
lw = 1;
alph = 0.2;

f = figure;
for i = 1:numel(meta)
    ax = nexttile;
    hold on;
    for j = 1:numel(cond2use)
        dat = obj(i).psth(:,:,cond2use(j));
        mu = mean(dat,2);
        sd = std(dat,[],2) ./ sqrt(size(dat,2));
        shadedErrorBar(obj(i).time,mu,sd,{'Color',c(j,:),'LineWidth',lw,'LineStyle',ls},alph,ax)
    end
    xlabel('Time from go cue (s)')
    ylabel('spks / s')
    title([meta(i).anm ' ' meta(i).date])
end
sgtitle('Grand Average PSTHs')


%% tongue trajectories for early and not early trials
close all

feat = 'tongue_length';
featix = find(ismember(kin(1).featLeg,feat));


cond2use = [1 2];

c = linspecer(numel(cond2use),'qualitative');
ls = '-';
lw = 0.1;
alph = 0.2;

xlims = [-2.3 2];

f = figure;
for i = 1:numel(meta)
    ax = nexttile;
    hold on;
    for j = 1:numel(cond2use)
        trix = params(i).trialid{cond2use(j)};
        dat = kin(i).dat(:,trix,featix);
        for t = 1:numel(trix)
            patchline(obj(i).time,dat(:,t),'EdgeAlpha',alph,'EdgeColor',c(j,:),'LineWidth',lw,'LineStyle',ls)
        end
    end
    xlim(xlims)
    xlabel('Time from go cue (s)')
    ylabel('Tongue length')
    title([meta(i).anm ' ' meta(i).date])
end


%% find half licks
% for every trial find
% time of half licks (not touching lickport)
% time of first 4 post go cue licks
% get neural data aligned to these times

feat = 'tongue_length';
% feat = 'top_tongue_xdisp_view2';
% feat = 'top_tongue_ydisp_view2';
featix = find(ismember(kin(1).featLeg,feat));

lt_thresh = 0.03; % lick duration must be at least 30 ms to count (half or full lick)
samp_thresh = floor(lt_thresh / params(1).dt);

clear rez 

for sessix = 1:numel(meta)
    nTrials = obj(sessix).bp.Ntrials;
    rez(sessix).fullLickStartIX    = cell(nTrials,1);
    rez(sessix).fullLickStartTimes = cell(nTrials,1);
    rez(sessix).halfLickStartIX    = cell(nTrials,1);
    rez(sessix).halfLickStartTimes = cell(nTrials,1);
    rez(sessix).fullLickEndIX      = cell(nTrials,1);
    rez(sessix).fullLickEndTimes   = cell(nTrials,1);
    rez(sessix).halfLickEndIX      = cell(nTrials,1);
    rez(sessix).halfLickEndTimes   = cell(nTrials,1);
    for trix = 1:nTrials
 
        gc = obj(sessix).bp.ev.goCue(trix);

        % get bpod lick times
        lickL = obj(sessix).bp.ev.lickL{trix};
        lickR = obj(sessix).bp.ev.lickR{trix};
        lt = sort([lickL lickR]);

        % align lick times
        lt = lt - gc;


        
  
        % get tongue traj
        tong = kin(sessix).dat(:,trix,featix);

        [istart, iend, len] = ZeroOnesCount(fillmissing(tong,'constant',0)); % start ix, end ix, and nSamples (len), for each tongue protrusion
        % -- if not using tong len
%         baseTong = median(tong);
%         temp = tong;
%         temp(tong==baseTong) = 0;
%         temp(tong~=baseTong) = 1;
%         [istart, iend, len] = ZeroOnesCount(temp); % start ix, end ix, and nSamples (len), for each tongue protrusion

        
        % delete tong protrusions that are less than a threshold long in
        % duration
        toDelete = find(len < samp_thresh);
        istart(toDelete) = [];
        iend(toDelete) = [];
        len(toDelete) = [];
        istarttm = obj(sessix).time(istart);
        iendtm = obj(sessix).time(iend);

        % skip trial if no tongue protrusions found
        if isempty(istarttm) || isempty(iendtm)
            continue
        end

        % for each tongue protrusion, determine if it's a half or full lick
        isFull = logical(zeros(size(istarttm)));
        for iTouch = 1:numel(lt)
            [~,ix] = min(abs(istarttm - lt(iTouch))); % find index in istarttm of the tongue protrusion closest in time to current touch time
            % check if the closest protrusion time contains the touch time in its range
            chk = (istarttm(ix) <= lt(iTouch)) && (  (iendtm(ix) >= lt(iTouch)) ||  (abs(iendtm(ix) - lt(iTouch)) <= 0.005) ); % buffer of 5 ms after protrusion end time since there may be imprecisions in kin
            if chk
                isFull(ix) = true;
            end
        end
        
        % now for current trial, we have when each protrusion starts and
        % ends, and if it is a full or incomplete lick
        rez(sessix).fullLickStartIX{trix}    = istart(isFull);
        rez(sessix).fullLickStartTimes{trix} = istarttm(isFull);
        rez(sessix).halfLickStartIX{trix}    = istart(~isFull);
        rez(sessix).halfLickStartTimes{trix} = istarttm(~isFull);
        rez(sessix).fullLickEndIX{trix}    = iend(isFull);
        rez(sessix).fullLickEndTimes{trix} = iendtm(isFull);
        rez(sessix).halfLickEndIX{trix}    = iend(~isFull);
        rez(sessix).halfLickEndTimes{trix} = iendtm(~isFull);

% %     find lick times pre and post go cue
% %     restrict post go cue to first 4 licks
%         nLicks = min(4,numel(postlt));
%         postlt = postlt(1:nLicks);
%         prelt = lt(lt<0);
%         postlt = lt(lt>0);

    clear isFull

    end
end

fullLickAll = [];
halfLickAll = [];
for i = 1:numel(meta)
    fullLickSession = cell2mat(rez(i).fullLickStartTimes')';
    halfLickSession = cell2mat(rez(i).halfLickStartTimes')';
    fullLickAll = [fullLickAll ; fullLickSession];
    halfLickAll = [halfLickAll ; halfLickSession];
end

%% PLOT (when does each lick type occur)

close all



% plot distribution of times when each lick type occurs
c = linspecer(2,'qualitative');
alph = 1;
nbins = 100;

f = figure;
t = tiledlayout('flow');
ax = nexttile; hold on;
histogram(fullLickAll,nbins,'facecolor',c(1,:),'edgecolor','none','facealpha',alph,'Normalization','probability')
title('full licks')
xline(mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.goCue) - mode(obj(1).bp.ev.goCue),'k--')

ax = nexttile; hold on;
histogram(halfLickAll,nbins,'facecolor',c(2,:),'edgecolor','none','facealpha',alph,'Normalization','probability')
title('incomplete licks')
xline(mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.goCue) - mode(obj(1).bp.ev.goCue),'k--')

ax = nexttile; hold on;
histogram(fullLickAll(fullLickAll<0),nbins/2,'facecolor',c(1,:),'edgecolor','none','facealpha',alph,'Normalization','probability')
title('full licks')
xline(mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.goCue) - mode(obj(1).bp.ev.goCue),'k--')

ax = nexttile; hold on;
histogram(halfLickAll(halfLickAll<0),nbins/2,'facecolor',c(2,:),'edgecolor','none','facealpha',alph,'Normalization','probability')
title('incomplete licks')
xline(mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue),'k--')
xline(mode(obj(1).bp.ev.goCue) - mode(obj(1).bp.ev.goCue),'k--')

ylabel(t,'prob')
xlabel(t,'Time from go cue (s)')


%% neural activity aligned to full or incomplete licks
close all

clear fulldat halfdat

win.pre = 40;
win.post = 40;
win.end = numel(obj(1).time);
win.N = win.pre + win.post + 1;

% first by session
for sessix = 1:numel(meta)
    dat = rez(sessix);
    trialdat = obj(sessix).trialdat;

    nFull = numel(cell2mat(dat.fullLickStartIX'));
    nHalf = numel(cell2mat(dat.halfLickStartIX'));

    rez(sessix).fulldat = nan(win.N, size(trialdat,2), nFull);
    rez(sessix).halfdat = nan(win.N, size(trialdat,2), nHalf);

    fct = 1;
    hct = 1;
    for trix = 1:size(trialdat,3)
        full = dat.fullLickStartIX{trix};
        half = dat.halfLickStartIX{trix};

        % only look at pre go cue stuff
        [~,gc] = min(abs(obj(sessix).time - 0));
        full = full(full<gc);
        half = half(half<gc);
        
        for i = 1:numel(full)
            if (full(i)+win.post)>win.end
                ix = full(i)-win.pre:win.end;
                tmp = trialdat(ix,:,trix);
                offset = win.N - size(tmp,1);
                tmp = cat(1,tmp,nan(offset,size(tmp,2)));
            elseif (full(i)-win.post)<1
                ix = 1:full(i)+win.post;
                tmp = trialdat(ix,:,trix);
                offset = win.N - size(tmp,1);
                tmp = cat(1,nan(offset,size(tmp,2)),tmp);
            else
                ix = full(i)-win.pre:full(i)+win.post;
                tmp = trialdat(ix,:,trix);
            end
 
            rez(sessix).fulldat(:,:,fct) = tmp;
%             fulldat(1:numel(ix),:,fct) = trialdat(ix,:,trix);
            fct = fct + 1;
        end

        for i = 1:numel(half)
            if (half(i)+win.post)>win.end
                ix = half(i)-win.pre:win.end;
                tmp = trialdat(ix,:,trix);
                offset = win.N - size(tmp,1);
                tmp = cat(1,tmp,nan(offset,size(tmp,2)));
            elseif (half(i)-win.post)<1
                ix = 1:half(i)+win.post;
                tmp = trialdat(ix,:,trix);
                offset = win.N - size(tmp,1);
                tmp = cat(1,nan(offset,size(tmp,2)),tmp);
            else
                ix = half(i)-win.pre:half(i)+win.post;
                tmp = trialdat(ix,:,trix);
            end
 
            rez(sessix).halfdat(:,:,fct) = tmp;
%             halfdat(1:numel(ix),:,fct) = trialdat(ix,:,trix);
            fct = fct + 1;
        end

    end
end


%
% f = figure;
fullmu = [];
halfmu = [];
for sessix = 1:numel(meta)

%     ax = nexttile;
%     hold on;

    temp = rez(sessix).fulldat;
    mu = nanmean(temp,3);
    % plot(mu,'k','LineWidth',0.2)
%     plot(nanmean(mu,2),'k','LineWidth',0.2)
    fullmu = cat(2,fullmu,mu);

    temp = rez(sessix).halfdat;
    mu = nanmean(temp,3);
    % plot(mu,'r','LineWidth',0.2)
%     plot(nanmean(mu,2),'r','LineWidth',0.2)
    halfmu = cat(2,halfmu,mu);

%     xlabel('Time from lick (ms)')
%     ylabel('spks / s')

%     title([meta(sessix).anm ' ' meta(sessix).date])

end

alph = 0.2;
c = linspecer(2,'qualitative');

f = figure;
ax = gca;
hold on;

time = (-win.pre:win.post) .* params(1).dt;

temp = fullmu;
mu = nanmean(temp,2);
sd = nanstd(temp,[],2) ./ sqrt(size(temp,2));
shadedErrorBar(time,mu,sd,{'Color',c(1,:),'LineWidth',1,'LineStyle','-'},alph,ax)


temp = halfmu;
mu = nanmean(temp,2);
sd = nanstd(temp,[],2) ./ sqrt(size(temp,2));
shadedErrorBar(time,mu,sd,{'Color',c(2,:),'LineWidth',1,'LineStyle','-'},alph,ax)

xline(0,'k--')

xlabel('Time from first tongue visible (s)')
ylabel('spks / s')

title('early lick aligned activity before go cue ')


















