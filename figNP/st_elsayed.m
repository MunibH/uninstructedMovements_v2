clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = '/Users/munib/Economo-Lab/code/uninstructedMovements_v3';
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
params.dt = 1/100;

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

datapth = '/Users/munib/Economo-Lab/data';

meta = [];

% --- ALM ---
% meta = loadJEB6_ALMVideo(meta,datapth);
% meta = loadJEB7_ALMVideo(meta,datapth); % selectivity in ME
% % % meta = loadEKH1_ALMVideo(meta,datapth); % selectivity in ME
% meta = loadEKH3_ALMVideo(meta,datapth); % selectivity in ME
meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB13_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth); % selectivity in ME % go cue is at 2.3 instead of 2.5 like all other sessions??
% meta = loadJEB15_ALMVideo(meta,datapth);

meta = meta(1);

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
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj, nullalltime, onlyAW, delayOnly);

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

%% plot coding dimensions from activity N/P spaces
close all

sav = 0;

plotmiss = 0;
plotno = 0;

titlestring = 'Null';
axnull = plotCDProj(cd_null_all,cd_null,sav,titlestring,plotmiss,plotno);
% plotCDVarExp(cd_null_all,sav,titlestring)
% plotSelectivity(cd_null_all,cd_null,sav,titlestring)
% plotSelectivityExplained(cd_null_all,cd_null,sav,titlestring)

titlestring = 'Potent';
axpotent = plotCDProj(cd_potent_all,cd_potent,sav,titlestring,plotmiss,plotno);
% plotCDVarExp(cd_potent_all,sav,titlestring)
% plotSelectivity(cd_potent_all,cd_potent,sav,titlestring)
% plotSelectivityExplained(cd_potent_all,cd_potent,sav,titlestring)

% match y axis limits across null and potent space CDs
for i = 1:numel(axnull)
    ys = [ min(axnull(i).YLim(1),axpotent(i).YLim(1)) max(axnull(i).YLim(2),axpotent(i).YLim(2)) ];
    axnull(i).YLim = ys;
    axpotent(i).YLim = ys;
end

titlestring = 'Null | Potent CDs';
% plotCDProj_NP(cd_potent_all,cd_null_all,cd_potent,cd_null,sav,titlestring,plotmiss)
% plotSelectivityExplained_NP(cd_potent_all,cd_null_all,cd_potent,cd_null,sav,titlestring)


%% selectivity in various epochs in cd late hits and errors
close all

plotCDSelectivity(meta,obj,cd_null_all,cd_null,cd_potent_all,cd_potent);


% plotCD_EpochSelectivity(meta,obj,cd_null_all,cd_null,cd_potent_all,cd_potent);


%% variance explained
close all

% plotVarExplained_Recon(rez); % R2 (same as VE) using reconstructed data
% plotVarianceExplained_NP(rez,meta);
plotVarianceExplained_DelayResponse(rez,meta);

%% single trial magnitude activity in N/P
close all

cond2use = 2:3;
plotME_NPMagnitude_singleTrials(meta,obj,params,me,rez,cond2use)

% corr_ME_NPMagnitude(meta,obj,params,me,rez);

%% trial-averaged magnitude of activity in N/P
close all

cond2use = [2 3];
plotSubpsaceMagnitude_TrialAvg(obj,meta,params,rez,cond2use)

% plotSubpsaceMagnitude_TrialAvg_byRT(obj,meta,params,rez,cond2use)

%% session-averaged magnitude of activity in N/P
close all

cond2use = [2 3];
plotSubpsaceMagnitude_SessionAvg(obj,meta,params,rez,me,cond2use)

% plotSubpsaceMagnitude_TrialAvg_byRT(obj,meta,params,rez,cond2use)

%% magnitude of activity in N/P aligned to movement bouts


stitch_dist = 0.05; % in seconds, stitch together movement bouts shorter than this % 0.025
purge_dist = 0.1; % in seconds, remove move bouts shorter than this value, after stitching complete % 0.1
tbout = 0.3; % move/non-move bout/transition required to be at least this long in seconds % 0.3
[dat.mdat,dat.mdat_leg,dat.qdat,dat.qdat_leg,newme] = nonGC_moveTransitions(obj,me,params,stitch_dist,purge_dist,tbout);

% plot

close all

ndims = 10;
sav = 0;
% plot_nonGC_moveTransitions_trialAvg_v6(dat,obj,newme,rez,params,meta,ndims) % sumsqmag across dimensions, move to quiet
plot_nonGC_moveTransitions_trialAvg_v7(dat,obj,newme,rez,params,meta,ndims) % sumsqmag across dimensions, quiet to move


% plotSubspaceActivityAlignedToMoveBoutEnd(dat,obj,me,rez,params,meta)



%% plot NP projections (trial-averaged and single-trials) TODO

close all

sav = 0;
% -----------------------------------------------------------------------
% -- Null and Potent Space Single Trial Projections --
% -----------------------------------------------------------------------

% % - projections showing move / quiet somehow
cond2use = [2 3]; % right hits, left hits
ndims = 4; % top ndims variance explaining dimensions
% plotSingleTrialNPHeatmaps(rez,params,me,ndims,cond2use,meta);

% -----------------------------------------------------------------------
% -- Null and Potent Space Trial-averaged Projections --
% -----------------------------------------------------------------------

ndims = 5; % how many n/p dimensions to plot, in order of variance explained
cond2plot = [1 2]; % right hit, left hit
% plot_NP_PSTH(rez,obj,params,ndims,cond2plot,meta)

%% PCA on trial-avg N/P activity

close all

cond2use = 1:2; % right/left hits, corresponding to rez.N_potent/null_psth
pcaNP(meta,obj,params,rez,cond2use)

% pause, press any key to cycle through sessions


%% selectivity in N/P activity
close all

cond2use_ta = [2 3]; % right and left hits, corresponding to trial-avg projs onto n/p
cond2use_st = [2 3]; % right and left hits, corresponding to single-trial projs onto n/p
subTrials = 35;
plotSelectivityNeurons(meta,obj,params,cond2use_ta,cond2use_st) % to get number of selective cells

% cond2use_ta = [1 2]; % right and left hits, corresponding to trial-avg projs onto n/p
% cond2use_st = [8 9]; % right and left hits, corresponding to single-trial projs onto n/p
% % plotSelectivityNP(meta,obj,params,rez,cond2use_ta,cond2use_st)
%
% % hits
% subTrials = 35;
% plotSelectivityNPPref(meta,obj,params,rez,cond2use_ta,cond2use_st, subTrials)

% misses
% subTrials = 10; % doesn't show anything interesting and also no error trials
% cond2use_ta = [3 4]; % right,left miss, corresponding to trial-avg projs onto n/p
% plotSelectivityNPPref(meta,obj,params,rez,cond2use_ta,cond2use_st, subTrials)

%% plot reconstructed PSTHS

plotReconstructedPSTHs(obj,rez,params)


%% decode choice from N/P

cv.nFolds = 4; % number of cv folds
cv.binSize = 50; % ms
cv.dt = floor(cv.binSize / (params(1).dt*1000)); % samples
cv.tm = obj(1).time(1:cv.dt:numel(obj(1).time));
cv.numT = numel(cv.tm);
cv.train = 1; % fraction of trials to use for training (1-train for testing)
cv.nShuffles = 2;

cond2use = [8 9]; % right and left hits
choiceDecoder_v2(obj,rez,params,cv,cond2use)


%% decode choice from CD N/P


cv.nFolds = 4; % number of cv folds
cv.binSize = 50; % ms
cv.dt = floor(cv.binSize / (params(1).dt*1000)); % samples
cv.tm = obj(1).time(1:cv.dt:numel(obj(1).time));
cv.numT = numel(cv.tm);
cv.train = 1; % fraction of trials to use for training (1-train for testing)
cv.nShuffles = 2;

cond2use = [8 9]; % right and left hits
cdix = 1; % late (not orthogonal to early)
choiceDecoder_CD_NP(obj,rez,cd_null,cd_potent,params,cv,cond2use,cdix)
title('Hits')

% how much better do we get from adding error trials
% CD_NP_late (not orth to early) would predict that you don't gain any
% decodability in the null space during delay period, but you would in the
% potent space
% cond2use = [8 9 10 11]; % right and left hits
% cdix = 1; % late (not orthogonal to early)
% choiceDecoder_CD_NP_withErrors(obj,rez,cd_null,cd_potent,params,cv,cond2use,cdix)
% title('Hits and Errors')



%% how much variance does NP_ramping_mode explain of the reconstructed data from NP

plotVarExpRampingNP(cd_null,cd_potent,obj,rez,params)


%% ME, Mag_NP, Mag_fullActivity
% looking at cortical disengagement stuff here

cond2use = [8 9 10 11]; % right, left hits, no early

plotMagnitudeCluPotent_ME(obj,rez,params,me,cond2use)



%% Reaction time stuff

rt = firstJawRT(obj);
% rt = firstTongueRT(obj);

% plotRT_ME_singleTrials(obj,params,me,rt)

nQuantiles = 3;
% cdSelectivityByRTQuartile(obj,cd_null,cd_potent,params,nQuantiles,rt);
[h,p] = npSelectivityByRTQuartile(meta,obj,rez,params,nQuantiles,rt,me);

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

    %     thisme = zscore(me(sessix).data(tix,trix));
    thisme = normalize(me(sessix).data(tix,trix));
    %     null = sum(rez(sessix).N_null(:,trix,:).^2,3);
    %     potent = sum(rez(sessix).N_potent(:,trix,:).^2,3);

    for t = 1:numel(trix)

        % null
        ndims = rez(sessix).dPrep;
        for d = 1:ndims
            null = normalize(rez(sessix).N_null(tix,trix(t),d).^2);
            [r.null{sessix}(:,t,d),lagtm] = xcorr(thisme(:,t),null,maxlag,'normalized'); % {session}(time,trial,dim)
        end

        % potent
        ndims = rez(sessix).dMove;
        for d = 1:ndims
            potent = normalize(rez(sessix).N_potent(tix,trix(t),d).^2);
            [r.potent{sessix}(:,t,d),lagtm] = xcorr(thisme(:,t),potent,maxlag,'normalized');
        end
    end


end

mu.null = cellfun(@(x) squeeze(nanmean(nanmean(x,3),2)),r.null,'UniformOutput',false); % mean across trials for each session
mu.potent = cellfun(@(x) squeeze(nanmean(nanmean(x,3),2)),r.potent,'UniformOutput',false);
nullcc = cell2mat(mu.null);
potentcc = cell2mat(mu.potent);

% plots
col = getColors();
lw = 2;
alph = 0.2;
f = figure;
f.Position = [680   694   383   284];
ax = gca;
hold on;
% shadedErrorBar(lagtm*dt,mean(nullcc,2),std(nullcc,[],2)./sqrt(numel(meta)),{'Color',col.null,'LineWidth',lw},alph,ax);
% shadedErrorBar(lagtm*dt,mean(potentcc,2),std(potentcc,[],2)./sqrt(numel(meta)),{'Color',col.potent,'LineWidth',lw},alph,ax);
shadedErrorBar(lagtm*dt,mean(nullcc,2),getCI(nullcc),{'Color',col.null,'LineWidth',lw},alph,ax);
shadedErrorBar(lagtm*dt,mean(potentcc,2),getCI(potentcc),{'Color',col.potent,'LineWidth',lw},alph,ax);
% xlim([-0.5 0.5])

xlabel('Time lag, movement and subspace activity (s)')
ylabel('Correlation')
title('st-elsayed')

%% onset of selectivity in potent space and motion energy (just exploring)

close all
clear null potent me_

cond2use_ta = [1 2]; % right and left hits, corresponding to trial-avg projs onto n/p
cond2use_st = [8 9]; % right and left hits, corresponding to single-trial projs onto n/p
% plotSelectivityNP(meta,obj,params,rez,cond2use_ta,cond2use_st)

% hits
subTrials = 35;

for i = 1:numel(obj)
    [null{i},potent{i},me_.right{i}, me_.left{i}] = plotSelectivityNPPrefOnset(meta(i),obj(i),params(i),me(i),rez(i),cond2use_ta,cond2use_st, subTrials);
end

sm = 21;
smtype = 'reflect';
null = mySmooth(cell2mat(null),sm,smtype);
potent = mySmooth(cell2mat(potent),sm,smtype);
me_.right = mySmooth(cell2mat(me_.right),sm,smtype);
me_.left = mySmooth(cell2mat(me_.left),sm,smtype);



xlims = [-2.4 0];

f = figure;

ax = nexttile;
imagesc(obj(1).time,1:numel(obj),null'); %colormap(linspecer);
colorbar;
xlim(xlims);

ax = nexttile;
imagesc(obj(1).time,1:numel(obj),potent'); %colormap(linspecer);
colorbar;
clim([0 4])
xlim(xlims);

ax = nexttile;
imagesc(obj(1).time,1:numel(obj),(me_.right - me_.left)'); %colormap(linspecer);
colorbar;
xlim(xlims);

ax = nexttile;
imagesc(obj(1).time,1:numel(obj),me_.right'); %colormap(linspecer);
xlim(xlims);

ax = nexttile;
imagesc(obj(1).time,1:numel(obj),me_.left'); %colormap(linspecer);
xlim(xlims);

%
% temp = me_.right - me_.left;
% for i = 1:numel(obj)
%     cc.null(i) = corr(temp(:,i),null(:,i));
%     cc.potent(i) = corr(temp(:,i),potent(:,i));
% end
%
% col = getColors;
% f = figure;
% ax = gca;
% hold on;
% nb = 10;
% histogram(cc.null,nb,'EdgeColor','none','FaceColor',col.null,'FaceAlpha',0.5);
% histogram(cc.potent,nb,'EdgeColor','none','FaceColor',col.potent,'FaceAlpha',0.5);


%% plot selectivity in AW trials (right - left selectivity) (for supp4s)
close all

sel = [];

sm = 5;
smtype = 'zeropad';

edges = [-0.5 0];
for i = 1:numel(edges)
    [~,ix(i)] = min(abs(obj(1).time - edges(i)));
end

cond2use = [13 14];
for sessix = 1:numel(meta)
    trix = params(sessix).trialid(cond2use);
    for c = 1:numel(trix)
        psth{c} = mySmooth(mean(obj(sessix).trialdat(:,:,trix{c}),3),sm,smtype);
    end

    % find pref
    temp = cell2mat(cellfun(@(x) nanmean(x(ix(1):ix(2),:),2), psth, 'UniformOutput',false));


    sel = cat(2,sel,psth{1}-psth{2});

    %     break

end

lw = 2;
alph = 0.2;

f = figure;
ax = gca;
mu = nanmean(sel,2);
% sd = nanstd(sel,[],2) ./ sqrt(size(sel,2)) * ;
sd = getCI(sel)*1.96;
shadedErrorBar(obj(1).time,mu,sd,{'Color','k','LineWidth',lw},alph,ax)
xlim([-2.3 2])

tstart = mode(obj(1).bp.ev.bitStart) - 2.5;
sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;
gc = 0;

xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--'); xline(gc,'k--');
xlabel('Time from go cue (s)')
ylabel('Selectivity (spikes / s)')


%% magnitude squared n/p

sav = 1;

mag = getMagnitudeAcrossDimsNP(rez,21);

cols = getColors;

tstart = mode(obj(1).bp.ev.bitStart) - 2.5;
sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;
gc = 0;


close all
for sessix = 1:numel(mag)

    f = figure;
    f.Position = [301         532        1138         308];
    t = tiledlayout('flow');

    ax = nexttile;
    imagesc(obj(sessix).time,1:obj(sessix).bp.Ntrials,mag(sessix).null')
    colormap(linspecer)
    title('Null')
    xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--'); xline(gc,'k--');

    ax = nexttile;
    imagesc(obj(sessix).time,1:obj(sessix).bp.Ntrials,mag(sessix).potent')
    colormap(linspecer)
    title('Potent')

    xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--'); xline(gc,'k--');



    ax = nexttile;
    hold on;
    nr = mean(mag(sessix).null(:,params(sessix).trialid{2}),2);
    pr = mean(mag(sessix).potent(:,params(sessix).trialid{2}),2);
    nl = mean(mag(sessix).null(:,params(sessix).trialid{3}),2);
    pl = mean(mag(sessix).potent(:,params(sessix).trialid{3}),2);
    plot(obj(sessix).time,nr,'Color',cols.rmiss,'LineWidth',2)
    plot(obj(sessix).time,nl,'Color',cols.lmiss,'LineWidth',2)
    plot(obj(sessix).time,pr,'Color',cols.rhit_aw,'LineWidth',2)
    plot(obj(sessix).time,pl,'Color',cols.lhit_aw,'LineWidth',2)
    xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--'); xline(gc,'k--');


    xlabel(t,'Time from go cue (s)')
    ylabel(t,'Sum sq. activity (a.u.)')
    sgtitle([meta(sessix).anm ' ' meta(sessix).date])


    if sav
        pth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3\figNP\figs\fig4\magnitudeNP';
        fn = [meta(sessix).anm '_' meta(sessix).date];
        mysavefig(f,pth,fn)
    end

end
















