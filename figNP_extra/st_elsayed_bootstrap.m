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

% for projections
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off (6)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off (7)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off (8)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off (9)

params.condition(end+1) = {'R&~stim.enable&~autowater&~early'}; % (10)
params.condition(end+1) = {'L&~stim.enable&~autowater&~early'}; % (11)

% for ramping
params.condition(end+1) = {'hit&~stim.enable&~autowater'};               % all hits, no stim, aw off (12)

% autowater
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};             % right hits, no stim, aw off (13)
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};             % left hits, no stim, aw off  (14)
params.condition(end+1) = {'R&miss&~stim.enable&autowater'};             % right hits, no stim, aw off (15)
params.condition(end+1) = {'L&miss&~stim.enable&autowater'};             % left hits, no stim, aw off  (16)

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

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth); % not usable b/c no usable left miss trials
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
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
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
%     kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end



%% Null and Potent Space

clearvars -except obj meta params me sav datapth kin rt

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -----------------------------------------------------------------------
disp('finding null and potent spaces')
for sessix = 1:numel(meta)
    disp([meta(sessix).anm ' ' meta(sessix).date])
    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix));

    % -- null and potent spaces
    cond2use = [2:5 13:16]; % right hit, left hit, right miss, left miss, 2afc and aw
    cond2proj = [2:5 13:16];
    nullalltime = 0; % use all time points to estimate null space if 1
    onlyAW = 0; % only use AW trials
    delayOnly = 0; % only use delay period
    responseOnly = 0;
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj, nullalltime, onlyAW, delayOnly, responseOnly);
end
disp('DONE')


%% BOOTSTRAP PARAMS

clear boot bootobj bootparams

boot.iters = 1000; % number of bootstrap iterations (most papers do 1000)

boot.N.anm = 5; % number of animals to sample w/ replacement
boot.N.sess = 2; % number of sessions to sample w/ replacement (if Nsessions for an animal is less than this number, sample Nsessions)
boot.N.trials.hit = 50; % number of hit trials to sample 
boot.N.trials.miss = 20; % number of miss trials to sample 
boot.N.clu = 20; % number of neurons to sample 



%%
clear samp cd_null cd_potent cd_null_shuf cd_potent_shuf cd_null_all cd_potent_all cd_null_all_shuf cd_potent_all_shuf

for iboot = 1:boot.iters
    disp(['Iteration ' num2str(iboot) '/' num2str(boot.iters)]);

    % randomly sample animals
    [objixs,uAnm] = groupSessionsByAnimal(meta);
    samp.anm = randsample(uAnm,boot.N.anm,true);

    % randomly sample sessions
    for ianm = 1:boot.N.anm
        objix = find(objixs{ismember(uAnm,samp.anm{ianm})});
        samp.sessix{ianm} = randsample(objix,boot.N.sess,true);
        %         if numel(objix) == 1 % only one session for the current animal, use that session
        %             samp.sessix{ianm} = objix;
        %         else % more than one session, sample with replacement
        %             samp.sessix{ianm} = randsample(objix,boot.N.sess,true);
        %         end
    end

    % randomly sample trials
    for ianm = 1:boot.N.anm
        for isess = 1:numel(samp.sessix{ianm})
            sessix = samp.sessix{ianm}(isess);
            trialid = params(sessix).trialid;
            for icond = 1:numel(trialid)
                cond = params(sessix).condition{icond};
                if contains(cond,{'R&hit','L&hit'})
                    nTrials2Sample = boot.N.trials.hit;
                elseif contains(cond,{'R&miss','L&miss'})
                    nTrials2Sample = boot.N.trials.miss;
                else
                    nTrials2Sample = boot.N.trials.hit;
                end
                if numel(trialid{icond}) >= nTrials2Sample % if there are enough trials, sample without replacement
                    samp.trialid{ianm}{isess}{icond} = randsample(trialid{icond},nTrials2Sample,true);
                else % if there are not enough trials, sample with replacement
                    try % catches no response trial types with 0 trials
                        samp.trialid{ianm}{isess}{icond} = randsample(trialid{icond},nTrials2Sample,true);
                    catch
                        samp.trialid{ianm}{isess}{icond} = randsample(trialid{icond},numel(trialid{icond}),true);
                    end
                end
            end
        end
    end

    % randomly sample neurons (from null and potent reconstructions)
    for ianm = 1:boot.N.anm
        for isess = 1:numel(samp.sessix{ianm})
            sessix = samp.sessix{ianm}(isess);
            cluid = 1:size(obj(sessix).psth,2); % index of cluster in psth,trialdat,etc -> different from bootstrapping in fig1 b/c we don't need to re-do getSeq, just going to grab data from rez(sessix)
            samp.cluid{ianm}{isess} = randsample(cluid,boot.N.clu,true);
        end
    end

    % build pseudopopulation from sampled stuff
    samp.trialdat.null = cell(numel(params(1).condition),1);
    samp.trialdat.potent = cell(numel(params(1).condition),1);
    samp.trialdat.shuf.null = cell(numel(params(1).condition),1);
    samp.trialdat.shuf.potent= cell(numel(params(1).condition),1);
    samp.me = cell(numel(params(1).condition),1);
    for ianm = 1:boot.N.anm
        for isess = 1:numel(samp.sessix{ianm})
            sessix = samp.sessix{ianm}(isess);
            bootparams.tmin = params(sessix).tmin;
            bootparams.tmax = params(sessix).tmax;
            bootparams.dt = params(sessix).dt;
            bootparams.cluid = samp.cluid{ianm}{isess};
            bootparams.trialid = samp.trialid{ianm}{isess};
            bootparams.condition = params(sessix).condition;
            bootparams.timeWarp = params(sessix).timeWarp;
            bootparams.smooth = params(sessix).smooth;
            bootparams.bctype = params(sessix).bctype;


            trialdat.null = rez(sessix).recon.null(:,:,samp.cluid{ianm}{isess});
            trialdat.potent = rez(sessix).recon.potent(:,:,samp.cluid{ianm}{isess});
            trialdat.me = me(sessix).data;
            trialid = bootparams.trialid;
            for icond = 1:numel(trialid) % -4 for aw only
                % control
                if numel(trialid{icond}) == 0
                    if ismember(icond,[4 5 8 9 15 16])
                        nt = 20;
                    else
                        nt = 50;
                    end
                    tocat = nan(245,20,80);
                    samp.trialdat.null{icond} = cat(3,samp.trialdat.null{icond},tocat);
                    samp.trialdat.potent{icond} = cat(3,samp.trialdat.potent{icond},tocat);
                else
                    samp.trialdat.null{icond} = cat(3,samp.trialdat.null{icond},trialdat.null(:,trialid{icond},:));
                    samp.trialdat.potent{icond} = cat(3,samp.trialdat.potent{icond},trialdat.potent(:,trialid{icond},:));
                end


                samp.me{icond} = cat(2,samp.me{icond},trialdat.me(:,trialid{icond}));


                %                 if icond == 2 || icond == 3 % create a shuffled trial type CD for hits
                %                     nTrials2SampleShuf = numel(trialid{icond});
                % %                     trixshuf = cell2mat(trialid(2:3)');
                %                     trixshuf = trialid{1};
                %                     trixshuf = randsample(trixshuf,nTrials2SampleShuf,true);
                %                 elseif icond == 4 || icond == 5 % create a shuffled trial type CD for misses
                %                     nTrials2SampleShuf = numel(trialid{icond});
                %                     trixshuf = cell2mat(trialid(4:5)');
                %                     trixshuf = randsample(trixshuf,nTrials2SampleShuf,true);
                %                 else
                %                     trixshuf = trialid{icond};
                %                 end
                %
                %                 % shuffled trial type labels
                %                 samp.trialdat.shuf.null{icond} = cat(3,samp.trialdat.shuf.null{icond},trialdat.null(:,trixshuf,:));
                %                 samp.trialdat.shuf.potent{icond} = cat(3,samp.trialdat.shuf.potent{icond},trialdat.potent(:,trixshuf,:));

            end
        end
    end

    for icond = 1:numel(samp.trialdat.null) % -4 for aw only
        try
            samp.psth.null(:,:,icond) = squeeze(mean(samp.trialdat.null{icond},2));
            samp.psth.potent(:,:,icond) = squeeze(mean(samp.trialdat.potent{icond},2));

            % samp.psth.shuf.null(:,:,icond) = squeeze(mean(samp.trialdat.shuf.null{icond},2));
            % samp.psth.shuf.potent(:,:,icond) = squeeze(mean(samp.trialdat.shuf.potent{icond},2));

            samp.meavg(:,icond) = mean(samp.me{icond},2);
        catch
            
        end


        
    end
    me_(iboot).trialdat = samp.me;
    me_(iboot).conddat  = samp.meavg;

    % coding dimensions

    cond2use = [2 3]; % right hits, left hits
    cond2proj = [2:5]; % right hits, left hits, right miss, left miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
    rampcond = 12; %

    cd_null(iboot) = getCodingDimensions(samp.psth.null,nan,nan,obj(1),params(1),cond2use,cond2use_trialdat, cond2proj, rampcond);
    cd_potent(iboot) = getCodingDimensions(samp.psth.potent,nan,nan,obj(1),params(1),cond2use,cond2use_trialdat, cond2proj, rampcond);
    % project single trials
    for icond = 1:numel(samp.trialdat.null)-4 % -4 for aw only
        cd_null(iboot).trialdat{icond} = tensorprod(samp.trialdat.null{icond},cd_null(iboot).cd_mode_orth,3,1);
        cd_potent(iboot).trialdat{icond} = tensorprod(samp.trialdat.potent{icond},cd_null(iboot).cd_mode_orth,3,1);
    end

    %     cd_null_shuf(iboot) = getCodingDimensions(samp.psth.shuf.null,nan,nan,obj(1),params(1),cond2use,cond2use_trialdat, cond2proj, rampcond);
    %     cd_potent_shuf(iboot) = getCodingDimensions(samp.psth.shuf.potent,nan,nan,obj(1),params(1),cond2use,cond2use_trialdat, cond2proj, rampcond);


    clear samp


end


cd_null_all = concatRezAcrossSessions(cd_null); % in *_bootstrap.m this is concatenating rez across bootstrap iterations (each iteration is a pseudo-session)
cd_potent_all = concatRezAcrossSessions(cd_potent);

% cd_null_all_shuf = concatRezAcrossSessions(cd_null_shuf); % in *_bootstrap.m this is concatenating rez across bootstrap iterations (each iteration is a pseudo-session)
% cd_potent_all_shuf = concatRezAcrossSessions(cd_potent_shuf);

me_all = zeros([size(me_(1).conddat),boot.iters]);
for iboot = 1:boot.iters
    me_all(:,:,iboot) = me_(iboot).conddat;
end


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
% plotCD_EpochSelectivity(meta,obj,cd_null_all_shuf,cd_null_shuf,cd_potent_all_shuf,cd_potent_shuf);


%% plot bootstrapped ME
close all

f = figure;
ax = gca;
hold on;

tstart = mode(obj(1).bp.ev.bitStart) - 2.5;
sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;
gc = 0;

lw = 1;
alph = 0.15;

c = getColors;

temp = squeeze(me_all(:,2,:));
mu = mean(temp,2);
stderr = std(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.rhit,'LineWidth',lw},alph,ax)

temp = squeeze(me_all(:,3,:));
mu = mean(temp,2);
stderr = std(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.lhit,'LineWidth',lw},alph,ax)

% f = figure;
% ax = gca;
% hold on;

temp = squeeze(me_all(:,4,:));
mu = nanmean(temp,2);
stderr = nanstd(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.rmiss,'LineWidth',lw,'LineStyle','--'},alph,ax)

temp = squeeze(me_all(:,5,:));
mu = mean(temp,2);
stderr = std(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.lmiss,'LineWidth',lw,'LineStyle','--'},alph,ax)
xlim([-2.3 2])

xlabel('Time from go cue (s)')
ylabel('Motion energy (a.u.)')
xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--'); xline(gc,'k--')

%% variability in me and CDs

cond2use = [2 3];
cdix = 1; % cd trial type

for iboot = 1:boot.iters

    nulldat = cd_null(iboot).trialdat(cond2use);
    potentdat = cd_potent(iboot).trialdat(cond2use);
    medat = me_(iboot).trialdat(cond2use);

    for icond = 1:numel(cond2use)
        v.null(:,icond,iboot) = squeeze(var(nulldat{icond}(:,:,cdix),[],2));
        v.potent(:,icond,iboot) = squeeze(var(potentdat{icond}(:,:,cdix),[],2));
        v.me(:,icond,iboot) = squeeze(var(medat{icond}(:,:,cdix),[],2));
    end
end

%
f = figure;

xlims = [-2.3 0];

tstart = mode(obj(1).bp.ev.bitStart) - 2.5;
sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;

lw = 1.5;
alph = 0.2;
c = getColors;

ax = nexttile;
hold on;

temp = squeeze(v.null(:,1,:));
mu = mean(temp,2);
stderr = std(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.rmiss,'LineWidth',lw},alph,ax);
temp = squeeze(v.null(:,2,:));
mu = mean(temp,2);
stderr = std(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.lmiss,'LineWidth',lw},alph,ax);
xlim(xlims)
xlabel('Time from go cue (s)')
ylabel('Variance (a.u.)')
title('null cd trial type')
xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--')

ax = nexttile;
hold on;

temp = squeeze(v.potent(:,1,:));
mu = mean(temp,2);
stderr = std(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.rhit_aw,'LineWidth',lw},alph,ax);
temp = squeeze(v.potent(:,2,:));
mu = mean(temp,2);
stderr = std(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.lhit_aw,'LineWidth',lw},alph,ax);
xlim(xlims)
xlabel('Time from go cue (s)')
ylabel('Variance (a.u.)')
title('potent cd trial type')
xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--')

ax = nexttile;
hold on;

temp = squeeze(v.me(:,1,:));
mu = mean(temp,2);
stderr = std(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.rhit,'LineWidth',lw},alph,ax);
temp = squeeze(v.me(:,2,:));
mu = mean(temp,2);
stderr = std(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.lhit,'LineWidth',lw},alph,ax);
xlim(xlims)
xlabel('Time from go cue (s)')
ylabel('Variance (a.u.)')
title('motion energy')
xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--')


for i = 1:numel(xlims)
    [~,ix(i)] = min(abs(obj(1).time - xlims(i)));
end

temp1 = squeeze(v.null(ix(1):ix(2),1,:));
temp2 = squeeze(v.me(ix(1):ix(2),1,:));
cc = corrcoef(temp1(:),temp2(:));
v.cc.null(1) = cc(1,2);
temp1 = squeeze(v.null(ix(1):ix(2),2,:));
temp2 = squeeze(v.me(ix(1):ix(2),2,:));
cc = corrcoef(temp1(:),temp2(:));
v.cc.null(2) = cc(1,2);


tempnull = squeeze(v.null(ix(1):ix(2),1,:));
temppotent = squeeze(v.potent(ix(1):ix(2),1,:));
tempme = squeeze(v.me(ix(1):ix(2),1,:)); 
for i = 1:boot.iters
    v.cc.null(i,1) = corr(tempnull(:,i),tempme(:,i));
    v.cc.potent(i,1) = corr(temppotent(:,i),tempme(:,i));
end
tempnull = squeeze(v.null(ix(1):ix(2),2,:));
temppotent = squeeze(v.potent(ix(1):ix(2),2,:));
tempme = squeeze(v.me(ix(1):ix(2),2,:)); 
for i = 1:boot.iters
    v.cc.null(i,2) = corr(tempnull(:,i),tempme(:,i));
    v.cc.potent(i,2) = corr(temppotent(:,i),tempme(:,i));
end

ax = nexttile;
hold on;
temp = v.cc.null(:,1);
bar(1,mean(temp),'edgecolor','none','facecolor',c.rmiss)
errorbar(1,mean(temp),std(temp),std(temp),'LineStyle','none','Color','k','LineWidth',1);

temp = v.cc.null(:,2);
bar(2,mean(temp),'edgecolor','none','facecolor',c.lmiss)
errorbar(2,mean(temp),std(temp),std(temp),'LineStyle','none','Color','k','LineWidth',1);

temp = v.cc.potent(:,1);
bar(4,mean(temp),'edgecolor','none','facecolor',c.rhit_aw)
errorbar(4,mean(temp),std(temp),std(temp),'LineStyle','none','Color','k','LineWidth',1);

temp = v.cc.potent(:,2);
bar(5,mean(temp),'edgecolor','none','facecolor',c.lhit_aw)
errorbar(5,mean(temp),std(temp),std(temp),'LineStyle','none','Color','k','LineWidth',1);

ax.XTick = [1.5 4.5];
ax.XTickLabel{1} = 'Null';
ax.XTickLabel{2} = 'Potent';
ylabel('Correlation')

sgtitle('Correlation of variance b/w ME and CD TrialType')









