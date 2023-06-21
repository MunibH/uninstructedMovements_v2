clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig3/')))
rmpath(genpath(fullfile(utilspth,'fig2/')))
rmpath(genpath(fullfile(utilspth,'figx/')))
rmpath(genpath(fullfile(utilspth,'mc_stim/')))
rmpath(genpath(fullfile(utilspth,'MotionMapper/')))
rmpath(genpath(fullfile(utilspth,'musall2919/')))

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
params.condition(1)     = {'(hit|miss|no)'};                             % all trials       (1)
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};      % right hits, 2afc (2)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};      % left hit, 2afc   (3)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};     % right miss, 2afc (4)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};     % left miss, 2afc  (5)
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};        % 2afc hits        (6)
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};         % aw hits          (7)
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};       % right hits, aw   (8)
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};       % left hits, aw    (9)
params.condition(end+1) = {'R&no&~stim.enable&~autowater&~early'};       % right no, 2afc   (10)
params.condition(end+1) = {'L&no&~stim.enable&~autowater&~early'};       % left no, 2afc    (11)
% for ramping
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % all hits, no stim, aw off (12)

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;
params.bctype = 'reflect'; % reflect, zeropad, none

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance


params.advance_movement = 0.0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
% meta = loadEKH3_ALMVideo(meta,datapth); % not usable b/c no usable left miss trials
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written



%% LOAD DATA OBJECTS

[obj,params] = loadSessionData(meta,params);

% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
%     kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end



%% BOOTSTRAP PARAMS

clear boot bootobj bootparams

boot.iters = 1000; % number of bootstrap iterations (most papers do 1000)

boot.N.anm = 5; % number of animals to sample w/ replacement
boot.N.sess = 2; % number of sessions to sample w/ replacement (if Nsessions for an animal is less than this number, sample Nsessions)
boot.N.trials.hit = 50; % number of hit trials to sample w/o replacement
boot.N.trials.miss = 20; % number of miss trials to sample w/o replacement
boot.N.clu = 20; % number of neurons to sample w/o replacement

%%

clear samp allrez

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
                    samp.trialid{ianm}{isess}{icond} = randsample(trialid{icond},nTrials2Sample,false);
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

    % randomly sample neurons 
    for ianm = 1:boot.N.anm
        for isess = 1:numel(samp.sessix{ianm})
            sessix = samp.sessix{ianm}(isess);
            cluid = 1:size(obj(sessix).psth,2); % index of cluster in psth,trialdat,etc 
            samp.cluid{ianm}{isess} = randsample(cluid,boot.N.clu,false);
        end
    end

    samp.trialdat = cell(numel(params(1).condition),1);
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
            


            trialdat = obj(sessix).trialdat(:,samp.cluid{ianm}{isess},:);
            me_trialdat = me(sessix).data;
            trialid = bootparams.trialid;
            for icond = 1:numel(trialid)
                samp.trialdat{icond} = cat(3,samp.trialdat{icond},trialdat(:,:,trialid{icond}));
                samp.me{icond} = cat(2,samp.me{icond},me_trialdat(:,trialid{icond}));
            end
        end
    end

    for icond = 1:numel(samp.trialdat)
        samp.psth(:,:,icond) = squeeze(mean(samp.trialdat{icond},3));
        samp.meavg(:,icond) = mean(samp.me{icond},2);
    end
    me_(iboot).trialdat = samp.me;
    me_(iboot).conddat  = samp.meavg;

    bootobj = rmfield(obj(1),'psth');
    bootobj.psth = samp.psth;
    bootparams.alignEvent = params(1).alignEvent;

    % coding dimensions

    % % 2afc (early, late, go)
    cond2use = [2 3]; % left hit, right hit
    cond2proj = [2 3 4 5 6 7 8 9 10 11];
    rampcond = 12;
    rez_2afc = getCodingDimensions_2afc(bootobj,bootparams,cond2use,cond2proj,rampcond);

    % % aw (context mode)
    cond2use = [6 7]; % hit 2afc, hit aw
    cond2proj = [2 3 4 5 6 7 8 9 10 11];
    rez_aw = getCodingDimensions_aw(bootobj,bootparams,cond2use,cond2proj);

    allrez(iboot) = concatRezAcrossSessions(rez_2afc,rez_aw);

    clear samp


end


% CONCAT BOOTSTRAP ITERATIONS

% mean across sessions for each bootstrap iteration
temp = allrez;
for i = 1:numel(temp)
    temp(i).cd_proj = nanmean(temp(i).cd_proj,4);
%     temp(i).cd_varexp = nanmean(temp(i).cd_varexp,1);
%     temp(i).cd_varexp_epoch = nanmean(temp(i).cd_varexp_epoch,1);
%     temp(i).selectivity_squared = nanmean(temp(i).selectivity_squared,3);
%     temp(i).selexp = nanmean(temp(i).selexp,3);
end

% concatenate bootstrap iterations
rez = temp(1);
for i = 2:numel(allrez)
    rez.cd_proj = cat(4,rez.cd_proj,temp(i).cd_proj);
%     rez.cd_varexp = cat(1,rez.cd_varexp,temp(i).cd_varexp);
%     rez.cd_varexp_epoch = cat(1,rez.cd_varexp_epoch,temp(i).cd_varexp_epoch);
%     rez.selectivity_squared = cat(3,rez.selectivity_squared,temp(i).selectivity_squared);
%     rez.selexp = cat(3,rez.selexp,temp(i).selexp);
end

me_all = zeros([size(me_(1).conddat),boot.iters]);
for iboot = 1:boot.iters
    me_all(:,:,iboot) = me_(iboot).conddat;
end



%% PLOTS
close all

sav = 0; % 1=save, 0=no_save

plotmiss = 0;
plotaw = 1;
plotCDProj(rez,obj(1),sav,plotmiss,plotaw,params(1).alignEvent,3)

%%
% center data for go mode
close all
alph = 0.2;
lw = 2;

rampcd = 3;
cond = [1 2];
rampproj = squeeze(rez.cd_proj(:,cond,rampcd,:));
rampproj = squeeze(mean(rampproj,2));
mu = mean(rampproj,2)*-1;
sd = std(rampproj,[],2);
f = figure;
f.Position = [698   436   343   230];
f.Renderer = 'painters';
ax = gca;
ax = prettifyPlot(ax);
shadedErrorBar(obj(1).time,mu,sd,{'Color','k','LineWidth',lw},alph,ax)
% xlim([obj(1).time(5);2])
% 
% title('ramp','FontSize',8)
% xlabel(['Time from go cue (s)'])
% ylabel('Activity (a.u.)')
% ax.FontSize = 10;
% ax = prettifyPlot(ax);
%%

% plotCDVarExp(allrez,sav)
% plotSelectivity(allrez,rez,sav)
% plotSelectivityExplained(allrez,rez,sav)

% plotCDContext_singleTrials(obj, params, rez_aw, rez_2afc)



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

temp = squeeze(me_all(:,10,:));
mu = nanmean(temp,2);
stderr = nanstd(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.rmiss,'LineWidth',lw,'LineStyle','--'},alph,ax)

temp = squeeze(me_all(:,11,:));
mu = mean(temp,2);
stderr = std(temp,[],2);
shadedErrorBar(obj(1).time,mu,stderr,{'Color',c.lmiss,'LineWidth',lw,'LineStyle','--'},alph,ax)
xlim([-2.3 2])

xlabel('Time from go cue (s)')
ylabel('Motion energy (a.u.)')
xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--'); xline(gc,'k--')







