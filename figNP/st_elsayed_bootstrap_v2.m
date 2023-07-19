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
% 2afc
params.condition(1)     = {'R&hit&~stim.enable&~autowater'};                 % right hits, no stim, aw off   (1)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off    (2)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off  (3)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off   (4)

% for projections
% 2afc
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off   (5)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off    (6)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off  (7)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off   (8)


% for ramping
params.condition(end+1) = {'hit&~stim.enable&~autowater'};               % all hits, no stim, aw off (9)

% wc
params.condition(end+1) = {'R&hit&~stim.enable&autowater'};             % right hits, no stim, aw on   (10)
params.condition(end+1) = {'L&hit&~stim.enable&autowater'};             % left hits, no stim, aw on    (11)
% params.condition(end+1) = {'R&miss&~stim.enable&autowater'};            % error right, no stim, aw on  (12)
% params.condition(end+1) = {'L&miss&~stim.enable&autowater'};            % error left, no stim, aw on   (13)

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
% meta = loadJEB6_ALMVideo(meta,datapth);
% meta = loadJEB7_ALMVideo(meta,datapth);
% % % meta = loadEKH1_ALMVideo(meta,datapth); % not usable b/c no usable left miss trials
% meta = loadEKH3_ALMVideo(meta,datapth);
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB13_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);
% meta = loadJEB15_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);

meta = meta(3);

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
    disp(['Loading kinematics ' num2str(sessix) '/' num2str(numel(meta))])
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
% -----------------------------------------------------------------------
disp('finding null and potent spaces')
for sessix = 1:numel(meta)
    disp([meta(sessix).anm ' ' meta(sessix).date])
    % -- input data
    % trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix));
    trialdat_zscored = permute(obj(sessix).trialdat, [1 3 2]);

    % -- null and potent spaces
    cond2use = [1:4 10:11]; % right hit, left hit, right miss, left miss, 2afc
    cond2proj = 1:9; % ~early versions
    nullalltime = 0; % use all time points to estimate null space if 1
    onlyAW = 0; % only use AW trials
    delayOnly = 0; % only use delay period
    responseOnly = 0;
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj, nullalltime, onlyAW, delayOnly, responseOnly);

    % -- coding dimensions
    cond2use = [5 6]; % right hits, left hits (corresponding to null/potent psths in rez)
    cond2proj = [1:9]; % right hits, left hits, right miss, left miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [5 6]; % for calculating selectivity explained in full neural pop
    rampcond = 9; % corresponding to cond2proj in null/potent analysis
    %     cd_null(sessix) = getCodingDimensions(rez(sessix).N_null_psth,rez(sessix).N_null,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,rampcond);
    %     cd_potent(sessix) = getCodingDimensions(rez(sessix).N_potent_psth,rez(sessix).N_potent,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,rampcond);

    cd_null(sessix) = getCodingDimensions(rez(sessix).recon_psth.null,trialdat_zscored,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj, rampcond);
    cd_potent(sessix) = getCodingDimensions(rez(sessix).recon_psth.potent,trialdat_zscored,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj, rampcond);
end
disp('DONE')


%%

% plotReconstructedPSTHs(obj,rez,params)


%% BOOTSTRAP PARAMS

clear boot bootobj bootparams

boot.iters = 1000; % number of bootstrap iterations (most papers do 1000)

boot.N.anm = 5; % number of animals to sample w/ replacement
boot.N.sess = 2; % number of sessions to sample w/ replacement (if Nsessions for an animal is less than this number, sample Nsessions)
boot.N.trials.hit = 50; % number of hit trials to sample
boot.N.trials.miss = 20; % number of miss trials to sample
boot.N.clu = 20; % number of neurons to sample

boot.cd.cond2use = [5 6]; % right hits, left hits
boot.cd.cond2proj = [5:8];
boot.cd.cond2use_trialdat = [5 6]; 
boot.cd.rampcond = 9; 

%%
clear bootrez samp cd_null cd_potent cd_null_shuf cd_potent_shuf cd_null_all cd_potent_all cd_null_all_shuf cd_potent_all_shuf trialdat sel

[objixs,uAnm] = groupSessionsByAnimal(meta);
nConds = numel(params(1).condition);
nSess = boot.N.anm*boot.N.sess;
nClu = nSess * boot.N.clu;

for iboot = 1:boot.iters
    disp(['Iteration ' num2str(iboot) '/' num2str(boot.iters)]);

    % randomly sample animals
    bootrez.anm = randsample(uAnm,boot.N.anm,true);

    % randomly sample sessions
    xs = 1:boot.N.sess:nSess;
    for ianm = 1:boot.N.anm
        objix = find(objixs{ismember(uAnm,bootrez.anm{ianm})});
        bootrez.sessix(xs(ianm):xs(ianm)+1) = randsample(objix,boot.N.sess,true);
    end

    % randomly sample trials
    for isess = 1:numel(bootrez.sessix)
        sessix = bootrez.sessix(isess);
        trialid = params(sessix).trialid;
        for icond = 1:nConds
            cond = params(sessix).condition{icond};
            nTrialsCond = numel(trialid{icond});
            if nTrialsCond == 0
                nTrials2Sample = 0;
                % error(['no trials found for ' meta(sessix).anm ' ' meta(sessix).date ' ' cond])
            else
                if contains(cond,{'hit'})
                    nTrials2Sample = boot.N.trials.hit;
                elseif contains(cond,{'miss'})
                    nTrials2Sample = boot.N.trials.miss;
                end
            end
            
            bootrez.trialid{isess,icond} = randsample(trialid{icond},nTrials2Sample,true);
        end
    end

    % randomly sample neurons (from null and potent reconstructions)
    for isess = 1:numel(bootrez.sessix)
        sessix = bootrez.sessix(isess);
        cluid = 1:size(obj(sessix).psth,2); % index of cluster in psth,trialdat,etc -> different from bootstrapping in fig1 b/c we don't need to re-do getSeq, just going to grab data from rez(sessix)
        bootrez.cluid(isess,:) = randsample(cluid,boot.N.clu,true); % (sessions,clusters)
    end

    % get single trial data and PSTHs based off randomly sampled parameters
    [bootrez.trialdat.null, bootrez.trialdat.potent] = deal(cell(nConds,1));
    for i = 1:nConds
        nTrialsCond = numel(bootrez.trialid{1,i}) * nSess;
        [bootrez.trialdat.null{i}, ...
        bootrez.trialdat.potent{i},...
        bootrez.trialdat.full{i}] = deal(nan(numel(obj(1).time), nTrialsCond , nClu));
    end
    [bootrez.psth.null, bootrez.psth.potent, bootrez.psth.full] = deal(nan(numel(obj(1).time), nClu , nConds));
    ct = 1;
    tct = 1;
    for isess = 1:numel(bootrez.sessix)
        sessix = bootrez.sessix(isess);
        cluid = bootrez.cluid(isess,:);
        for icond = 1:nConds
            trix = bootrez.trialid{isess,icond};
            nTrialsCond = numel(trix);

            % temp = permute(obj(sessix).trialdat(:,cluid,trix), [1 3 2]);
            % dims = size(temp); % (time,trials,neurons)
            % temp = reshape(temp,dims(1)*dims(2),dims(3));
            % temp2 = zscore(temp);
            % trialdat.full = reshape(temp2,dims(1),dims(2),dims(3));
            trialdat.full = permute(obj(sessix).trialdat(:,cluid,trix), [1 3 2]);
            bootrez.psth.full(:,ct:ct+(nClu/nSess)-1,icond) = squeeze(mean(trialdat.full,2));

            trialdat.null = rez(sessix).recon.null(:,trix,cluid);
            bootrez.psth.null(:,ct:ct+(nClu/nSess)-1,icond) = squeeze(mean(trialdat.null,2));
            trialdat.potent = rez(sessix).recon.potent(:,trix,cluid);
            bootrez.psth.potent(:,ct:ct+(nClu/nSess)-1,icond) = squeeze(mean(trialdat.potent,2));

            bootrez.trialdat.full{i}(:,tct:tct+nTrialsCond-1,ct:ct+(nClu/nSess)-1) = trialdat.full;
            bootrez.trialdat.null{i}(:,tct:tct+nTrialsCond-1,ct:ct+(nClu/nSess)-1) = trialdat.null;
            bootrez.trialdat.potent{i}(:,tct:tct+nTrialsCond-1,ct:ct+(nClu/nSess)-1) = trialdat.potent;

            
        end
        ct = ct + (nClu/nSess);
    end

    % coding dimensions
    % cd_null(iboot) = getCodingDimensions(bootrez.psth.null,nan,nan,obj(1),params(1),boot.cd.cond2use,boot.cd.cond2use_trialdat, boot.cd.cond2proj, boot.cd.rampcond);
    % cd_potent(iboot) = getCodingDimensions(bootrez.psth.potent,nan,nan,obj(1),params(1),boot.cd.cond2use,boot.cd.cond2use_trialdat, boot.cd.cond2proj, boot.cd.rampcond);

    % find coding dimensions from full population, project null/potent data on these CDs
    % % 2afc (early, late, go)
    cond2use = [5 6]; % right hit, left hit (~early)
    cond2proj = [1:9];
    rampcond = 9;
    rez_2afc = getCodingDimensions_boot(obj(1),params(1),bootrez.psth.full,cond2use,cond2proj,rampcond);

    cd_full(iboot).cd_proj = permute(tensorprod(bootrez.psth.full,rez_2afc.cd_mode_orth,2,1),[1 3 2]); % (time,cd,cond)
    cd_null(iboot).cd_proj = permute(tensorprod(bootrez.psth.null,rez_2afc.cd_mode_orth,2,1),[1 3 2]); % (time,cd,cond)
    cd_potent(iboot).cd_proj = permute(tensorprod(bootrez.psth.potent,rez_2afc.cd_mode_orth,2,1),[1 3 2]); % (time,cd,cond)

    % project single trials onto CDs
    for icond = 1:nConds
        % cd_null(iboot).trialdat{icond} = tensorprod(bootrez.trialdat.null{icond},cd_null(iboot).cd_mode_orth,3,1);
        % cd_potent(iboot).trialdat{icond} = tensorprod(bootrez.trialdat.potent{icond},cd_null(iboot).cd_mode_orth,3,1);

        % if using full pop cds for projs
        cd_null(iboot).trialdat{icond} = tensorprod(bootrez.trialdat.null{icond},rez_2afc.cd_mode_orth,3,1);
        cd_potent(iboot).trialdat{icond} = tensorprod(bootrez.trialdat.potent{icond},rez_2afc.cd_mode_orth,3,1);
    end

    % %% latency to selectivity onset in cd choice
    % cdix = 1;
    % itiix = [-2.39 -2.2];
    % ix = findTimeIX(obj(1).time,itiix);
    % ix = ix(1):ix(2);
    % null1 = squeeze(cd_null(iboot).cd_proj(:,cdix,5));
    % null2 = squeeze(cd_null(iboot).cd_proj(:,cdix,6));
    % sel.null = null1 - null2;
    % iti = mean(sel.null(ix));
    % sel.null = sel.null - iti;
    % itisd = std(sel.null(ix));
    % thresh = itisd * 3; % when selectivity increases baseline + 0.3*iti_stdev (already baseline subtracted tho)
    % % find when selectivity emerges post-sample only
    % % to do this, i'm finding last time selectivity is not ever below
    % % threshold
    % sampleix = [-2.15 0]; % 50 ms after sample start
    % ix = findTimeIX(obj(1).time,sampleix);
    % ix = ix(1):ix(2);
    % temp = flip(sel.null(ix));
    % onsetix = find(temp<=(thresh),1,'first');
    % onsetix = numel(temp) - onsetix; % from sample onset
    % onsetix = onsetix + ix(1); % add back presample
    % sel.ix.null(iboot) = onsetix;
    % 
    % 
    % 
    % 
    % % f = figure;
    % % ax = gca;
    % % hold on;
    % % plot(sel.null,'color',cols.null,'linewidth',2)
    % % xline(sel.ix.null(iboot))
    % 
    % cdix = 1;
    % itiix = [-2.39 -2.2];
    % ix = findTimeIX(obj(1).time,itiix);
    % ix = ix(1):ix(2);
    % potent1 = squeeze(cd_potent(iboot).cd_proj(:,cdix,5));
    % potent2 = squeeze(cd_potent(iboot).cd_proj(:,cdix,6));
    % sel.potent = potent1 - potent2;
    % iti = mean(sel.potent(ix));
    % sel.potent = sel.potent - iti;
    % itisd = std(sel.potent(ix));
    % thresh = itisd * 3; % when selectivity increases baseline + 0.3*iti_stdev (already baseline subtracted tho)
    % % find when selectivity emerges post-sample only
    % % to do this, i'm finding last time selectivity is not ever below
    % % threshold
    % sampleix = [-2.15 0]; % 50 ms after sample start
    % ix = findTimeIX(obj(1).time,sampleix);
    % ix = ix(1):ix(2);
    % temp = flip(sel.potent(ix));
    % onsetix = find(temp<=(thresh),1,'first');
    % onsetix = numel(temp) - onsetix; % from sample onset
    % onsetix = onsetix + ix(1); % add back presample
    % sel.ix.potent(iboot) = onsetix;

    
    % figure(1); 
    % plot(sel.null); hold on; plot(sel.potent)
    % pause
    % clf(figure(1))

    clear bootrez

end
cd_full_all = concatRezAcrossSessions(cd_full);
cd_null_all = concatRezAcrossSessions(cd_null); % in *_bootstrap.m this is concatenating rez across bootstrap iterations (each iteration is a pseudo-session)
cd_potent_all = concatRezAcrossSessions(cd_potent);


%% plot coding dimensions from activity N/P spaces
close all

sav = 0;

plotmiss = 0;
plotno = 0;

titlestring = 'Null';
axnull = plotCDProj_v2(cd_null_all,cd_null,sav,titlestring,plotmiss,plotno,rez_2afc,obj);
% plotCDVarExp(cd_null_all,sav,titlestring)
% plotSelectivity(cd_null_all,cd_null,sav,titlestring)
% plotSelectivityExplained(cd_null_all,cd_null,sav,titlestring)

titlestring = 'Potent';
axpotent = plotCDProj_v2(cd_potent_all,cd_potent,sav,titlestring,plotmiss,plotno,rez_2afc,obj);
% plotCDVarExp(cd_potent_all,sav,titlestring)
% plotSelectivity(cd_potent_all,cd_potent,sav,titlestring)
% plotSelectivityExplained(cd_potent_all,cd_potent,sav,titlestring)

% titlestring = 'Full';
% axpotent = plotCDProj_v2(cd_full_all,cd_full,sav,titlestring,plotmiss,plotno,rez_2afc,obj);
% plotCDVarExp(cd_potent_all,sav,titlestring)
% plotSelectivity(cd_potent_all,cd_potent,sav,titlestring)
% plotSelectivityExplained(cd_potent_all,cd_potent,sav,titlestring)

% % match y axis limits across null and potent space CDs
for i = 1:numel(axnull)
    ys = [ min(axnull(i).YLim(1),axpotent(i).YLim(1)) max(axnull(i).YLim(2),axpotent(i).YLim(2)) ];
    axnull(i).YLim = ys;
    axpotent(i).YLim = ys;
end

% % titlestring = 'Null | Potent CDs';
% % plotCDProj_NP(cd_potent_all,cd_null_all,cd_potent,cd_null,sav,titlestring,plotmiss)
% % plotSelectivityExplained_NP(cd_potent_all,cd_null_all,cd_potent,cd_null,sav,titlestring)


%% selectivity in various epochs in cd late hits and errors
close all
clc

plotCDSelectivity_v2(meta,obj,cd_null_all,cd_null,cd_potent_all,cd_potent);

% plotCD_EpochSelectivity(meta,obj,cd_null_all,cd_null,cd_potent_all,cd_potent);
% plotCD_EpochSelectivity(meta,obj,cd_null_all_shuf,cd_null_shuf,cd_potent_all_shuf,cd_potent_shuf);



%% latency plots

close all

nix = sel.ix.null;
pix = sel.ix.potent;
nix([find(isnan(nix)) find(isnan(pix))]) = [];
pix([find(isnan(pix)) find(isnan(pix))]) = [];

% f = figure;
% ax = gca;
% f.Renderer = 'painters';
% ax = prettifyPlot(ax);
% hold on;
% nbins = 50;
% histogram(obj(1).time(nix)+2.5,nbins,'linewidth',2,'edgecolor',cols.null,'facecolor','none','Normalization','cdf','displaystyle','stairs')
% histogram(obj(1).time(pix)+2.5,nbins,'linewidth',2,'edgecolor',cols.potent,'facecolor','none','Normalization','cdf','displaystyle','stairs')
% % histogram(log(sel.ix.null),nbins,'edgecolor','none','facecolor',cols.null)
% % histogram(log(sel.ix.potent),nbins,'edgecolor','none','facecolor',cols.potent)
% xlabel('time from sample (s)')
% ylabel('CDF')



f = figure;
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;
Y = [obj(1).time(nix)'+2.5 obj(1).time(pix)'+2.5];
Z = Y(:,1) - Y(:,2); % time difference in selectivity onset (null - potent) +ve means potent first, -ve means null first
nullZ = Z(Z<=0);
potentZ = Z(Z>=0);
h = histogram(Z,nbins,'edgecolor','none','Normalization','probability','Visible','off');
bars = h.Values;
binedges = h.BinEdges;

x = find(binedges<=0);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.null;

x = find(binedges>=0);
x = x(1:end-1);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.potent;


prob.nullfirst = sum(Z<0) ./ numel(nix);
prob.potentfirst = sum(Z>0) ./ numel(nix);
disp(['Null first: ' num2str(prob.nullfirst)])
disp(['Potent first: ' num2str(prob.potentfirst)])
xline(0,'k--')


% [h,p,ci,stats] = my_ttest(x,m,varargin)
[h,p,ci,stats] = my_ttest(Y(:,1),Y(:,2));
title(['paired t-test, p = ' num2str(p)])
xlabel('Difference in selectivity onset (s)')
ylabel('Probability')






