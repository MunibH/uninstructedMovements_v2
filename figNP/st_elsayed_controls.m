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
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off (2)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off (3)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off (4)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off (5)
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
params.tmax = 2.1;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;
params.bctype = 'reflect';

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
% params.quality = {'Excellent','Great','Good','Fair','Multi'};
params.quality = {'Excellent','Great','Good','Fair'};

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
% % % % meta = loadEKH1_ALMVideo(meta,datapth); % selectivity in ME
meta = loadEKH3_ALMVideo(meta,datapth); % selectivity in ME
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB13_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth); % selectivity in ME % go cue is at 2.3 instead of 2.5 like all other sessions??
% meta = loadJEB15_ALMVideo(meta,datapth);
% meta = loadJEB19_ALMVideo(meta,datapth);

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
    disp(['Loading kinematics - Session ' num2str(sessix) '/' num2str(numel(meta))])
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
    responseOnly = 0; % only use response period
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj, nullalltime, onlyAW, delayOnly, responseOnly);

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
% axnull = plotCDProj(cd_null_all,cd_null,sav,titlestring,plotmiss,plotno,cd_null,obj);
% plotCDVarExp(cd_null_all,sav,titlestring)
% plotSelectivity(cd_null_all,cd_null,sav,titlestring)
% plotSelectivityExplained(cd_null_all,cd_null,sav,titlestring)

titlestring = 'Potent';
% axpotent = plotCDProj(cd_potent_all,cd_potent,sav,titlestring,plotmiss,plotno,cd_null,obj);
% plotCDVarExp(cd_potent_all,sav,titlestring)
% plotSelectivity(cd_potent_all,cd_potent,sav,titlestring)
% plotSelectivityExplained(cd_potent_all,cd_potent,sav,titlestring)

% % % match y axis limits across null and potent space CDs
% for i = 1:numel(axnull)
%     ys = [ min(axnull(i).YLim(1),axpotent(i).YLim(1)) max(axnull(i).YLim(2),axpotent(i).YLim(2)) ];
%     axnull(i).YLim = ys;
%     axpotent(i).YLim = ys;
% end

titlestring = 'Null | Potent CDs';
plotCDProj_NP(cd_potent_all,cd_null_all,cd_potent,cd_null,sav,titlestring,plotmiss,obj)
% plotSelectivityExplained_NP(cd_potent_all,cd_null_all,cd_potent,cd_null,sav,titlestring)


%% generate plots

close all

% plotVarianceExplained_DelayResponse(rez,meta);

cond2use = 2:3;
% [nmag,pmag] = 
plotME_NPMagnitude_singleTrials(meta,obj,params,me,rez,cond2use); % sessions 18 and 3
% corr_ME_NPMagnitude(meta,obj,params,me,rez);


% NP_Magnitude_PSTH_bySession_v2(obj,params,rez,meta)

% doAndPlotXCorr_ME_NP

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

sm = 31;
smtype = 'zeropad';

edges = [-0.5 0];
for i = 1:numel(edges)
    [~,ix(i)] = min(abs(obj(1).time - edges(i)));
end

cond2use = [15 16];
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

xx = xline(tstart,'k--'); uistack(xx,'bottom')
xx = xline(sample,'k--');  uistack(xx,'bottom')
xx = xline(delay,'k--');  uistack(xx,'bottom')
xx = xline(gc,'k--'); uistack(xx,'bottom')
xlabel('Time from go cue (s)')
ylabel('Selectivity (spikes / s)')




%% cds_dr_wc.mat vs cds_wc

temp = load('cds_dr_wc.mat');
cds_dr_wc = temp.cds;
temp = load('cds_wc.mat');
cds_wc = temp.cds;

clear temp

labels = {'TrialType','Go','Ramping'};
cond2use = [1 2];
for i = 1:size(cds_dr_wc.null,3)
    temp.dr = squeeze(cds_dr_wc.potent(:,cond2use,i,:)); % (time,cond,session)
    temp.dr = reshape(temp.dr,size(temp.dr,1)*size(temp.dr,2),size(temp.dr,3));
    temp.dr = temp.dr(:,[1:7,end-3:end]); % only AW sessions
    temp.wc = squeeze(cds_wc.potent(:,cond2use,i,:)); % (time,cond,session)
    temp.wc = reshape(temp.wc,size(temp.wc,1)*size(temp.wc,2),size(temp.wc,3));

    cc_ = corr(temp.dr,temp.wc);
    cc{i} = diag(cc_);

end

%%

close all

f=figure;
f.Position = [661   463   266   273];
ax = gca;
hold on;
div = 1;
xs = 1:numel(labels);
for i = 1:numel(xs)
    ci = cc{i};
    h(i) = bar(xs(i),mean(ci));
    h(i).FaceColor = [0.5 0.5 0.5];
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 1;
    h(i).BarWidth = 0.7;
    scatter(xs(i)*ones(size(ci)),ci,25,'MarkerFaceColor','k', ...
            'MarkerEdgeColor','w','LineWidth',1,'XJitter','randn','XJitterWidth',.15, ...
            'MarkerFaceAlpha',1)
    errorbar(h(i).XEndPoints,mean(ci),std(ci)./sqrt(numel(ci)),'LineStyle','none','Color','k','LineWidth',1);
end

% ylim([-0.001 ax.YLim(2)])
ax.XTick = xs;
xticklabels({'TrialType','Go','Ramping'})
ylabel('Corr(DR-WC , WC)')
ax.FontSize = 10;








