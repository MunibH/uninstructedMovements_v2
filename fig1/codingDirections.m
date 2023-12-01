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
params.condition(1)     = {'(hit|miss|no)'};                                              % all trials       (1)
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};      % right hits, 2afc (2)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};      % left hit, 2afc   (3)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};     % right miss, 2afc (4)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};     % left miss, 2afc  (5)
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};        % 2afc hits        (6)
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};         % aw hits          (7)
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};       % right hits, aw   (8)
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};       % left hits, aw    (9)
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};        % ramping          (10)

% params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % left hits, aw    (9)


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
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

%% POPULATION SELECTIVITY CORRELATION MATRIX

cond2use = [2 3]; % left hit, right hit
sel_corr_mat = getSelectivityCorrelationMatrix(obj,cond2use);

%% CODING DIMENSIONS

clearvars -except obj meta params sel_corr_mat

% % 2afc (early, late, go)
cond2use = [2 3]; % left hit, right hit
cond2proj = [2 3 4 5 6 7 8 9];      
rampcond = 10;
rez_2afc = getCodingDimensions_2afc(obj,params,cond2use,cond2proj,rampcond);

% % aw (context mode)
cond2use = [6 7]; % hit 2afc, hit aw
cond2proj = [2 3 4 5 6 7 8 9];
rez_aw = getCodingDimensions_aw(obj,params,cond2use,cond2proj);


allrez = concatRezAcrossSessions(rez_2afc,rez_aw);

%% PLOTS
close all

sav = 0; % 1=save, 0=no_save

% plotSelectivityCorrMatrix(obj(1),sel_corr_mat,params(1).alignEvent,sav)

plotmiss = 0;
plotaw = 0;
plotCDProj(allrez,obj(1),sav,plotmiss,plotaw,params(1).alignEvent)

% plotCDVarExp(allrez,sav)
% plotSelectivity(allrez,rez,sav)
% plotSelectivityExplained(allrez,rez,sav)

% plotCDContext_singleTrials(obj, params, rez_aw, rez_2afc)

%% plot avg trace for each session of cd choice
close all

cdix = 1; % cdchoice
cond = [1 2]; % right hit, left hit (DR)
dat = squeeze(allrez.cd_proj(:,cond,cdix,:)); % (time,cond,session)
datsel = squeeze(dat(:,1,:) - dat(:,2,:));


% cols = getColors;
% c{1} = cols.rhit;
% c{2} = cols.lhit;

darkRed = [110, 0, 2]./255;
lightRed = [255, 115, 117]./255;
darkBlue = [4, 13, 89]./255;
lightBlue = [85, 74, 255]./255;
cmap{1} = createcolormap(numel(meta), darkBlue, lightBlue);
cmap{2} = createcolormap(numel(meta), darkRed, lightRed);

sm = 11;

baselineix = 1:30;


% % projs
% f = figure;
% ax = gca;
% ax = prettifyPlot(ax);
% hold on;
% for isess = 1:numel(meta)
%     % ax = nexttile;
%     % ax = prettifyPlot(ax);
%     % hold on;
%     for icond = 1:numel(cond)
%         toplot = mySmooth(squeeze(dat(:,icond,isess)),sm,'reflect');
%         % toplot = normalize(toplot,'range',[0 1]);
%         toplot = toplot - mean(toplot(baselineix));
%         plot(obj(1).time,toplot,'Color',cmap{icond}(isess,:))
%         % patchline(obj(1).time,toplot,'EdgeColor',c{icond},'EdgeAlpha',0.5)
%     end
% end


ix = findTimeIX(obj(1).time,[-1 0]);
ix = ix(1):ix(2);
selmeans = mean(datsel(ix,:));
[~,sortix] = sort(selmeans);
datsel = datsel(:,sortix);

c1 = [0 0 0]./255;
c2 = [0.8 0.8 0.8];
cmap_ = flip(createcolormap(numel(meta), c1, c2));

% selectivity 
f = figure;
f.Position = [698   436   343   230];
f.Renderer = 'painters';
ax = gca;
ax = prettifyPlot(ax);
hold on;
for isess = 1:numel(meta)
    % ax = nexttile;
    % ax = prettifyPlot(ax);
    % hold on;
    toplot = mySmooth(datsel(:,isess),sm,'reflect');
    % toplot = normalize(toplot,'range',[0 1]);
    toplot = toplot - mean(toplot(baselineix));
    plot(obj(1).time,toplot,'Color',cmap_(isess,:),'LineWidth',1.5)
end
yline(0,'k-')
plot(obj(1).time,mean(datsel,2),'Color','r')
xlabel('Time from go cue (s)')
ylabel('Selectivity (spks/s)')
evtimes = getEventTimes(obj(1).bp.ev,{'bitStart','sample','delay','goCue'},params(1).alignEvent);
plotEventTimes(ax,evtimes)
xlim([-2.4;2])

%%


alph = 0.2;
cols = getColors;


f = figure;
f.Renderer = 'painters';
ax = gca;
ax = prettifyPlot(ax);
hold on;
icd = 3;
cond = [1 2]; 
lw = 1.5;
ls = '-';
temp = squeeze(allrez.cd_proj(:,cond,icd,:)) * -1 - 7;
mu = nanmean(temp,3);
sd = nanstd(temp,[],3) ./ sqrt(numel(obj)) / 2;
shadedErrorBar(obj(1).time,mu(:,1),sd(:,1),{'Color',cols.rhit,'LineWidth',lw,'LineStyle',ls},alph,ax)
shadedErrorBar(obj(1).time,mu(:,2),sd(:,2),{'Color',cols.lhit,'LineWidth',lw,'LineStyle',ls},alph,ax)
cond = [7 8]; 
lw = 1.5;
ls = '--';
temp = squeeze(allrez.cd_proj(:,cond,icd,:)) * -1 - 10;
mu = nanmean(temp,3);
sd = nanstd(temp,[],3) ./ sqrt(numel(obj)) / 2;
shadedErrorBar(obj(1).time,mu(:,1),sd(:,1),{'Color',cols.rhit,'LineWidth',lw,'LineStyle',ls},alph,ax)
shadedErrorBar(obj(1).time,mu(:,2),sd(:,2),{'Color',cols.lhit,'LineWidth',lw,'LineStyle',ls},alph,ax)




