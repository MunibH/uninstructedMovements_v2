clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig3/')))
rmpath(genpath(fullfile(utilspth,'mc_stim/')))
rmpath(genpath(fullfile(utilspth,'MotionMapper/')))


% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                                              % all trials       (1)
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};      % right hits, 2afc (2)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};      % left hit, 2afc   (3)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};     % right miss, 2afc (4)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};     % left miss, 2afc  (5)
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};        % 2afc hits        (6)
params.condition(end+1) = {'hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};         % aw hits          (7)
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % right hits, aw   (8)
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % left hits, aw    (9)


params.tmin = -3;
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
cond2proj = [2 3 4 5 8 9];
rez_2afc = getCodingDimensions_2afc(obj,params,cond2use,cond2proj);

% % aw (context mode)
cond2use = [6 7]; % hit 2afc, hit aw
cond2proj = [2 3 4 5 8 9];
rez_aw = getCodingDimensions_aw(obj,params,cond2use,cond2proj);


allrez = concatRezAcrossSessions(rez_2afc,rez_aw);

%% PLOTS
close all

sav = 0; % 1=save, 0=no_save

% plotSelectivityCorrMatrix(obj(1),sel_corr_mat,params(1).alignEvent,sav)

plotmiss = 0;
plotaw = 1;
plotCDProj(allrez,rez_2afc,sav,plotmiss,plotaw,params(1).alignEvent)
% set(gcf,'Position',[249   558   349   238])
% plotCDVarExp(allrez,sav)
% plotSelectivity(allrez,rez,sav)
% plotSelectivityExplained(allrez,rez,sav)





%% heatmap for context mode
close all

trialStart = mode(obj(1).bp.ev.bitStart - rez_2afc(1).align);
sample = mode(obj(1).bp.ev.sample - rez_2afc(1).align);
for sessix = 5 %1:numel(rez_aw)
%     trix = sort(cell2mat(params(sessix).trialid([6 7])')); % 2afc and aw hits
    trix = 1:obj(sessix).bp.Ntrials;
    trialdat = obj(sessix).trialdat(:,:,trix);
    trialdat = permute(trialdat, [1 3 2]); % time trials clu
    dims = size(trialdat);
    temp = reshape(trialdat, dims(1)*dims(2), dims(3));
    temp = zscore(temp);

    proj = temp * rez_aw(sessix).cd_mode_orth;
    proj = reshape(proj, dims(1), dims(2));

    proj(end+1:end+20,:) = repmat((obj(sessix).bp.autowater * max(max(proj)))', 20,1);
    tm = obj(sessix).time;
    tm = [tm tm(end)+(1:20)*mode(diff(tm))];

    figure; imagesc(1:numel(trix), tm, proj)
    set(gca,'YDir','normal')
    xlabel('Trials')
    ylabel(['Time (s) from ' params(sessix).alignEvent])
    title('Context Mode')
    colormap(flipud(linspecer))
    c = colorbar;
%     c.Limits(1) = c.Limits(1) ./ 3;
%     ylim([trialStart sample])
    xlim([1 obj(sessix).bp.Ntrials - 40])

    ax = gca;
    ax.FontSize = 15;

end




















