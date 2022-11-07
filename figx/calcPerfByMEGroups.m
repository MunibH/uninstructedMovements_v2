clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig2/')))
rmpath(genpath(fullfile(utilspth,'fig1/')))
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
params.condition(1)     = {'(hit|miss)&~stim.enable&~autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-40))'};             % right hits, no stim, aw off


params.tmin = -2.5;
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
% meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

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


%%

% for each session
% sort trials by average motion energy in delay period
% divide trials into N sequential (sorted) groups
% calculate performance for each of those groups

taxis = obj(1).time;
[~,ix1] = min(abs(obj(1).time - -0.5));
[~,ix2] = min(abs(obj(1).time - -0.02));

nGroups = 5;

clear perf
figure; hold on;
for sessix = 1:numel(meta)
    trix = params(sessix).trialid{1};
    %     trix = cell2mat(params(sessix).trialid(1:2)');
    me_ = me(sessix).data(ix1:ix2,trix);
    [~,sortix] = sort(mean(me_,1),'descend');
    trix_sorted = trix(sortix);
    plot(trix_sorted,'.')
    %     temp = me.data(:,trix_sorted);
    %     figure; imagesc(temp')

    nTrials = numel(trix);
    dt = floor(nTrials / nGroups);

    ct = 1;
    ii = 1:dt:nTrials;
    for i = 1:nGroups
        trix_ = trix_sorted(ii(i):ii(i)+dt-1);
        perf(sessix,ct) = sum(obj(sessix).bp.hit(trix_)) ./ numel(trix_);
        ct = ct + 1;
    end
end


%% plot
close all


h = 0;
s = 0;
v = 0.8;
xs = 1:5;
for i = 1:numel(xs)
    cols(i,:) = hsv2rgb([h s v-(0.1)*i]);
end

div = 1.3;

f = figure; hold on;
% f.Position = [680   755   244   223];

xs = [1 2 3 4 5];
for i = 1:numel(xs)
    temp = perf(:,i);
    b(i) = bar(xs(i),nanmean(temp));
    b(i).FaceColor = cols(i,:);
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.8;
    vs(i) = scatter(xs(i)*ones(size(temp)),temp,60,'MarkerFaceColor',cols(i,:)./div,...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25);
    errorbar(b(i).XEndPoints,nanmean(temp),nanstd(temp),'LineStyle','none','Color','k','LineWidth',1)


end

lab = 'ME Group ';
for i = 1:numel(xs)
    labels{i} = [lab num2str(i)];
end

xticklabels([' ' labels])
ylabel("Performance")
% ylim([0,1])
ax = gca;
ax.FontSize = 12;




%%
clc,close all
% want to predict performance up to a given trial given trial number
% then regress that out,
% then see how residual perf changes as a function of ME split into 5 groups

nGroups = 2;
figure; hold on;
for sessix = 1:numel(meta)
    clear perf
    % get running performance
    trix = sort(params(sessix).trialid{1});
    nTrials = numel(trix);
    trialsPerGroup = floor(nTrials / nGroups);
    for i = 1:nTrials
        perf(i) = sum(obj(sessix).bp.hit(trix(1:i))) ./ numel(1:i);
    end
    if perf(1) == 0
        perf(1) = eps;
    end
%     histogram(perf)
%     X = eye(nTrials); % a 1 for the trial number, 0s everywhere else, ordinal variable
%     X = [1:nTrials]';
    X = zeros(nTrials,1);
    X(1:trialsPerGroup) = 1;
    X(trialsPerGroup+1:trialsPerGroup*2) = 2;
%     X(trialsPerGroup*2+1:trialsPerGroup*3) = 3;
%     X(trialsPerGroup*3+1:trialsPerGroup*4) = 4;
%     X(trialsPerGroup*4+1:end) = 5;


    [beta, ~, output] = glmfit(X, perf', 'normal');
    plot(perf,'k')
    plot(output.resid,'r')

    yhat = glmval(beta,X,'identity');
    plot(yhat,'b')

    slope = gradient(yhat);
    plot(slope,'g')

end












