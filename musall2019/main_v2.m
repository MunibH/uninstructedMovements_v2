clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
% remove paths we don't need
rmpath(genpath(fullfile(utilspth,'fig1/')));
rmpath(genpath(fullfile(utilspth,'mc_stim/')));
rmpath(genpath(fullfile(utilspth,'fig3/')));
rmpath(genpath(fullfile(utilspth,'fig2/')));
rmpath(genpath(fullfile(utilspth,'MotionMapper/')));
rmpath(genpath(fullfile(utilspth,'figx/')));

% add paths for figure specific functions
addpath(genpath(pwd))

clc

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'R&(hit|miss)&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&(hit|miss)&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
params.condition(end+1) = {'(hit|miss)&~stim.enable&~autowater'};

params.tmin = -1.3;
params.tmax = 1;
params.dt = 1/200;
params.viddt = 0.0025;

% smooth with causal gaussian kernel
params.smooth = 31;
params.bctype = 'reflect';

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% svd params
params.reload_Cam0 = true; % redo video svd for Cam0
params.reload_Cam1 = true; % redo video svd for Cam1
params.nrDims = 200; % number of SVs to keep, 200 is probably more than needed, but that's ok

% regression params
params.frameRate = 1/params.dt / 10;
params.sPostTime = ceil(6 * params.frameRate);   % follow stim events for sPostStim in frames (used for eventType 2)
params.mPreTime = ceil(0.5 * params.frameRate);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
params.mPostTime = ceil(1 * params.frameRate);   % follow motor events for mPostStim in frames (used for eventType 3)
params.folds = 10;                               % number of cross validation folds


%% SPECIFY DATA TO LOAD

% --------------------------------------------------------- %
% TODO: make this code work for multiple sessions
% this code only works for 1 session at the moment
% --------------------------------------------------------- %
objpth = 'C:\Users\munib\Documents\Economo-Lab\data\';
vidpth = 'D:\Economo-Lab\data\Video\';

meta = [];

% jeb6 4-18
% jeb7 4-29
% jeb15 7-26

% --- ALM ---
% meta = loadJEB6_ALMVideo(meta,objpth);
meta = loadJEB7_ALMVideo(meta,objpth);
% meta = loadEKH1_ALMVideo(meta,objpth);
% meta = loadEKH3_ALMVideo(meta,objpth);
% meta = loadJGR2_ALMVideo(meta,objpth);
% meta = loadJGR3_ALMVideo(meta,objpth);
% meta = loadJEB14_ALMVideo(meta,objpth);
% meta = loadJEB15_ALMVideo(meta,objpth);

meta = meta(1);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);

params.framesPerTrial = numel(obj.time);       % nr. of frames per trial


%% trials

% here, we specify trials we want to use
% this is important for the analysis, but also just for computational
% reasons
% the more trials, the more video data, and we might not have enough memory
% to load all the video in memory and perform svd
% We only load/analyze the time points we need from the video data, and we crop
% the pixels, so it *might* be possible to get away with loading all
% trials, but I haven't tested the limits yet.

% we will use all hit or miss trials that aren't preceded by a no response
% trial - this is because the features we want to use are binary, would
% have to one hot encode to include a no reponse choice feature and the
% design matrix would get even more out of hand

trix = params.trialid{3};
params.trials2use = trix(~ismember(trix, find(obj.bp.no)+1)); % finding (hit|miss) trials where prev trial was not a no response
trials.right = params.trials2use(ismember(params.trials2use,params.trialid{1}));
trials.left = params.trials2use(ismember(params.trials2use,params.trialid{2}));

% sample equal number of left and right trials
params.nTrials2use = min(numel(trials.right),numel(trials.left));
params.trials2use = [randsample(trials.right,params.nTrials2use,false) ; randsample(trials.left,params.nTrials2use,false)];
trials.right = sort(params.trials2use(1:params.nTrials2use));
trials.rightix = 1:params.nTrials2use;
trials.left = sort(params.trials2use((params.nTrials2use+1):end));
trials.leftix = (params.nTrials2use+1):numel(params.trials2use);
trials.N = numel(trials.left) + numel(trials.right);

%%
clearvars -except meta params obj objpth vidpth trials

%% load raw video data if needed (or if reload set to true)

%--- Cam0 ---%
cam = 'Cam0';
vidout.cam0 = getVideoMeta(cam,vidpth,objpth,meta,params);

% crop pixels (helps with computation and also need to remove artifacts from droplets and empty space)
xcrop = 80; % just specifying how much to crop from left end of video, will always keep from xcrop:(xwidth)
ycrop = 1:200;
[grayframes_Cam0, vidSize_Cam0] = getVideoFrames(vidout.cam0,cam,params,obj,xcrop,ycrop); % (x,y,nFrames*trials)

% plotMovie(grayframes_Cam0)

%--- Cam1 ---%
cam = 'Cam1';
vidout.cam1 = getVideoMeta(cam,vidpth,objpth,meta,params);

% crop pixels (helps with computation and also need to remove artifacts from droplets and empty space)
xcrop = 230:370; % crop pixels (helps with computation and also need to remove artifacts from droplets and empty space)
ycrop = 55:200;
[grayframes_Cam1, vidSize_Cam1] = getVideoFrames(vidout.cam1,cam,params,obj,xcrop,ycrop); % (x,y,nFrames*trials)

% plotMovie(grayframes_Cam1)

%% check for video SVD and create if missing (or if reload set to true)

%--- Cam0 ---%
cam = 'Cam0';
sav = 1; % if 1 - save video SVs, 0 - don't save (saves to same directory as Cam0/Cam1)
vidR_Cam0 = performVideoSVD(cam,vidout.cam0,params,grayframes_Cam0,vidSize_Cam0,sav); % (time,trials,singular vectors)
absVidR_Cam0 = performME_SVD(cam,vidout.cam0,params,grayframes_Cam0,vidSize_Cam0,sav);

%--- Cam0 ---%
cam = 'Cam1';
sav = 1; % if 1 - save video SVs, 0 - don't save (saves to same directory as Cam0/Cam1)
vidR_Cam1 = performVideoSVD(cam,vidout.cam1,params,grayframes_Cam1,vidSize_Cam1,sav); % (time,trials,singular vectors)
absVidR_Cam1 = performME_SVD(cam,vidout.cam1,params,grayframes_Cam1,vidSize_Cam1,sav);


%% synchronize neural and video data

clutm = obj.time;

vidtm = params.tmin:params.viddt:params.tmax;

vidR_Cam0 = syncTimes(vidtm,vidR_Cam0,clutm);
vidR_Cam1 = syncTimes(vidtm,vidR_Cam1,clutm);

absVidR_Cam0 = syncTimes(vidtm,absVidR_Cam0,clutm);
absVidR_Cam1 = syncTimes(vidtm,absVidR_Cam1,clutm);

%% plot video SVDs
close all

% figure; imagesc(squeeze(vidR_Cam0(:,:,5))')
% figure; imagesc(squeeze(vidR_Cam1(:,:,5))')

% figure; imagesc(squeeze(absVidR_Cam0(:,1,:))')
% figure; imagesc(squeeze(absVidR_Cam0(:,2,:))')
% figure; imagesc(squeeze(absVidR_Cam1(:,1,:))')
% figure; imagesc(squeeze(absVidR_Cam1(:,2,:))')

% dim = 5;
% figure; imagesc(squeeze(absVidR_Cam0(:,:,dim))')
% figure; imagesc(squeeze(absVidR_Cam1(:,:,dim))')

right = vidR_Cam1(:,trials.rightix,:);
left = vidR_Cam1(:,trials.leftix,:);

temp = mean(vidR_Cam1(:,:,1:200),3);
figure; imagesc(temp'); colorbar;  colormap(linspecer); caxis([0 1]);


%% ephys data

N = permute(obj.trialdat,[1 3 2]); % single trial binned neural data (time,trials,clusters)

% trim trials
N = N(:,params.trials2use,:);

sz.N = size(N);
sz.vidR_Cam0 = size(vidR_Cam0);
sz.vidR_Cam1 = size(vidR_Cam1);
sz.absVidR_Cam0 = size(absVidR_Cam0);
sz.absVidR_Cam1 = size(absVidR_Cam1);
sz.leg = {'time','trials','clu/dims'};

%% %% video data
% instructed video is average post-gocue video 
% uninstructed video is the remaining (variability that is not captured
% with trial-averaged data)

[~,goCueIx] = min(abs(obj.time - 0));
toPlot = 0;

ix = 1; % right trials
[instructed_vidR_Cam0{ix}, uninstructed_vidR_Cam0{ix}]       = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.rightix,:),goCueIx,toPlot);
[instructed_vidR_Cam1{ix}, uninstructed_vidR_Cam1{ix}]       = getInstructedAndUninstructedSVDs(vidR_Cam1(:,trials.rightix,:),goCueIx,toPlot);
[instructed_absVidR_Cam0{ix}, uninstructed_absVidR_Cam0{ix}] = getInstructedAndUninstructedSVDs(absVidR_Cam0(:,trials.rightix,:),goCueIx,toPlot);
[instructed_absVidR_Cam1{ix}, uninstructed_absVidR_Cam1{ix}] = getInstructedAndUninstructedSVDs(absVidR_Cam1(:,trials.rightix,:),goCueIx,toPlot);

ix = 2; % left trials
[instructed_vidR_Cam0{ix}, uninstructed_vidR_Cam0{ix}]       = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.leftix,:),goCueIx,toPlot);
[instructed_vidR_Cam1{ix}, uninstructed_vidR_Cam1{ix}]       = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.leftix,:),goCueIx,toPlot);
[instructed_absVidR_Cam0{ix}, uninstructed_absVidR_Cam0{ix}] = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.leftix,:),goCueIx,toPlot);
[instructed_absVidR_Cam1{ix}, uninstructed_absVidR_Cam1{ix}] = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.leftix,:),goCueIx,toPlot);

%% reshape video and ephys data

instructed_vidR_Cam0 = catTimeTrialsSVDs(instructed_vidR_Cam0); % (time*trials, ndims)
uninstructed_vidR_Cam0 = catTimeTrialsSVDs(uninstructed_vidR_Cam0);
instructed_vidR_Cam1 = catTimeTrialsSVDs(instructed_vidR_Cam1);
uninstructed_vidR_Cam1 = catTimeTrialsSVDs(uninstructed_vidR_Cam1);
instructed_absVidR_Cam0 = catTimeTrialsSVDs(instructed_absVidR_Cam0);
uninstructed_absVidR_Cam0 = catTimeTrialsSVDs(uninstructed_absVidR_Cam0);
instructed_absVidR_Cam1 = catTimeTrialsSVDs(instructed_absVidR_Cam1);
uninstructed_absVidR_Cam1 = catTimeTrialsSVDs(uninstructed_absVidR_Cam1);


N = reshape(N,size(N,1)*size(N,2),size(N,3)); % (time*trials,clu)
N = bsxfun(@minus, N, mean(N)) ./ std(N); % make zero-mean, unit variance

% figure; imagesc((1:size(N,1)).*params.dt,1:size(N,2),mySmooth(N,21,'reflect')'); colormap(linspecer); colorbar;
% xlabel('time in session (s)')
% ylabel('clusters')

%% task regressors

%-- choice
% lick right = 0, left = 1
choice = nan(obj.bp.Ntrials,1);
for trix = 1:obj.bp.Ntrials
    if obj.bp.no(trix)
        continue
    end
    if obj.bp.R(trix) && obj.bp.hit(trix)
        choice(trix) = 0;
    elseif obj.bp.R(trix) && obj.bp.miss(trix)
        choice(trix) = 1;
    elseif obj.bp.L(trix) && obj.bp.hit(trix)
        choice(trix) = 1;
    elseif obj.bp.L(trix) && obj.bp.miss(trix)
        choice(trix) = 0;
    end
end

%-- previous choice (every trial after a left choice trial)
% if licked right = -1, left = 1
prevchoice = [0 ; choice(1:end-1)];

% subset choice and previous choice to just the trials we want
choice = choice(params.trials2use);
prevchoice = prevchoice(params.trials2use);

%-- success
success = obj.bp.hit;

%-- previous success
prevsuccess = [0 ; success(1:end-1)];

% subset trials
success = success(params.trials2use);
prevsuccess = prevsuccess(params.trials2use);

%-- reward (all frames after reward given)
reward = obj.bp.ev.reward - obj.bp.ev.goCue;
reward = reward(params.trials2use);
rewardix = nan(size(reward));
for i = 1:numel(reward)
    if ~isnan(reward(i))
        [~,rewardix(i)] = min(abs(obj.time - reward(i)));
    end
end

%-- now shape them accordingly

% choice
temp = zeros(numel(choice),numel(obj.time));
temp(:,1) = choice;
choice = temp'; % (time,trials)
choice = choice(:);
taskEvents(:,1) = choice;

% previous choice
temp = zeros(numel(prevchoice),numel(obj.time));
temp(:,1) = prevchoice;
prevchoice = temp'; % (time,trials)
prevchoice = prevchoice(:);
taskEvents(:,2) = prevchoice;
% figure; plot((1:39560).*params.dt,prevchoice)

% success
temp = zeros(numel(success),numel(obj.time));
temp(:,1) = success;
success = temp'; % (time,trials)
success = success(:);
taskEvents(:,3) = success;

% previous success
temp = zeros(numel(prevsuccess),numel(obj.time));
temp(:,1) = prevsuccess;
prevsuccess = temp'; % (time,trials)
prevsuccess = prevsuccess(:);
taskEvents(:,4) = prevsuccess;

% reward
temp = zeros(numel(reward),numel(obj.time));
for i = 1:numel(rewardix)
    if ~isnan(rewardix(i))
        temp(i,rewardix(i)) = 1; % 1 where reward was given on trial i
    end
end
reward = temp';
reward = reward(:);
taskEvents(:,5) = reward;

% figure; imagesc((1:size(taskEvents,1)).*params.dt,1:size(taskEvents,2),taskEvents'); colormap(gray); colorbar;

% MAKE DESIGN MATRIX
taskLabels = {'choice' 'prevchoice' 'success' 'prevsuccess' 'reward'}; %some task variables
taskEventType = [1 1 1 1 1]; %different type of events.

[taskR, taskIdx] = makeDesignMatrix(taskEvents, taskEventType, params); %make design matrix for task variables

% choiceR = taskR(:,taskIdx==1);
% rewardR = taskR(:,taskIdx==5); figure; plot(rewardR(1:460*10,:)) % plot 10 trials

%% movement regressors

% just using licks. licks before go cue are uninstructed licks
% licks after go cue are instructed licks

licks.L = obj.bp.ev.lickL(params.trials2use);
licks.R = obj.bp.ev.lickR(params.trials2use);
lick = zeros(numel(params.trials2use),numel(obj.time));
for trix = 1:numel(params.trials2use)
    L = licks.L{trix} - obj.bp.ev.goCue(params.trials2use(trix));
    R = licks.R{trix} - obj.bp.ev.goCue(params.trials2use(trix));
    
    L = L(L>=params.tmin & L<=params.tmax);
    R = R(R>=params.tmin & R<=params.tmax);

    licktm = sort([L R])';
    for i = 1:numel(licktm)
        [~,ix] = min(abs(obj.time - licktm(i)));
        lick(trix,ix) = 1; % 1 at time point where a lick occurred for each trial
    end
end
lick = lick'; % (time, trials) binary matrix, 1 if lick at that time, 0 o/w

moveEvents(:,1) = lick(:);


% MAKE DESIGN MATRIX

moveLabels = {'licks'}; % some movement variables
moveEventType = [3]; %different type of events.

[moveR, moveIdx] = makeDesignMatrix(moveEvents, moveEventType, params); %make design matrix for task variables


% figure; plot(moveR(:,moveIdx==1),'r')
% hold on
% plot(moveR(:,moveIdx==2),'g')

%% make full design matrix

clear fullR moveLabels_ regLabels regIdx

% not using instructed_vidR since it's rank-deficient

fullR = [taskR, moveR, ...
         instructed_absVidR_Cam0, uninstructed_absVidR_Cam0, ...
         instructed_absVidR_Cam1, uninstructed_absVidR_Cam1, ...
         instructed_vidR_Cam0, uninstructed_vidR_Cam0, ...
         instructed_vidR_Cam1, uninstructed_vidR_Cam1]; %make new, single design matrix
% fullR = [taskR, moveR, instructed_absVidR, uninstructed_absVidR, uninstructed_vidR]; %make new, single design matrix


fullR = bsxfun(@minus, fullR, mean(fullR, 1));

moveLabels_ = [moveLabels, {'iME_Cam0'}, {'uiME_Cam0'}, {'iME_Cam1'}, {'uiME_Cam1'},  {'iVideo_Cam0'}, {'uiVideo_Cam0'},  {'iVideo_Cam1'}, {'uiVideo_Cam1'}];
% moveLabels = [moveLabels, {'iME'}, {'uiME'}, {'uiVideo'}];


% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
regLabels = cat(2,taskLabels,moveLabels_);

% % index to reconstruct different response kernels
nrDims = params.nrDims;
regIdx = [
            ones(sum(taskIdx==1),1) * find(ismember(regLabels,'choice')); ...
            ones(sum(taskIdx==2),1) * find(ismember(regLabels,'prevchoice')); ...
            ones(sum(taskIdx==3),1) * find(ismember(regLabels,'success')); ...
            ones(sum(taskIdx==4),1) * find(ismember(regLabels,'prevsuccess')); ...
            ones(sum(taskIdx==5),1) * find(ismember(regLabels,'reward')); ...
            ones(sum(moveIdx==1),1) * find(ismember(regLabels,'licks')); ...
            ones(nrDims,1) * find(ismember(regLabels,'iME_Cam0')); ...
            ones(nrDims,1) * find(ismember(regLabels,'uiME_Cam0')); ...
            ones(nrDims,1) * find(ismember(regLabels,'iME_Cam1')); ...
            ones(nrDims,1) * find(ismember(regLabels,'uiME_Cam1')); ...
            ones(nrDims,1) * find(ismember(regLabels,'iVideo_Cam0')); ...
            ones(nrDims,1) * find(ismember(regLabels,'uiVideo_Cam0')); ...
            ones(nrDims,1) * find(ismember(regLabels,'iVideo_Cam1')); ...
            ones(nrDims,1) * find(ismember(regLabels,'uiVideo_Cam1')); ...
         ];

% regIdx = [
%             ones(sum(taskIdx==1),1) * find(ismember(regLabels,'choice')); ...
%             ones(sum(taskIdx==2),1) * find(ismember(regLabels,'prevchoice')); ...
%             ones(sum(taskIdx==3),1) * find(ismember(regLabels,'success')); ...
%             ones(sum(taskIdx==4),1) * find(ismember(regLabels,'prevsuccess')); ...
%             ones(sum(taskIdx==5),1) * find(ismember(regLabels,'reward')); ...
%             ones(sum(moveIdx==1),1) * find(ismember(regLabels,'licks')); ...
%             ones(nrDims*2,1) * find(ismember(regLabels,'iME')); ...
%             ones(nrDims,1) * find(ismember(regLabels,'uiME')); ...
%             ones(nrDims*2,1) * find(ismember(regLabels,'iVideo')); ...
%             ones(nrDims,1) * find(ismember(regLabels,'uiVideo')) ...
%          ];

%% run QR and check for rank-defficiency. This will show whether a given regressor is highly collinear with other regressors in the design matrix.
% The resulting plot ranges from 0 to 1 for each regressor, with 1 being
% fully orthogonal to all preceeding regressors in the matrix and 0 being
% fully redundant. Having fully redundant regressors in the matrix will
% break the model, so in this example those regressors are removed. In
% practice, you should understand where the redundancy is coming from and
% change your model design to avoid it in the first place!

rejIdx = false(1,size(fullR,2));
[~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize normalized design matrix
figure; plot(abs(diag(fullQRR)),'linewidth',2); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
axis square; ylabel('Norm. vector angle'); xlabel('Regressors');
if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
    temp = ~(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1)));
    fprintf('Design matrix is rank-defficient. Removing %d/%d additional regressors.\n', sum(temp), sum(~rejIdx));
    rejIdx(~rejIdx) = temp; %reject regressors that cause rank-defficint matrix
end

%% run full model fit for single trials

[~, dimBeta] = ridgeMML(N, fullR, true); %make model fit
fullFit = (fullR - mean(fullR, 1)) * dimBeta;

disp('Full model fit completed');

uninstructedLabels = {'uiME_Cam0','uiVideo_Cam0','uiME_Cam1','uiVideo_Cam1'};
lbIdx = find(ismember(regLabels,uninstructedLabels));
uiIdx = ismember(regIdx ,lbIdx) ; % logical of uninstructed vars in cols of fullR

instructedLabels = {'licks','iME_Cam0','iVideo_Cam0','iME_Cam1','iVideo_Cam1'};
lbIdx = find(ismember(regLabels,instructedLabels));
iIdx = ismember(regIdx ,lbIdx) ; % logical of instructed vars in cols of fullR

lbIdx = find(ismember(regLabels,taskLabels));
taskIdx = ismember(regIdx ,lbIdx) ; % logical of task vars in cols of fullR

ui_fit = (fullR(:, uiIdx) - mean(fullR(:,uiIdx),1)) * dimBeta(uiIdx, :); % uninstructed vars fit 
i_fit = (fullR(:, iIdx) - mean(fullR(:,iIdx),1)) * dimBeta(iIdx,:);      % instructed vars fit 
task_fit = (fullR(:, taskIdx) - mean(fullR(:,taskIdx),1)) * dimBeta(taskIdx,:);   % instructed vars fit 

%% variance explained (https://economictheoryblog.com/2014/11/05/proof/)

R = corrcoef(N,fullFit);
ve_full = R(1,2)^2;

R = corrcoef(N,i_fit);
ve_i = R(1,2)^2;

R = corrcoef(N,ui_fit);
ve_ui = R(1,2)^2;

R = corrcoef(N,task_fit);
ve_task = R(1,2)^2;

ve = [ve_full,ve_ui,ve_i,ve_task];

%
figure;
X = categorical({'Full','Unin. Move','In. Move','Task'});
h = histogram('Categories', X, 'BinCounts', ve, 'EdgeColor','none');
ylabel('Frac. VE')
ax = gca;
ax.FontSize = 13;

%% run cross-validation
% %full model - this will take a moment
[Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels, fullR2] = crossValModel(fullR, N, regLabels, regIdx, regLabels, params.folds);
save(fullfile(vidout.cam0.svpth, 'cvFull.mat'), 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels', 'fullR2'); %save some results

% %task model alone - this will take a moment
[Vtask, taskBeta, taskR, taskIdx, taskRidge, taskLabels, taskR2] = crossValModel(fullR, N, taskLabels, regIdx, regLabels, params.folds);
save(fullfile(vidout.cam0.svpth, 'cvTask.mat'), 'Vtask', 'taskBeta', 'taskR', 'taskIdx', 'taskRidge', 'taskLabels', 'taskR2'); %save some results

% %movement model alone - this will take a moment
[Vmove, moveBeta, moveR, moveIdx, moveRidge, moveLabels, moveR2] = crossValModel(fullR, N, moveLabels_, regIdx, regLabels, params.folds);
save(fullfile(vidout.cam0.svpth, 'cvMove.mat'), 'Vmove', 'moveBeta', 'moveR', 'moveIdx', 'moveRidge', 'moveLabels', 'moveR2'); %save some results

%% plot R2

r2 = [mean(fullR2),mean(taskR2),mean(moveR2)];

%
figure;
X = categorical({'Full','Task','Move'});
h = histogram('Categories', X, 'BinCounts', r2, 'EdgeColor','none');
ylabel('R2')
ax = gca;
ax.FontSize = 13;


%% reconstruct PSTHs 

% N is original data in (time*trials,neurons)
% fullFit is reocnstructed data of same size

N_ = reshape(N,sz.N(1),sz.N(2),sz.N(3));
fullFit = reshape(fullFit,sz.N(1),sz.N(2),sz.N(3));

% pick k random neurons
k = 20;

%
% close all
% clus = randsample(size(N_,3),k,false);
clus = 1:size(N_,3);
f = figure;
f.Position = [203         332        1418         635];
t = tiledlayout('flow');
alph = 0.1;
for i = 1:numel(clus)
    clu = clus(i);
    
    tempN = N_(:,:,clu);
%     tempNhat = fullFit(:,:,clu);
    tempNhat = ui_fit(:,:,clu);

    ax = nexttile; hold on;
    plot(obj.time,mean(tempN,2),'Color','k','LineWidth',3);
    plot(obj.time,mean(tempNhat,2),'Color',[108, 212, 108]./255,'LineWidth',1);

%     shadedErrorBar(obj.time,mean(tempN,2),std(tempN,[],2)./sqrt(size(tempN,2)), {'Color','k','LineWidth',6},0,ax)
%     shadedErrorBar(obj.time,mean(tempNhat,2),std(tempNhat,[],2)./sqrt(size(tempNhat,2)), {'Color',[41, 128, 64]./255,'LineWidth',1},0,ax)

%     patchline(obj.time,psthN,'EdgeColor',[0 0 0],'LineWidth',5,'EdgeAlpha',0.3);
%     patchline(obj.time,psthFullFit,'EdgeColor',[41, 128, 64]./255,'LineWidth',2,'EdgeAlpha',1)
    xlim([obj.time(10) obj.time(end)])
    ax.FontSize = 13;
    xline(0,'k--')
    ax.XTick = [];
    ax.YTick = [];
end
title(t, 'PSTHs and reconstruction')
xlabel(t, 'Time (s) from go cue')
ylabel(t, 'Centered Firing Rate (Hz)')








%% Helper functions

function vidout = getVideoMeta(cam,vidpth,datapth,meta,params)
vidpth = fullfile(vidpth,meta.anm,meta.date,cam); % assumes that data is stored in datapth/Video/<anm>/<date>/Cam0/Cam1
vids = dir(vidpth);
vids = {vids.name}';
vids = patternMatchCellArray(vids,{meta.anm,meta.date},'all');
vids = natsortfiles(vids);

vidout.vidpth = vidpth;
vidout.vids = vids(params.trials2use); % trim trials
vidout.vidtm = params.tmin:params.viddt:params.tmax;
vidout.nFrames = numel(vidout.vidtm);
vidout.svpth = fullfile(datapth,'Video',meta.anm,meta.date);
end




function [grayframes_Cam, vidSize_Cam] = getVideoFrames(vidout,cam,params,obj,xcrop,ycrop)

reload = params.(['reload_' cam]); % if true, loads video data even if SVs already exist

fpth = vidout.svpth;
vids = vidout.vids;
nFrames = vidout.nFrames;
if ~exist(fullfile(fpth ,['vidR_' cam   '.mat']), 'file') || ~exist(fullfile(fpth ,['absVidR_' cam   '.mat']), 'file') || reload

    for vidix = 1:numel(vids) % vidix is same as trial since vids are sorted by timestamp
        disp(['Loading vid for trial: ' num2str(vidix) '/' num2str(numel(vids))])
        vid = vids{vidix};

        % video time
        tm = obj.traj{1}(params.trials2use(vidix)).frameTimes - 0.5 - obj.bp.ev.(params.alignEvent)(params.trials2use(vidix));
        % find idx of tmin and tmax in video time
        ix1 = find(tm >= params.tmin, 1, 'first');
        ix2 = find(tm <= params.tmax, 1, 'last');
        nFramesTrial = ix2 - ix1;
        %         t1 = tm(ix1);
        %         t2 = tm(ix2);
        if nFramesTrial ~= nFrames
            ix2 = ix2 + (nFrames-nFramesTrial);
        end

        Cnt = 0;
        v = VideoReader(fullfile(vidout.vidpth, vid));
        %     cVideo = zeros(v.Width,v.Height,nFrames, 'uint8');

        frames = read(v,[ix1 ix2]);

        if numel(xcrop) == 1
            xcrop = xcrop:size(frames,2); % crop pixels (helps with computation and also need to remove artifacts from droplets and empty space)
        end

        for fix = 1:nFrames
            mov = im2gray(frames(ycrop,xcrop,:,fix));
            grayframes(:,:,fix,vidix) = mov; % (x,y,frames,trials)
        end

    end
    vidSize_Cam = size(grayframes);
    grayframes_Cam = reshape(grayframes,vidSize_Cam(1),vidSize_Cam(2),vidSize_Cam(3)*vidSize_Cam(4)); % (x,y,frames*trials)
else
    disp('video and motion energy SVDs already exist and params.reload set to not reload video data - skipping')
    grayframes_Cam = nan;
    vidSize_Cam = nan;
    return
end

end


function vidR = performVideoSVD(cam,vidout,params,grayframes,vidSize,sav)

reload = params.(['reload_' cam]); % if true, computes video SVs even if they already exist

fpth = vidout.svpth;

nrDims = params.nrDims; % svd dims to keep
if ~exist(fullfile(fpth, ['vidR_' cam   '.mat']), 'file') || reload
    cVideo = reshape(grayframes, [], size(grayframes,3)); % merge all pixels (x*y,nFrames*nTrials)
    %     cVideo = (cVideo - uint8(mean(cVideo,2))); % mean center

    [~,s,Vr] = fsvd(single(cVideo),nrDims,1,0,0);
    vidR = s*Vr';
    vidR = bsxfun(@minus, vidR, mean(vidR, 2));
    vidR = bsxfun(@rdivide, vidR, std(vidR, [], 2))'; %this is the video regressor in the model
    vidR = reshape(vidR,vidSize(3),vidSize(4),nrDims); % (time,trials,videodims)
    if sav
        savedir = fpth;
        if ~exist(savedir,'dir')
            mkdir(savedir)
        end
        save(fullfile(savedir, ['vidR_' cam   '.mat']),'vidR');
    end
else
    temp = load(fullfile(fpth ,['vidR_' cam   '.mat']),'vidR');
    vidR = temp.vidR; clear temp;
end

end


function absVidR = performME_SVD(cam,vidout,params,grayframes,vidSize,sav)

reload = params.(['reload_' cam]); % if true, computes video SVs even if they already exist

fpth = vidout.svpth;
nrDims = params.nrDims;
if ~exist(fullfile(fpth ,['absVidR_' cam   '.mat']), 'file') || reload
    cVideo = grayframes; %merge all pixels
    cVideo = cat(3,cVideo(:,:,1),cVideo); %duplicate first frame
    cVideo = abs(diff(cVideo,[],3)); %compute absolute of temporal derivative (motion energy)
    cVideo = reshape(cVideo, [], size(cVideo,3)); %merge all pixels
    %     cVideo = (cVideo - uint8(mean(cVideo,2))); % mean center

    [~, s, Vr] = fsvd(single(cVideo), nrDims, 1, 0, 0);
    absVidR = s * Vr';
    absVidR = bsxfun(@minus, absVidR, mean(absVidR, 2));
    absVidR = bsxfun(@rdivide, absVidR, std(absVidR, [], 2))'; %this is the motion energy regressor in the model
    absVidR = reshape(absVidR,vidSize(3),vidSize(4),nrDims); % (time,trials,videodims)
    if sav
        savedir = fpth;
        if ~exist(savedir,'dir')
            mkdir(savedir)
        end
        save(fullfile(savedir, ['absVidR_' cam   '.mat']),'absVidR');
    end
else
    temp = load(fullfile(fpth ,['absVidR_' cam   '.mat']),'absVidR');
    absVidR = temp.absVidR; clear temp;
end

end



























