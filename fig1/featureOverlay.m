clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig3/')))
rmpath(genpath(fullfile(utilspth,'fig2/')))
rmpath(genpath(fullfile(utilspth,'figx/')))
rmpath(genpath(fullfile(utilspth,'figNP/')))
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
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};      % right hits, 2afc (2)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};      % left hit, 2afc   (3)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};     % right miss, 2afc (4)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};     % left miss, 2afc  (5)
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};        % 2afc hits        (6)
params.condition(end+1) = {'hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};         % aw hits          (7)
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % right hits, aw   (8)
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % left hits, aw    (9)


params.tmin = -1.4;
params.tmax = 0.5;
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

params.behav_only = 1;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM ---
% meta = loadJEB6_ALMVideo(meta,datapth);
% meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth);
% meta = loadEKH3_ALMVideo(meta,datapth);
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

meta = meta(2); % jeb15 07-27

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

%%

% set session data to use
anm = meta.anm;
date = meta.date;

dat.obj = obj;
dat.params = params;
dat.meta = meta;

vidpth = ['C:\Users\munib\Documents\Economo-Lab\data\Video\JEB15\2022-07-27\Cam' num2str(view-1)];
vids = dir(vidpth);
vids = vids(~[vids(:).isdir]);
vids = {vids(:).name};
vidfns = natsortfiles(vids)';

%% side cam

close all

view = 1;
vidpth = ['C:\Users\munib\Documents\Economo-Lab\data\Video\JEB15\2022-07-27\Cam' num2str(view-1)];
vids = dir(vidpth);
vids = vids(~[vids(:).isdir]);
vids = {vids(:).name};
vidfns = natsortfiles(vids)';


feats = [4 6]; % jaw,nose
trix = 1:obj.bp.Ntrials;
cols = {[1 0 0], [0 1 0]};
ds = 20;

time2use = [-0.9 0];
ix = findTimeIX(obj.time,time2use);
ix = ix(1):ix(2);

traj = [];
for itrix = 1:numel(trix)
    trial = trix(itrix);
    % get obj trajectories
    frametimes = dat.obj.traj{view}(trial).frameTimes - 0.5;
    delay = dat.obj.bp.ev.delay(trial);
    gocue = dat.obj.bp.ev.goCue(trial);
    % find frametimes b/w delay and go cue
    mask = frametimes>=delay & frametimes <=gocue;
    traj = cat(1,traj,obj.traj{view}(trial).ts(mask,[1 2],feats));
end
traj = downsample(traj,ds);

% load single frame
trial = 25;
v = VideoReader(fullfile(vidpth,vidfns{trial}));
im = read(v,100);

% Set up figure, axes to hold image, image, and plot
% This allows us to keep everything positioned in the proper order, and
% just change the image data and the plot location every loop iteration
hFig = figure('MenuBar','none',...
    'Units','pixels','Color',[0 0 0],...
    'Position',[658 404 v.Width*1.15 v.Height*1.25]);

hAx = axes('Parent',hFig,...
    'Units','pixels',...
    'Position',[v.Width*0.05 v.Height*0.05 v.Width v.Height ],...
    'NextPlot','add',...
    'Visible','off',...
    'XTick',[],...
    'YTick',[]);
%     hAx = axes('Parent',hFig,...
%         'Units','pixels',...
%         'Position',[v.Width v.Height v.Width v.Height ],...
%         'NextPlot','add',...
%         'Visible','off',...
%         'XTick',[],...
%         'YTick',[]);
hIm = image(uint8(zeros(v.Height,v.Width,3)),...
    'Parent',hAx);


hIm.CData = flipud(im);
temp = traj(:,:,1);
temp = mySmooth(temp,11,'reflect');
plot(temp(:,1)+30,size(im,1) - temp(:,2)- 5,'Color',cols{1});
temp = traj(:,:,2);
temp = mySmooth(temp,11,'reflect');
plot(temp(:,1),size(im,1) - temp(:,2)- 5,'Color',cols{2});

saveas(hFig,"figs/sidecam.svg",'svg')



%% bottom cam

close all

view = 2;
vidpth = ['C:\Users\munib\Documents\Economo-Lab\data\Video\JEB15\2022-07-27\Cam' num2str(view-1)];
vids = dir(vidpth);
vids = vids(~[vids(:).isdir]);
vids = {vids(:).name};
vidfns = natsortfiles(vids)';

feats = [5 6 8 9 10]; % toppaw,botpaw,jaw,topnose,botnose
trix = 1:obj.bp.Ntrials;
cols = {[0 0 1], [0 0 1], [1 0 0], [0.9216    0.9216    0.2039], [0.9686    0.5490    0.1255]};
ds = 20;

time2use = [-0.9 0];
ix = findTimeIX(obj.time,time2use);
ix = ix(1):ix(2);

traj = [];
for itrix = 1:numel(trix)
    trial = trix(itrix);
    % get obj trajectories
    frametimes = dat.obj.traj{view}(trial).frameTimes - 0.5;
    delay = dat.obj.bp.ev.delay(trial);
    gocue = dat.obj.bp.ev.goCue(trial);
    % find frametimes b/w delay and go cue
    mask = frametimes>=delay & frametimes <=gocue;
    traj = cat(1,traj,obj.traj{view}(trial).ts(mask,[1 2],feats));
end
traj = downsample(traj,ds);

% load single frame
trial = 25;
v = VideoReader(fullfile(vidpth,vidfns{trial}));
im = read(v,100);

% Set up figure, axes to hold image, image, and plot
% This allows us to keep everything positioned in the proper order, and
% just change the image data and the plot location every loop iteration
hFig = figure('MenuBar','none',...
    'Units','pixels','Color',[0 0 0],...
    'Position',[658 404 v.Width*1.15 v.Height*1.25]);

hAx = axes('Parent',hFig,...
    'Units','pixels',...
    'Position',[v.Width*0.05 v.Height*0.05 v.Width v.Height ],...
    'NextPlot','add',...
    'Visible','off',...
    'XTick',[],...
    'YTick',[]);
%     hAx = axes('Parent',hFig,...
%         'Units','pixels',...
%         'Position',[v.Width v.Height v.Width v.Height ],...
%         'NextPlot','add',...
%         'Visible','off',...
%         'XTick',[],...
%         'YTick',[]);
hIm = image(uint8(zeros(v.Height,v.Width,3)),...
    'Parent',hAx);


hIm.CData = flipud(im);
temp = traj(:,:,1);
temp = mySmooth(temp,21,'reflect');
plot(temp(:,1),size(im,1) - temp(:,2)- 5,'Color',cols{1});
temp = traj(:,:,2);
temp = mySmooth(temp,5,'reflect');
plot(temp(:,1),size(im,1) - temp(:,2)- 5,'Color',cols{2});
temp = traj(:,:,3);
temp = mySmooth(temp,5,'reflect');
plot(temp(:,1),size(im,1) - temp(:,2)- 5,'Color',cols{3});
temp = traj(:,:,4);
temp = mySmooth(temp,5,'reflect');
plot(temp(:,1),size(im,1) - temp(:,2)- 5,'Color',cols{4});
temp = traj(:,:,5);
temp = mySmooth(temp,5,'reflect');
plot(temp(:,1),size(im,1) - temp(:,2)- 5,'Color',cols{5});

saveas(hFig,"figs/bottomcam.svg",'svg')
















