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


params.tmin = -2.35;
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

params.behav_only = 1;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM --- 
% meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth);
% meta = loadEKH3_ALMVideo(meta,datapth);
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);
% meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

meta = meta(1); % jeb15 07-27

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth, 1);
end

%%

% set session data to use
anm = meta.anm;
date = meta.date;

dat.obj = obj;
dat.params = params;
dat.meta = meta;

%%

% trix = 1:obj.bp.Ntrials;



% cols = linspecer(4,'qualitative');



feats = 'motion_energy'; % nose
view = 1; % side cam
% cols = cols(2,:);
cols = getColors;
trix = [24 26 34 44 46 49 67 130 293 324 327 ];
trajoff = 120;
ms = 27;



vidpth = ['C:\Users\munib\Documents\Economo-Lab\data\Video\JEB7\2021-04-29\Cam' num2str(view-1)];
vids = dir(vidpth);
vids = vids(~[vids(:).isdir]);
vids = {vids(:).name};
vidfns = natsortfiles(vids)';

for itrix = 1:numel(trix)
    
    trial = trix(itrix);


    clear v hFig hAx hIm
    v = VideoReader(fullfile(vidpth,vidfns{trial}));


    % get obj trajectories

    frametimes = dat.obj.traj{view}(trial).frameTimes - 0.5;
    sample = dat.obj.bp.ev.sample(trial);
    delay = dat.obj.bp.ev.delay(trial);
    gocue = dat.obj.bp.ev.goCue(trial);


    % resmaple obj.time to frametimes, but first need to get frametimes from
    % -2.5 from go cue to 2.5 to go cue
    tmin = dat.params.tmin;
    tmax = dat.params.tmax;

    [~,ix1] = min(abs(frametimes - gocue - tmin));
    [~,ix2] = min(abs(frametimes - gocue - tmax));

    frames = read(v,[ix1 ix2]);
    frametimes = frametimes(ix1:ix2);


    % trajectories
    sm = 11;
    traj = mySmooth(me.data{trial}(ix1:ix2),sm,'reflect');
    traj = normalize(traj,'range',[10 80]);

    close all

    nFrames = size(frames,4);
    s(nFrames) = struct('cdata',[],'colormap',[]);


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

    hIm = image(uint8(zeros(v.Height,v.Width,3)),...
        'Parent',hAx);


    % loop through frames
    xx = linspace(1,size(frames,2),size(frames,4));
    for iframe = 1:size(frames,4)
        im = frames(:,:,:,iframe);
        hIm.CData = flipud(im);

        % draw traj
        toplot = traj(1:iframe);
        null = toplot;
        null(null>me.moveThresh) = nan;        
        % if traj(iframe)<=me.moveThresh
        %     col = cols.null;
        % else
        %     col = cols.potent;
        % end
        plot(xx(1:iframe),toplot,'Color',cols.potent,'linewidth',1.5);
        plot(xx(1:iframe),null,'Color',cols.null,'linewidth',1.5);
        
        % epoch annotation
        t = '';
        if frametimes(iframe) >= sample && frametimes(iframe) < delay
            t = 'sample';
        elseif frametimes(iframe) >= delay && frametimes(iframe) < gocue
            delete(tt);
            t = 'delay';
        elseif frametimes(iframe) >= gocue && frametimes(iframe) <= (gocue + 0.35)
            t = 'go cue';
        end
        tt = text(hAx,15,220,t,'Color','w','FontSize',12);

        drawnow
        % Save the frame in structure for later saving to video file
        s(iframe) = getframe(hAx);
%         b = getframe(hTraceAx);
%         s(iframe).cdata = cat(1,a.cdata,b.cdata);
%         s(iframe).colormap = [];
%         s(iframe).cdata(236:238,:,:) = 0;

        delete(tt);

    end



    %
    % Remove any unused structure array elements

    % Write to a video file
    % This could be done within the original loop, but I wanted to show it
    % separately
    fnout = [anm '_' date '_trial' num2str(trial) '.mp4'];
    vOut = VideoWriter(fullfile('meVideos',fnout),'MPEG-4');
    vOut.FrameRate = 100;
    vOut.Quality = 100;
    open(vOut)
    for k = 1:numel(s)
        writeVideo(vOut,s(k))
    end
    close(vOut)
    
    cla(hAx);

end























