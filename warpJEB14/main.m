clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')));
rmpath(genpath(fullfile(utilspth,'fig2/')));
rmpath(genpath(fullfile(utilspth,'figNP/')));
rmpath(genpath(fullfile(utilspth,'mc_stim/')));

% add paths for figure specific functions
addpath(genpath(pwd))

clc

%%

% JEB14 delay period = 0.9
% all other sessions delay period = 0.7

% warping delay period for all jeb14 trials to 0.9 seconds, will only use
% ~early trials, however b/c of this

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

meta = loadJEB14_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

objs = loadObjs(meta);


%% warp

for sessix = 1:numel(objs)

    obj = objs(sessix);

    warp.delayOld = 0.7;
    warp.delayNew = 0.9;
    warp.offset = 0.2; % need to add 0.2 seconds to everything that happens post-delay

    p = zeros(obj.bp.Ntrials,2);

    for trix = 1:obj.bp.Ntrials
        temp = obj.bp.ev.goCue(trix)-obj.bp.ev.delay(trix);
        if ~isNearValue(temp,0.7)
            continue
        end
        x = [obj.bp.ev.delay(trix) obj.bp.ev.goCue(trix)]; % original
        y = [obj.bp.ev.delay(trix) obj.bp.ev.delay(trix)+warp.delayNew]; % warped

        % fit original time to warped time
        p(trix,:) = polyfit(x,y,1);

        %     yhat = polyval(p,x);
    end



    % warp obj.bp.ev fields
    obj = warpEvents(obj,warp);

    % warp spikes
    obj = warpSpikes(obj,warp,p);

    % warp video
    obj = warpVideo(obj,warp,p);

    % save object
    datapth = meta(sessix).datapth;
    pthsplit = strsplit(datapth,'\');
    fn = pthsplit{end};
    pth = strjoin(pthsplit(1:8),'\');

    fnsplit = strsplit(fn,'.');
    fn = [fnsplit{1} '_delayWarped.mat'];

    svpth = fullfile(pth,fn);
    save(svpth,'obj','-v7.3')

end

%% save object

for i = 1:numel(obj)
    
end


%% helper functions

function obj = warpEvents(obj,warp)
for trix = 1:obj.bp.Ntrials
    temp = obj.bp.ev.goCue(trix)-obj.bp.ev.delay(trix);
    if ~isNearValue(temp,0.7)
        continue
    end
    obj.bp.ev.goCue(trix) = obj.bp.ev.goCue(trix) + warp.offset;
    obj.bp.ev.reward(trix) = obj.bp.ev.reward(trix) + warp.offset;
    lickL = obj.bp.ev.lickL{trix};
    lickL(lickL>=(obj.bp.ev.goCue(trix) - warp.offset)) = lickL(lickL>=(obj.bp.ev.goCue(trix) - warp.offset)) + warp.offset;
    obj.bp.ev.lickL{trix} = lickL;
    lickR = obj.bp.ev.lickR{trix};
    lickR(lickR>=(obj.bp.ev.goCue(trix) - warp.offset)) = lickR(lickR>=(obj.bp.ev.goCue(trix) - warp.offset)) + warp.offset;
    obj.bp.ev.lickR{trix} = lickR;
end
end


function obj = warpSpikes(obj,warp,pfit)


for prbix = 1:numel(obj.clu)
    for cluix = 1:numel(obj.clu{prbix})
        for trix = 1:obj.bp.Ntrials
            temp = obj.bp.ev.goCue(trix)-obj.bp.ev.delay(trix);
            if ~isNearValue(temp,0.9)
                continue
            end

            % find spike times for current trial
            spkmask = ismember(obj.clu{prbix}(cluix).trial,trix);
            spktm = obj.clu{prbix}(cluix).trialtm(spkmask);
            warptm = spktm;

            if isempty(spktm)
                continue
            end

            % warp spikes between delay and gocue-0.2
            % add 0.2 to spikes >= gocue
            delay = obj.bp.ev.delay(trix);
            gc = obj.bp.ev.goCue(trix) - warp.offset;

            mask = (spktm>=delay) & (spktm<=gc );
            tm = spktm(mask);
            if ~isempty(tm)
                % warp spikes
                warptm(mask) = polyval(pfit(trix,:),tm);
            end

            mask = (spktm>=gc);
            tm = spktm(mask);
            if ~isempty(tm)
                % add 0.2 to spks
                warptm(mask) = tm + 0.2;
            end

            obj.clu{prbix}(cluix).trialtm(spkmask) = warptm; % after warping

        end
    end
end
end % warpSpikes



function obj = warpVideo(obj,warp,pfit)

for view = 1:numel(obj.traj)

    for trix = 1:obj.bp.Ntrials
        temp = obj.bp.ev.goCue(trix)-obj.bp.ev.delay(trix);
        if ~isNearValue(temp,0.9)
            continue
        end

        % warp frameTimes between delay and gocue-0.2
        % add 0.2 to frameTimes >= gocue
        delay = obj.bp.ev.delay(trix);
        gc = obj.bp.ev.goCue(trix) - warp.offset;

        ft = obj.traj{view}(trix).frameTimes - 0.5;
        warptm = ft;

        mask = (ft>=delay) & (ft<=gc );
        tm = ft(mask);
        if ~isempty(tm)
            % warp frametimes
            warptm(mask) = polyval(pfit(trix,:),tm);
        end

        mask = (ft>=gc);
        tm = ft(mask);
        if ~isempty(tm)
            % add 0.2 to spks
            warptm(mask) = tm + 0.2;
        end


        obj.traj{view}(trix).frameTimes = warptm + 0.5;



    end
end
end % warpVideo
















