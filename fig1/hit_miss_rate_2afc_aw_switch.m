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
params.condition(1)     = {'R&hit&~stim.enable&~autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&hit&~stim.enable&autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&autowater&((1:Ntrials)>20)&((1:Ntrials)<(Ntrials-20))'};             % left hits, no stim, aw off

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

meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB15_ALMVideo(meta,datapth);


params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params);


%%

% for each session
% find autowater blocks
% remove last autowater block if it leads to end of session
% find shortest aw block

% find shortest aw block across sessions

for sessix = 1:numel(meta)
    awtrials = obj(sessix).bp.autowater;
    startlocs = find(diff(awtrials) == 1) + 1; %location of 1st high value in each run
    endlocs = find(diff(awtrials) == -1);  %location of last high value in each run. Guaranteed to be different from startlocs.
    if (numel(startlocs) ~= numel(endlocs)) && (endlocs(end) < startlocs(end)) % if last aw block leads to end of session
        startlocs(end) = [];
    end

    %     figure; hold on
    %     plot(awtrials)
    %     plot(startlocs,0.5*ones(size(startlocs)),'g.','MarkerSize',20)
    %     plot(endlocs,0.5*ones(size(endlocs)),'r.','MarkerSize',20)


    awlens(sessix) = min(endlocs - startlocs);
end

minawblock = min(awlens);

trials.before = minawblock;
trials.after = minawblock;


%%
                                                             
clear hr mr nr
% for each session
% find autowater blocks
% remove last autowater block if it leads to end of session
% for each aw block, calculate hit rate,miss rate, nr rate trials.before and trials.after 

for sessix = 1:numel(meta)
    awtrials = obj(sessix).bp.autowater;
    startlocs = find(diff(awtrials) == 1) + 1; %location of 1st high value in each run
    endlocs = find(diff(awtrials) == -1);  %location of last high value in each run. Guaranteed to be different from startlocs.
    if (numel(startlocs) ~= numel(endlocs)) && (endlocs(end) < startlocs(end)) % if last aw block leads to end of session
        startlocs(end) = [];
    end

    % startlocs is trial num of the start of an aw block
    % for each startlocs, calculate performance up to that point (with
    % trials.before and trials.after buffer)
    for s = 1:numel(startlocs)
        awstarttrix = startlocs(s);

        trixs = awstarttrix-trials.before:awstarttrix+trials.after;
        for i = 1:numel(trixs)
            hr{sessix}(s,i) = sum(obj(sessix).bp.hit(1:trixs(i))) ./ trixs(i); % hit rate (aw block,trial leading/after)
            mr{sessix}(s,i) = sum(obj(sessix).bp.miss(1:trixs(i))) ./ trixs(i); % miss rate(aw block,trial leading/after)
            nr{sessix}(s,i) = sum(obj(sessix).bp.no(1:trixs(i))) ./ trixs(i); % no resp rate (aw block,trial leading/after)
        end
    end

end

% combine across sessions
c = hr;
rows = size(c,1);
cols = size(c,2);
m = cell(rows,1);
% Concatenate one dim first
for n=1:rows
    m{n} = cat(1,c{n,:});
end
% Now concatenate the single column of cells into a matrix
hr = cat(1,m{:}); % (aw blocks all sessions, trials)

c = mr;
rows = size(c,1);
cols = size(c,2);
m = cell(rows,1);
% Concatenate one dim first
for n=1:rows
    m{n} = cat(1,c{n,:});
end
% Now concatenate the single column of cells into a matrix
mr = cat(1,m{:}); % (aw blocks all sessions, trials)

c = nr;
rows = size(c,1);
cols = size(c,2);
m = cell(rows,1);
% Concatenate one dim first
for n=1:rows
    m{n} = cat(1,c{n,:});
end
% Now concatenate the single column of cells into a matrix
nr = cat(1,m{:}); % (aw blocks all sessions, trials)


%% plot

figure; 
ax = gca;
hold on
shadedErrorBar(-trials.before:trials.after,mean(hr,1),std(hr,[],1)./sqrt(numel(meta)), {'Color',[45, 156, 71]./255, 'LineWidth',2},0.5,ax)
shadedErrorBar(-trials.before:trials.after,mean(mr,1),std(mr,[],1)./sqrt(numel(meta)), {'Color',[0 0 0], 'LineWidth',2},0.5,ax)
shadedErrorBar(-trials.before:trials.after,mean(nr,1),std(nr,[],1)./sqrt(numel(meta)), {'Color',[150 150 150]./255, 'LineWidth',2},0.5,ax)

xline(0,'k:','LineWidth',2)
xlabel('Trials to switch')
ylabel('Fraction of trials')
ax.FontSize = 13;


h(1) = plot(nan,nan,'Color',[45, 156, 71]./255,'LineWidth',2);
h(2) = plot(nan,nan,'Color',[0 0 0],'LineWidth',2);
h(3) = plot(nan,nan,'Color',[150 150 150]./255,'LineWidth',2);

legend(h,{'hit rate','miss rate','no resp rate'},'Location','best')




