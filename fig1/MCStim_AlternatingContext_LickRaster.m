% Quantifying behavioral performance and cortical dependence in the
% Alternating Context Task
clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig3/')))
rmpath(genpath(fullfile(utilspth,'mc_stim/')))
rmpath(genpath(fullfile(utilspth,'MotionMapper/')))


% add paths for figure specific functions
addpath(genpath(pwd))

clc

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% set conditions to calculate behavioral performance for 
% right
params.condition(1) =     {'R&~stim.enable&~autowater'};             % right DR
params.condition(end+1) = {'R&stim.enable&~autowater'};              % right DR stim

params.condition(end+1) = {'R&~stim.enable&autowater'};              % right WC
params.condition(end+1) = {'R&stim.enable&autowater'};               % right WC stim


% left
params.condition(end+1) = {'L&~stim.enable&~autowater'};             % left DR
params.condition(end+1) = {'L&stim.enable&~autowater'};              % left DR stim

params.condition(end+1) = {'L&~stim.enable&autowater'};              % left WC
params.condition(end+1) = {'L&stim.enable&autowater'};               % left WC stim



params.alignEvent = 'goCue';
params.tmin = -2.4; 
params.tmax = 2.5;
params.dt = 1/100;

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw'}};
params.smooth = 15;
params.advance_movement = 0;
params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this mu

params.behav_only = 1;
%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';


meta = [];
meta = loadMAH13_MCGCStim(meta,datapth);
meta = loadMAH14_MCGCStim(meta,datapth);

meta = meta(2);
%%
% ----------------------------------------------
% -- Behavioral and Video Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadBehavSessionData(meta,params);
for sessix = 1:length(obj)
    obj(sessix).time = params(sessix).tmin:params(sessix).dt:params(sessix).tmax;
end
% %% Get kinematics
% for sessix = 1:length(meta)
%     kin(sessix) = getKinematics_NoME(obj(sessix),params(sessix));
% end

%% conds from top - right DR ctrl, right DR stim, right WC ctrl, right WC stim

close all

sessix = 1; 

cols = getColors();
clrs.rhit = cols.rhit;
clrs.lhit = cols.lhit;

conds = [1 2 3 4];

plotLickRasterv2(sessix,clrs,obj,params,'2AFC',conds);

%% conds from top - left DR ctrl, left DR stim, left WC ctrl, left WC stim

% close all

sessix = 1; 

cols = getColors();
clrs.rhit = cols.rhit;
clrs.lhit = cols.lhit;

conds = [5 6 7 8];

plotLickRasterv2(sessix,clrs,obj,params,'2AFC',conds);
%%

% % % %% 2AFC trials: top subplot = left trials (control on top, stim on bottom); bottom subplot = right trials (control on top, stim on bottom)
% % % close all
% % % 
% % % sessix = 1; 
% % % 
% % % cols = getColors();
% % % clrs.rhit = cols.rhit;
% % % clrs.lhit = cols.lhit;
% % % 
% % % conds = [1 2 3 4];
% % % 
% % % figure(1)
% % % plotLickRaster(sessix,clrs,obj,params,'2AFC',conds);
% % % %% AW trials: top subplot = left trials (control on top, stim on bottom); bottom subplot = right trials (control on top, stim on bottom)
% % % clrs.rhit = cols.rhit_aw;
% % % clrs.lhit = cols.lhit_aw;
% % % 
% % % conds = [5 6 7 8];
% % % 
% % % figure(2)
% % % plotLickRaster(sessix,clrs,obj,params,'AW',conds);
% % % %% For supplement: show what tongue looks like during the stim period
% % % sessix = 2;
% % % 
% % % trix2plot = 13;
% % % smooth = 31;
% % % offset = 5;
% % % 
% % % stim.stimstart = 0;
% % % stim.stimstop = 1;
% % % stim.stimepoch = 'Go cue';
% % % 
% % % kinfeat = 'tongue_length';
% % % featix = find(strcmp(kin(sessix).featLeg,kinfeat));
% % % 
% % % cond2plot = [1 2 3 4];
% % % condition = '2AFC';
% % % figure();
% % % plotKinTracking_CtrlvsStim(params,obj,sessix,cols,kin,offset,smooth,kinfeat,featix,trix2plot,cond2plot,stim,condition);
% % % 
% % % cond2plot = [5 6 7 8];
% % % condition = 'AW';
% % % figure();
% % % plotKinTracking_CtrlvsStim(params,obj,sessix,cols,kin,offset,smooth,kinfeat,featix,trix2plot,cond2plot,stim,condition);
% % % %% Calculate percentage of time where tongue is visible during stim period
% % % cond2plot = 1:8;
% % % 
% % % stimepoch = 0;
% % % stim.stimstart = 0;
% % % stim.stimstop = 0.6;
% % % 
% % % for sessix = 1:length(meta)
% % %     pctTime = NaN(1,length(cond2plot));
% % % 
% % %     for cond = 1:length(cond2plot)
% % %         stimstart = stim.stimstart;
% % %         startix = find(obj(sessix).time>stimstart,1,'first');
% % %         stimstop = stim.stimstop;
% % %         stopix = find(obj(sessix).time<stimstop,1,'last');
% % % 
% % %         condtrix = params(sessix).trialid{cond2plot(cond)};
% % %         tempTime = NaN(1,length(condtrix));
% % %         for t = 1:length(condtrix)
% % %             currtrial = condtrix(t);
% % %             featix = find(strcmp(kin(sessix).featLeg,'tongue_length'));
% % %             tonguedat = kin(sessix).dat_std(startix:stopix,currtrial,featix);
% % %             dur = length(tonguedat);
% % %             tonguevis = length(find(~isnan(tonguedat)));
% % %             tempTime(t) = tonguevis/dur;
% % % 
% % %             %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %             %         offset = 5;
% % %             %         figure(cond)
% % %             %         toplot = offset*t+tonguedat;
% % %             %         if cond==1||cond==2
% % %             %             plottime = 0.5+obj(sessix).time(startix:stopix);
% % %             %         else
% % %             %             plottime = 1+obj(sessix).time(startix:stopix);
% % %             %         end
% % %             %         plot(plottime,toplot); hold on
% % %             %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         end
% % %         pctTime(cond) = mean(tempTime,'omitnan');
% % %         %title(num2str(pctTime(cond)))
% % %     end
% % %     stimEffects(sessix).pctTime = pctTime;
% % % end
% % % %% Make bar plot comparing pct of time with tongue visible for each condition
% % % pctTime = [];
% % % for sess = 1:length(stimEffects)
% % %     pctTime = [pctTime;stimEffects(sess).pctTime];
% % % end
% % % 
% % % cols = getColors();
% % % anmNames = {'MAH13','MAH13','MAH13','MAH13','MAH13','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14'};
% % % sigcutoff = 0.05;
% % % 
% % % subplot(1,2,1)
% % % condition = '2AFC';
% % % starheight = 0.4;
% % % temp = NaN(length(stimEffects),4);
% % % temp(:,1) = pctTime(:,1);
% % % temp(:,2) = pctTime(:,3);
% % % temp(:,3) = pctTime(:,2);
% % % temp(:,4) = pctTime(:,4);
% % % CtrlvsStimBarPlot(cols,temp,anmNames,sigcutoff,starheight,condition)
% % % ylabel('Fraction of time')
% % % ylim([0 0.42])
% % % title('2AFC trials')
% % % clear temp
% % % 
% % % subplot(1,2,2)
% % % condition = 'AW';
% % % starheight = 0.4;
% % % temp = NaN(length(stimEffects),4);
% % % temp(:,1) = pctTime(:,5);
% % % temp(:,2) = pctTime(:,7);
% % % temp(:,3) = pctTime(:,6);
% % % temp(:,4) = pctTime(:,8);
% % % CtrlvsStimBarPlot(cols,temp,anmNames,sigcutoff,starheight,condition)
% % % ylabel('Fraction of time')
% % % ylim([0 0.42])
% % % title('Autowater trials')
% % % clear temp
% % % sgtitle(['Frac of time with tongue visible--from ' num2str(stim.stimstart) ' to ' num2str(stim.stimstop) ' (s)'])
