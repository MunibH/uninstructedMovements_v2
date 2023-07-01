function [null,potent,rightm,leftm] = plotSelectivityNPPrefOnset(meta,obj,params,me,rez,cond2use_ta,cond2use_st,subTrials)

cols = getColors;
lw = 2;
alph = 0.2;

mainCond = 1; % right trials should always be higher

%%

% % Sort by delay selectivity
edges = [mode(obj.bp.ev.delay) mode(obj.bp.ev.goCue)] - mode(obj.bp.ev.goCue);
% edges = [-0.4 0]; % (s) relative to go cue

% sort by late-sample selectivity
% edges = [mode(obj(1).bp.ev.delay)-0.5 mode(obj(1).bp.ev.delay)] - mode(obj.bp.ev.goCue);

% % sort by presample-selectivity
% edges = [mode(obj(1).bp.ev.sample)-0.3 mode(obj(1).bp.ev.sample)] - mode(obj.bp.ev.goCue);

% % sort by response-selectivity
% edges = [mode(obj(1).bp.ev.goCue)+0.02 mode(obj(1).bp.ev.goCue)+0.5] - mode(obj.bp.ev.goCue);

modparams.subTrials = subTrials;

[plotsel.null, plotsel.nullDims] = getSortedSelectivityPref(obj,params,rez,edges,cond2use_ta,cond2use_st,modparams,'null');
[plotsel.potent, plotsel.potentDims] = getSortedSelectivityPref(obj,params,rez,edges,cond2use_ta,cond2use_st,modparams,'potent');

%%
trix{1} = params.trialid{cond2use_st(1)};
trix{2} = params.trialid{cond2use_st(2)};
m{1} = me.data(:,trix{1});
m{2} = me.data(:,trix{2});


%%

null = nanmean(plotsel.null(:,plotsel.nullDims),2);
potent = nanmean(plotsel.potent(:,plotsel.potentDims),2);
rightm = nanmean(m{1},2);
leftm = nanmean(m{2},2);


%% mean squared selectivity

% xlims = [-2.4 2];
% 
% col = getColors;
% sm = 41;
% lw = 2;
% alph = 0.1;
% 
% smtype = 'reflect';
% 
% sample = mode(obj(1).bp.ev.sample - 2.5);
% delay = mode(obj(1).bp.ev.delay - 2.5);
% tstart = mode(obj(1).bp.ev.bitStart - 2.5);
% 
% dat.null = mySmooth(plotsel.null(:,plotsel.nullDims),sm,smtype);
% dat.potent = mySmooth(plotsel.potent(:,plotsel.potentDims),sm,smtype);
% % dat.null = mySmooth(plotsel.null,sm,smtype);
% % dat.potent = mySmooth(plotsel.potent,sm,smtype);
% 
% clrs{1} = col.null;
% clrs{2} = col.potent;
% 
% fns = fieldnames(dat);
% f = figure; 
% f.Position = [88   603   404   256];
% ax = gca;
% hold on;
% yyaxis(ax,'left');
% for i = 1:numel(fns)
%     temp = dat.(fns{i});
%     mu = nanmean(temp,2);
%     sig = nanstd(temp,[],2) ./ sqrt(size(temp,2));
% %     shadedErrorBar(obj(1).time,mu,getCI(temp),{'Color',clrs{i},'LineWidth',2,'LineStyle','-'},alph,ax);
%     plot(obj(1).time,mu,'Color',clrs{i},'LineWidth',2,'LineStyle','-');
%     xline(0,'k--')
%     xline(delay,'k--')
%     xline(sample,'k--')
% end
% yline(0,'k-')
% xlim(xlims)
% xlabel('Time from go cue (s)')
% ylabel('Selectivity (a.u.)')
% 
% clrs{1} = col.rhit;
% clrs{2} = col.lhit;
% 
% yyaxis(ax,'right');
% temp = me_{1};
% mu = nanmean(temp,2);
% sig = nanstd(temp,[],2) ./ sqrt(size(temp,2));
% % shadedErrorBar(obj(1).time,mu,getCI(temp),{'Color',clrs{1},'LineWidth',2,'LineStyle','-'},alph,ax);
% plot(obj(1).time,mu,'Color',clrs{1},'LineWidth',2,'LineStyle','-');
% temp = me_{2};
% mu = nanmean(temp,2);
% sig = nanstd(temp,[],2) ./ sqrt(size(temp,2));
% % shadedErrorBar(obj(1).time,mu,getCI(temp),{'Color',clrs{2},'LineWidth',2,'LineStyle','-'},alph,ax);
% plot(obj(1).time,mu,'Color',clrs{2},'LineWidth',2,'LineStyle','-');
% ylabel('Motion energy (a.u.)')

end












