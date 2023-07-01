function plotSelectivityNP(meta,obj,params,rez,cond2use_ta,cond2use_st)

cols = getColors;
lw = 2;
alph = 0.2;

mainCond = 1; % right trials should always be higher

%%

% Sort by delay selectivity
edges = [mode(obj(1).bp.ev.delay) mode(obj(1).bp.ev.goCue)] - 2.5;
% edges = [-0.8 0]; % (s) relative to go cue

% sort by late-sample selectivity
% edges = [mode(obj(1).bp.ev.delay)-0.5 mode(obj(1).bp.ev.delay)] - 2.5;

% % sort by presample-selectivity
% edges = [mode(obj(1).bp.ev.sample)-0.3 mode(obj(1).bp.ev.sample)] - 2.5;

% % sort by response-selectivity
% edges = [mode(obj(1).bp.ev.goCue)+0.02 mode(obj(1).bp.ev.goCue)+0.5] - 2.5;

modparams.subTrials = 35;

[plotsel.null, plotsel.null_sep, plotsel.null_sel_dims] = getSortedSelectivity(obj,params,rez,edges,cond2use_ta,cond2use_st,modparams,'null');
[plotsel.potent, plotsel.potent_sep, plotsel.potent_sel_dims] = getSortedSelectivity(obj,params,rez,edges,cond2use_ta,cond2use_st,modparams,'potent');

%% plot

xlims = [-2.3 2.5];

f = figure; f.Position = [794   516   302   384];
ax = gca;
hold on;
imagesc(obj(1).time, 1:size(plotsel.null,2),plotsel.null'); 
c = colorbar;
colormap(flipud(redblue))
xlabel('Time (s) from go cue')
ylabel('Dimensions')
c.Label.String = 'Selectivity (a.u.)';
title('null')
ax.YDir = 'reverse';
xlim(xlims)
ylim([1 size(plotsel.null,2)])
yline(plotsel.null_sep(1),'k--','LineWidth',2)
yline(plotsel.null_sep(2),'k--','LineWidth',2)


f = figure; f.Position = [794   516   302   384];
ax = gca;
hold on;
imagesc(obj(1).time, 1:size(plotsel.potent,2),plotsel.potent'); 
c = colorbar;
colormap(flipud(redblue))
xlabel('Time (s) from go cue')
ylabel('Dimensions')
c.Label.String = 'Selectivity (a.u.)';
title('potent')
ax.YDir = 'reverse';
xlim(xlims)
ylim([1 size(plotsel.potent,2)])
yline(plotsel.potent_sep(1),'k--','LineWidth',2)
yline(plotsel.potent_sep(2),'k--','LineWidth',2)



%% mean squared selectivity

col = getColors;
sm = 31;
lw = 2;
alph = 0.1;

smtype = 'zeropad';

sample = mode(obj(1).bp.ev.sample - 2.5);
delay = mode(obj(1).bp.ev.delay - 2.5);
tstart = mode(obj(1).bp.ev.bitStart - 2.5);

null = mySmooth(plotsel.null(:,1:plotsel.null_sep(2)),sm,smtype);
potent = mySmooth(plotsel.potent(:,1:plotsel.potent_sep(2)),sm,smtype);

% % null = normalize(null,'zscore');
% % potent = normalize(potent,'zscore');
% 
% % null = mySmooth(plotsel.null,sm);
% % potent = mySmooth(plotsel.potent,sm);
% 
% % % squared selectivity
% % null = abs(null);
% % potent = abs(potent);
% % null = null.^2;
% % potent = potent.^2;
% 
% mu.null = mean(null,2);
% mu.potent = mean(potent,2);
% 
% sig.null = std(null,[],2) ./ sqrt(size(null,2));
% sig.potent = std(potent,[],2) ./ sqrt(size(potent,2));
% 
% f = figure; 
% ax = gca;
% hold on;
% shadedErrorBar(obj(1).time,mu.null,sig.null,{'Color',col.null,'LineWidth',2,'LineStyle','-'},alph,ax);
% shadedErrorBar(obj(1).time,mu.potent,sig.potent,{'Color',col.potent,'LineWidth',2,'LineStyle','-'},alph,ax);
% xlim(xlims)
% xline(0,'k--')
% xline(delay,'k--')
% xline(sample,'k--')
% yline(0,'k-')
% 
% xlabel('Time (s) from go cue')
% ylabel('Selectivity (a.u.)')


%%
clear null potent mu sig

xlims = [tstart 1];

sm = 31;

null.right = mySmooth(plotsel.null(:,1:plotsel.null_sep(1)),sm,smtype);
null.left = mySmooth(plotsel.null(:,plotsel.null_sep(1)+1:plotsel.null_sep(2)),sm,smtype);

potent.right = mySmooth(plotsel.potent(:,1:plotsel.potent_sep(1)),sm,smtype);
potent.left = mySmooth(plotsel.potent(:,plotsel.potent_sep(1)+1:plotsel.potent_sep(2)),sm,smtype);

mu.null.right = mean(null.right,2);
mu.null.left = mean(null.left,2);

mu.potent.right = mean(potent.right,2);
mu.potent.left = mean(potent.left,2);

sig.null.right = std(null.right,[],2) ./ sqrt(size(null.right,2));
sig.null.left = std(null.left,[],2) ./ sqrt(size(null.left,2));

sig.potent.right = std(potent.right,[],2) ./ sqrt(size(potent.right,2));
sig.potent.left = std(potent.left,[],2) ./ sqrt(size(potent.left,2));

% gradients
dt = params(1).dt;
mu.grad.null.right = gradient(mu.null.right,dt);
mu.grad.null.left = gradient(mu.null.left,dt);

mu.grad.potent.right = gradient(mu.potent.right,dt);
mu.grad.potent.left = gradient(mu.potent.left,dt);


% plot

f = figure; 
ax1 = subplot(2,1,1);
hold on;
shadedErrorBar(obj(1).time,mu.null.right,sig.null.right,{'Color',col.rhit,'LineWidth',2,'LineStyle','-'},alph,ax1);
shadedErrorBar(obj(1).time,mu.null.left,sig.null.left,{'Color',col.lhit,'LineWidth',2,'LineStyle','-'},alph,ax1);
xlim(xlims)
xline(0,'k--')
xline(delay,'k--')
xline(sample,'k--')
yline(0,'k-')
title('Null | Right and Left Selective Dims')
ylabel('Activity (a.u.)')

% ax2 = subplot(4,1,2);
% hold on;
% plot(obj(1).time,mu.grad.null.right,'Color',col.rhit,'LineWidth',2);
% plot(obj(1).time,mu.grad.null.left,'Color',col.lhit,'LineWidth',2);
% xlim(xlims)
% xline(0,'k--')
% xline(delay,'k--')
% xline(sample,'k--')
% yline(0,'k-')
% ylabel('dA/dt (a.u.)')
% xlabel('Time (s) from go cue')

ax3 = subplot(2,1,2);
hold on;
shadedErrorBar(obj(1).time,mu.potent.right,sig.potent.right,{'Color',col.rhit,'LineWidth',2,'LineStyle','-'},alph,ax3);
shadedErrorBar(obj(1).time,mu.potent.left,sig.potent.left,{'Color',col.lhit,'LineWidth',2,'LineStyle','-'},alph,ax3);
xlim(xlims)
xline(0,'k--')
xline(delay,'k--')
xline(sample,'k--')
yline(0,'k-')
title('Potent | Right and Left Selective Dims')
ylabel('Activity (a.u.)')

% ax4 = subplot(4,1,4);
% hold on;
% plot(obj(1).time,mu.grad.potent.right,'Color',col.rhit,'LineWidth',2);
% plot(obj(1).time,mu.grad.potent.left,'Color',col.lhit,'LineWidth',2);
% xlim(xlims)
% xline(0,'k--')
% xline(delay,'k--')
% xline(sample,'k--')
% yline(0,'k-')
% title('Potent | Derivative')
% ylabel('dA/dt (a.u.)')
% xlabel('Time (s) from go cue')

% match y axis limits across null and potent space selectivity
ys = [ min(ax1.YLim(1),ax3.YLim(1)) max(ax1.YLim(2),ax3.YLim(2)) ];
ax1.YLim = ys;
ax3.YLim = ys;

% ys = [ min(ax2.YLim(1),ax4.YLim(1)) max(ax2.YLim(2),ax4.YLim(2)) ];
% ax2.YLim = ys;
% ax4.YLim = ys;

%% selectivity correlation matrix

% sel = plotsel.null;
% 
% sel_corr_mat.null = zeros(size(sel,1),size(sel,1));
% for i = 1:size(sel_corr_mat.null,1)
%     for j = 1:size(sel_corr_mat.null,1)
%         temp = corrcoef(sel(i,:),sel(j,:));
%         sel_corr_mat.null(i,j) = temp(1,2);
%     end
% end
% 
% sel = plotsel.potent;
% 
% sel_corr_mat.potent = zeros(size(sel,1),size(sel,1));
% for i = 1:size(sel_corr_mat.potent,1)
%     for j = 1:size(sel_corr_mat.potent,1)
%         temp = corrcoef(sel(i,:),sel(j,:));
%         sel_corr_mat.potent(i,j) = temp(1,2);
%     end
% end
% 
% 
% 
% f = figure;
% f.Position = [400   400   332   250];
% ax = gca;
% hold on;
% imagesc(obj(1).time,obj(1).time,sel_corr_mat.null)
% title('Null')
% c = colorbar;
% xlim(xlims)
% ylim(xlims)
% colormap(summer)
% xline(0,'k--')
% xline(delay,'k--')
% xline(sample,'k--')
% yline(0,'k--')
% yline(delay,'k--')
% yline(sample,'k--')
% 
% f = figure;
% f.Position = [680   400   332   250];
% ax = gca;
% hold on;
% imagesc(obj(1).time,obj(1).time,sel_corr_mat.potent)
% title('Potent')
% c = colorbar;
% xlim(xlims)
% ylim(xlims)
% colormap(spring)
% xline(0,'k--')
% xline(delay,'k--')
% xline(sample,'k--')
% yline(0,'k--')
% yline(delay,'k--')
% yline(sample,'k--')



end












