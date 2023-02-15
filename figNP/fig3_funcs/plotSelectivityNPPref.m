function plotSelectivityNPPref(meta,obj,params,rez,cond2use_ta,cond2use_st,subTrials)

cols = getColors;
lw = 2;
alph = 0.2;

mainCond = 1; % right trials should always be higher

%%

% % Sort by delay selectivity
for i = 1:numel(obj)
    edges{i} = [mode(obj(i).bp.ev.delay) mode(obj(i).bp.ev.goCue)] - mode(obj(i).bp.ev.goCue);
end


modparams.subTrials = subTrials;

[plotsel.null, plotsel.nullDims] = getSortedSelectivityPref(obj,params,rez,edges,cond2use_ta,cond2use_st,modparams,'null');
[plotsel.potent, plotsel.potentDims] = getSortedSelectivityPref(obj,params,rez,edges,cond2use_ta,cond2use_st,modparams,'potent');

%% plot hits

xlims = [-2.3 2.5];

temp = plotsel.null;
maxSelect =  max(abs(temp),[],1);                    % Max magnitude of selectivity for each dim
temp = temp ./ maxSelect;

f = figure; f.Position = [794   516   302   384];
ax = gca;
hold on;
imagesc(obj(1).time, 1:size(temp,2),temp'); 
c = colorbar;
colormap((redblue))
xlabel('Time (s) from go cue')
ylabel('Dimensions')
c.Label.String = 'Normalized Selectivity (a.u.)';
title('null')
ax.YDir = 'reverse';
xlim(xlims)
ylim([1 size(plotsel.null,2)])
yline(plotsel.nullDims(end),'k--','LineWidth',1)

temp = plotsel.potent;
maxSelect =  max(abs(temp),[],1);                    % Max magnitude of selectivity for each dim
temp = temp ./ maxSelect;

f = figure; f.Position = [794   516   302   384];
ax = gca;
hold on;
imagesc(obj(1).time, 1:size(temp,2),temp'); 
c = colorbar;
colormap((redblue))
xlabel('Time from go cue (s)')
ylabel('Dimensions')
c.Label.String = 'Normalized Selectivity (a.u.)';
title('potent')
ax.YDir = 'reverse';
xlim(xlims)
ylim([1 size(plotsel.potent,2)])
yline(plotsel.potentDims(end),'k--','LineWidth',1)




%% mean squared selectivity

col = getColors;
sm = 41;
lw = 2;
alph = 0.1;

smtype = 'reflect';

sample = mode(obj(1).bp.ev.sample - 2.5);
delay = mode(obj(1).bp.ev.delay - 2.5);
tstart = mode(obj(1).bp.ev.bitStart - 2.5);

% dat.null = mySmooth(plotsel.null(:,plotsel.nullDims),sm,smtype);
% dat.potent = mySmooth(plotsel.potent(:,plotsel.potentDims),sm,smtype);
dat.null = mySmooth(plotsel.null,sm,smtype);
dat.potent = mySmooth(plotsel.potent,sm,smtype);

clrs{1} = col.null;
clrs{2} = col.potent;

fns = fieldnames(dat);
f = figure; 
f.Position = [88   603   404   256];
ax = gca;
hold on;
for i = 1:numel(fns)
    temp = dat.(fns{i});
    mu = nanmean(temp,2);
    sig = nanstd(temp,[],2) ./ sqrt(size(temp,2));
    shadedErrorBar(obj(1).time,mu,getCI(temp),{'Color',clrs{i},'LineWidth',2,'LineStyle','-'},alph,ax);
    xline(0,'k--')
    xline(delay,'k--')
    xline(sample,'k--')
%     yline(0,'k-')
end
xlim(xlims)
xlabel('Time from go cue (s)')
ylabel('Selectivity (a.u.)')


end












