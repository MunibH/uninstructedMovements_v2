function plotMagnitudeCluPotent_ME(obj,rez,params,me,cond2use)


X = obj(1).time;
e = [0.1 2]; % (s)
for i = 1:numel(e)
    [~,ix(i)] = min(abs(X - e(i)));
end
X = X(ix(1):ix(2));

normType = 'zscore';
for sessix = 1:numel(rez)
    trix = cell2mat(params(sessix).trialid(cond2use)');

    %     allme(:,sessix) = (nanmean(me(sessix).data(:,trix),2));
    %     allpotent(:,sessix) = (nanmean(sum(rez(sessix).N_potent(:,trix,:),3).^2,2));
    %     allclu(:,sessix) = (mean(sum(obj(sessix).trialdat(:,:,trix),2).^2,3));

    allme(:,sessix) = normalize(nanmean(me(sessix).data(:,trix),2),normType);
    allpotent(:,sessix) = normalize(nanmean(sum(rez(sessix).N_potent(:,trix,:),3).^2,2),normType);
    allclu(:,sessix) = normalize(mean(sum(obj(sessix).trialdat(:,:,trix),2).^2,3),normType);


    decay.me(sessix) = fitNegExp(X,zscore(allme(ix(1):ix(2),sessix)));
    decay.potent(sessix) = fitNegExp(X,zscore(allpotent(ix(1):ix(2),sessix)));
    decay.clu(sessix) = fitNegExp(X,zscore(allclu(ix(1):ix(2),sessix)));
end


f = figure;
cols = getColors;
cols.clu = [153, 46, 199] ./ 255;
alph = 0.2;
tm = obj(1).time;
xlims = [-2.4 2];

sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;
gc = mode(obj(1).bp.ev.goCue) - 2.5;

ax = subplot(3,1,1);
dat = allme;
shadedErrorBar(tm, mean(dat,2), std(dat,[],2)./sqrt(numel(rez)), {'Color','k','LineWidth',1},alph,ax); xlim(xlims)
xline(sample,'k--','LineWidth',1); xline(delay,'k--','LineWidth',1); xline(gc,'k--','LineWidth',1)
ylabel('Motion Energy (a.u.)')
ax = subplot(3,1,2);
dat = allpotent;
shadedErrorBar(tm, mean(dat,2), std(dat,[],2)./sqrt(numel(rez)), {'Color',cols.potent,'LineWidth',1},alph,ax); xlim(xlims)
xline(sample,'k--','LineWidth',1); xline(delay,'k--','LineWidth',1); xline(gc,'k--','LineWidth',1)
ylabel('Potent - Magnitude (a.u.)')
ax = subplot(3,1,3);
dat = allclu;
shadedErrorBar(tm, mean(dat,2), std(dat,[],2)./sqrt(numel(rez)), {'Color',cols.clu,'LineWidth',1},alph,ax); xlim(xlims)
xline(sample,'k--','LineWidth',1); xline(delay,'k--','LineWidth',1); xline(gc,'k--','LineWidth',1)
ylabel('Full Neural Pop - Magnitude (a.u.)')
xlabel('Time from go cue (s)')


xs = [ 1 2 3 ];
clrs{1} = [0.5 0.5 0.5];
clrs{2} = cols.potent;
clrs{3} = cols.clu;

fns = fieldnames(decay);

f=figure;
ax = gca;
hold on;
xs = [1 2 3];
for i = 1:numel(fns)
    temp = decay.(fns{i});
    h(i) = bar(xs(i),mean(temp)); % mean across dims and sessions
    cix = i;
    h(i).FaceColor = clrs{i};
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.5;
    scatter(xs(i)*ones(size(temp)),temp,20,'MarkerFaceColor',[0 0 0], ...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
        'MarkerFaceAlpha',0.7)
    errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);
end
ax.XTick = xs;
xticklabels(fns);

ylabel('Decay term (neg exp fit)')
ax.FontSize = 12;


f = figure;
alph = 0.2;
tm = obj(1).time;

ax = subplot(3,1,1);
dat = allme;
plot(tm, dat, 'Color','k','LineWidth',1); xlim(xlims)
xline(sample,'k--','LineWidth',1); xline(delay,'k--','LineWidth',1); xline(gc,'k--','LineWidth',1)
ylabel('Motion Energy (a.u.)')
ax = subplot(3,1,2);
dat = allpotent;
plot(tm, dat, 'Color',cols.potent,'LineWidth',1); xlim(xlims)
xline(sample,'k--','LineWidth',1); xline(delay,'k--','LineWidth',1); xline(gc,'k--','LineWidth',1)
ylabel('Potent (a.u.)')
ax = subplot(3,1,3);
dat = allclu;
plot(tm, dat, 'Color',cols.clu,'LineWidth',1); xlim(xlims)
xline(sample,'k--','LineWidth',1); xline(delay,'k--','LineWidth',1); xline(gc,'k--','LineWidth',1)
ylabel('Full Pop (a.u.)')





end