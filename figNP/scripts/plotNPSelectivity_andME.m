
close all

dedges = [-0.7 0];


xlims = [-2.3 0];

cond2use = [8 9];
cols = getColors;
col{1} = cols.rhit;
col{2} = cols.lhit;
alph = 0.2;
f = figure;
t = tiledlayout('flow');
for i = 1:numel(obj)
    s = i;

    sample = mode(obj(s).bp.ev.sample) - mode(obj(s).bp.ev.goCue);
    delay = mode(obj(s).bp.ev.delay) - mode(obj(s).bp.ev.goCue);
    gc = 0;

    for j = 1:numel(dedges)
        [~,dix(j)] = min(abs(obj(s).time - dedges(j)));
    end

    tempme{1} = me_.right(:,s);
    tempme{2} = me_.left(:,s);

    memu{1} = nanmean(tempme{1}(dix(1):dix(2)),1);
    memu{2} = nanmean(tempme{2}(dix(1):dix(2)),1);
    if memu{1} > memu{2}
        pref = 1;
        nonpref = 2;
    else
        pref = 2;
        nonpref = 1;
    end


    ax = nexttile;
    yyaxis(ax,'left')
    hold on;
    plot(obj(1).time,tempme{pref} - tempme{nonpref},'Color','k','LineWidth',1,'LineStyle','-');
%     for c = 1:numel(cond2use)
%         trix = params(s).trialid{cond2use(c)};
%         temp = me(s).data(:,trix);
%         mu = nanmean(temp,2);
%         sig = nanstd(temp,[],2)./sqrt(numel(trix));
% %         shadedErrorBar(obj(1).time,mu,sig,{'Color',col{c},'LineWidth',1,'LineStyle','-'},alph,ax)
%         plot(obj(1).time,mu,'Color',col{c},'LineWidth',1,'LineStyle','-');
%     end
    yyaxis(ax,'right');
    plot(obj(1).time,mySmooth(potent(:,s),21),'Color',cols.potent,'LineWidth',1.5,'LineStyle','-');
%     plot(obj(1).time,mySmooth(null(:,s),21),'Color',cols.null,'LineWidth',1.5,'LineStyle','-');
    yline(0,'k--','LineWidth',1.5)

    xline(sample,'k--')
    xline(delay,'k--')
    xline(gc,'k--')
    xlim(xlims)
    title([meta(s).anm ' ' meta(s).date])
end

xlabel(t,'Time from go cue (s) ')

















