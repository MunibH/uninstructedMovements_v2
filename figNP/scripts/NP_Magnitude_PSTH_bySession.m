close all
clear col cols
cols = getColors;
col.null{1} = [3, 3, 173] ./ 255;
col.null{2} = [173, 3, 3] ./ 255;

col.potent{1} = [115, 169, 245] ./ 255;
col.potent{2} = [240, 110, 138] ./ 255;

alph = 0.2;
lw = 2;

sm = 31;

xlims = [-2.3,2];

sample = (mode(obj(1).bp.ev.sample) - 2.5);
delay = (mode(obj(1).bp.ev.delay) - 2.5);


cond2use = [8 9];
for sessix = 1:numel(meta)
    f = figure;
    ax = gca;
    hold on;
    for c = 1:numel(cond2use)
        trix = params(sessix).trialid{cond2use(c)};

        null = mySmooth(sum(rez(sessix).N_null(:,trix,:).^2,3),sm,'reflect');
        mu = nanmean(null,2);
        sig = nanstd(null,[],2) ./ sqrt(numel(trix));
        shadedErrorBar(obj(1).time,mu,getCI(null),{'Color',col.null{c},'LineWidth',lw},alph,ax)

        potent = mySmooth(sum(rez(sessix).N_potent(:,trix,:).^2,3),sm,'reflect');
        mu = nanmean(potent,2);
        sig = nanstd(potent,[],2) ./ sqrt(numel(trix));
        shadedErrorBar(obj(1).time,mu,getCI(potent),{'Color',col.potent{c},'LineWidth',lw},alph,ax)
    end
    title([meta(sessix).anm ' ' meta(sessix).date])
    ax.FontSize = 12;
    xline(sample,'k--')
    xline(delay,'k--')
    xline(0,'k--')
    xlim(xlims)
    xlabel('Time from go cue (s)')
    ylabel('Activity')
end





