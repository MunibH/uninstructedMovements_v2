[objix,uAnm] = groupSessionsByAnimal(meta);

sessixs = find(objix{3}); % jeb13
xlims = [-2.4 2];


cond2use = [15 16];
cols = getColors;
col{1} = cols.rhit_aw;
col{2} = cols.lhit_aw;
alph = 0.2;
f = figure;
for i = 1:numel(obj)
    s = i;

    sample = mode(obj(s).bp.ev.sample) - mode(obj(s).bp.ev.goCue);
    delay = mode(obj(s).bp.ev.delay) - mode(obj(s).bp.ev.goCue);
    gc = 0;


    ax = nexttile;
    hold on;
    for c = 1:numel(cond2use)
        trix = params(s).trialid{cond2use(c)};
        temp = me(s).data(:,trix);
        mu = nanmean(temp,2);
        sig = nanstd(temp,[],2)./sqrt(numel(trix));
        shadedErrorBar(obj(1).time,mu,sig,{'Color',col{c},'LineWidth',1,'LineStyle','-'},alph,ax)
%         plot(obj(1).time,mu,'Color',col{c},'LineWidth',1);
        
    end
    xline(sample,'k--')
    xline(delay,'k--')
    xline(gc,'k--')
    xlim(xlims)
    title([meta(s).anm ' ' meta(s).date])
end








