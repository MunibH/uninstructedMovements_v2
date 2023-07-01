
% cond2plot{1} = [8 10]; % right hit, miss
% cond2plot{2} = [9 11]; % left hit, miss

cols = getColors;
col{1} = cols.rhit;
col{2} = cols.lhit;
col{3} = cols.rhit * 0.5;
col{4} = cols.lhit * 0.5;

cond2use = 8:11;
for sessix = 1:numel(meta)
    me_ = me(sessix).data;
    for c = 1:numel(cond2use)
        trix = params(sessix).trialid{cond2use(c)};
        dat(:,c,sessix) = nanmean(me_(:,trix),2); % (time,cond,sess)
    end

end

%%

sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;


lw = 1;
alph = 0.2;

xlims = [-2.3 2];

f = figure;
ax = gca;
hold on;
for i = 1:2
    dat_ = squeeze(dat(:,i,:));
    shadedErrorBar(obj(1).time,nanmean(dat_,2), nanstd(dat_,[],2)./sqrt(size(dat_,2)), {'LineWidth',lw,'Color',col{i}}, alph, ax )
end
xlim(xlims)
xline(0,'k--');
xline(sample,'k--');
xline(delay,'k--');
ylabel('ME')

f = figure;
ax = gca;
hold on;
for i = 3:4
    dat_ = squeeze(dat(:,i,:));
    shadedErrorBar(obj(1).time,nanmean(dat_,2), nanstd(dat_,[],2)./sqrt(size(dat_,2)), {'LineWidth',lw,'Color',col{i}}, alph, ax )
end
xlim(xlims)
xline(0,'k--');
xline(sample,'k--');
xline(delay,'k--');
ylabel('ME')























