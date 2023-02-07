close all;

f = figure;

cond2use = 8:11;

col = getColors;
cols{1} = col.rhit;
cols{2} = col.lhit;
cols{3} = col.rmiss;
cols{4} = col.lmiss;

lw = 1;
alph = 0.17;

sm = 21;

for i = 1:numel(meta)
    ax = nexttile;
    hold on;
    for c = 1:numel(cond2use)
        trix = params(i).trialid{cond2use(c)};
        temp = mySmooth(me(i).data(:,trix),sm,'reflect');
        plot(obj(1).time,nanmean(temp,2),'Color',cols{c},'LineWidth',lw)
        %         shadedErrorBar(obj(1).time,nanmean(temp,2),nanstd(temp,[],2)./sqrt(size(temp,2)), {'Color',cols{c},'LineWidth',lw},alph,ax)
        allme(:,c,i) = normalize(nanmean(temp,2),'range',[0 1]);
    end
end

f = figure;
ax = gca;
hold on;
alph = 0.2;
for c = 1:2
    temp = squeeze(allme(:,c,:));
    mu = nanmean(temp,2);
    sigma = nanstd(temp,[],2) ./ sqrt(numel(obj));
    shadedErrorBar(obj(1).time,mu, sigma,{ 'Color',cols{c},'LineWidth',2},alph,ax)
end
xlim([-2.35 1])

f = figure;
ax = gca;
hold on;
alph = 0.2;
for c = 3:4
    temp = squeeze(allme(:,c,:));
    mu = nanmean(temp,2);
    sigma = nanstd(temp,[],2) ./ sqrt(numel(obj));
    shadedErrorBar(obj(1).time,mu, sigma,{ 'Color',cols{c},'LineWidth',2},alph,ax)
end
xlim([-2.35 1])

%%

f = figure;

cond2use = 2:5;
col = getColors;
cols{1} = col.rhit;
cols{2} = col.lhit;
cols{3} = col.rmiss;
cols{4} = col.lmiss;

lw = 1;
alph = 0.17;

sm = 21;

for sessix = 1:numel(meta)
    ax = nexttile;
    hold on;


    t1 = -0.02; t2 = 0;
    [~,ix1] = min(abs(obj(sessix).time - t1));
    [~,ix2] = min(abs(obj(sessix).time - t2));

    me_ = nanmean(me(sessix).data(ix1:ix2,:),1);

    % right and left hit trials
    for i = 1:numel(cond2use)
        trix = params(sessix).trialid{cond2use(i)};
        temp = me_(trix);
        h(i) = bar(i,mean(temp),'EdgeColor','none','FaceColor',cols{i});
        errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);

    end

    

end


%% group me and CDPotent by animal

[ix,anms] = groupSessionsByAnimal(meta);

close all;

cond2use = 8:11;

xlims = [-2.3 1.5];

col = getColors;
cols{1} = col.rhit;
cols{2} = col.lhit;
cols{3} = col.rmiss;
cols{4} = col.lmiss;

lw = 1.5;
alph = 0.17;

sm = 35;

sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;

for i = 1:numel(ix) % for each animal
    f = figure;
    f.Position = [680   422   381   556];

    ax = subplot(2,1,1); hold on;

    sessions = find(ix{i});
    for sessix = 1:numel(sessions)
        sess = sessions(sessix);
        thispotent(:,:,sessix) = squeeze(cd_potent_all.cd_proj(:,:,1,sess)); % (time,cond) % cd late
        thispotent(:,:,sessix) = mySmooth(thispotent(:,:,sessix),sm,'reflect');

        thisnull(:,:,sessix) = squeeze(cd_null_all.cd_proj(:,:,1,sess)); % (time,cond) % cd late
        thisnull(:,:,sessix) = mySmooth(thisnull(:,:,sessix),sm,'reflect');

        for c = 1:numel(cond2use)
            trix = params(sess).trialid{cond2use(c)};
            thisme(:,c,sessix) = nanmean(me(sess).data(:,trix),2);
            thisme(:,c,sessix) = mySmooth(thisme(:,c,sessix),sm,'reflect');
        end
    end
    % mean across sessions for current animal
    thispotent = nanmean(thispotent,3);
    thisnull = nanmean(thisnull,3);
    thisme = nanmean(thisme,3);

    % smooth error trials some more
    thispotent(:,3:4) = mySmooth(thispotent(:,3:4),sm,'reflect');
    thisnull(:,3:4) = mySmooth(thisnull(:,3:4),sm,'reflect');
    thisme(:,3:4) = mySmooth(thisme(:,3:4),sm,'reflect');

    ax1 = subplot(3,1,1); hold on;
    ax2 = subplot(3,1,2); hold on;
    ax3 = subplot(3,1,3); hold on;
    for c = 1:numel(cond2use)
        plot(ax1,obj(1).time,thisnull(:,c),'Color',cols{c},'LineWidth',lw)
        plot(ax2,obj(1).time,thispotent(:,c),'Color',cols{c},'LineWidth',lw)
        plot(ax3,obj(1).time,thisme(:,c),'Color',cols{c},'LineWidth',lw)
    end
    xlim(ax1,xlims)
    xlim(ax2,xlims)
    xlim(ax3,xlims)
    xline(ax1,0,'k--'); xline(ax2,0,'k--'); xline(ax3,0,'k--')
    xline(ax1,sample,'k--'); xline(ax2,sample,'k--'); xline(ax3,sample,'k--')
    xline(ax1,delay,'k--'); xline(ax2,delay,'k--'); xline(ax3,delay,'k--')
    

    sgtitle(anms{i})

end

