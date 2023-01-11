function plotSubpsaceMagnitude_TrialAvg(obj,meta,params,rez,cond2use)

close all

cols = getColors;
clrs{1} = cols.rhit;
clrs{2} = cols.lhit;

lw = 2;
alph = 0.2;

for sessix = 1:numel(meta)

    % Null space
    f = figure;
    f.Position = [698   436   343   230];
    ax = gca;
    hold on;
    for cond = 1:numel(cond2use)
        trix = params(sessix).trialid{cond2use(cond)};
        tempdat = rez(sessix).N_null(:,trix,:);
        tempdat = sum(tempdat.^2,3); % sum across dims
        datmean = mean(tempdat,2);   % mean across trials
        datstd = std(tempdat,[],2) ./ sqrt(size(tempdat,2));
        shadedErrorBar(obj(sessix).time,datmean,datstd,{'Color',clrs{cond},'LineWidth',2},alph,ax)

        nulldat(:,sessix,cond) = datmean;
    end
    sample = mode(obj(sessix).bp.ev.sample) - 2.5;
    delay = mode(obj(sessix).bp.ev.delay) - 2.5;
    xline(0,'k--')
    xline(sample,'k:')
    xline(delay,'k:')
    xlim([-2.5 2.5])
    title(['Null | ' meta(sessix).anm ' ' meta(sessix).date])
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    ax.FontSize = 10;
    axnull = ax;

    % Potent space
    f = figure;
    f.Position = [1041         437         343         230];
    ax = gca;
    hold on;
    for cond = 1:numel(cond2use)
        trix = params(sessix).trialid{cond2use(cond)};
        tempdat = rez(sessix).N_potent(:,trix,:);
        tempdat = sum(tempdat.^2,3); % sum across dims
        datmean = mean(tempdat,2);   % mean across trials
        datstd = std(tempdat,[],2) ./ sqrt(size(tempdat,2));
        shadedErrorBar(obj(sessix).time,datmean,datstd,{'Color',clrs{cond},'LineWidth',2},alph,ax)

        potentdat(:,sessix,cond) = datmean;
    end
    sample = mode(obj(sessix).bp.ev.sample) - 2.5;
    delay = mode(obj(sessix).bp.ev.delay) - 2.5;
    xline(0,'k--')
    xline(sample,'k--')
    xline(delay,'k--')
    xlim([-2.5 2.5])
    title(['Potent | ' meta(sessix).anm ' ' meta(sessix).date])
    xlabel('Time (s) from go cue')
    ylabel('Activity (a.u.)')
    ax.FontSize = 10;
    axpotent = ax;

    % match y axis limits across null and potent space CDs
    %     ys = [ min(axnull.YLim(1),axpotent.YLim(1)) max(axnull.YLim(2),axpotent.YLim(2)) ];
    %     axnull.YLim = ys;
    %     axpotent.YLim = ys;


end

%% plot heatmaps of sessions
dat = cat(2, nulldat(:,:,1), nulldat(:,:,2));
f = figure;
f.Position = [318   531   348   415];
ax = gca;
hold on;
imagesc(obj(1).time, 1:size(dat,2), dat')
ax.XLim(1) = -2.5;
ax.XLim(2) = 2.5;
ylim([0.5 size(dat,2)+0.5])
xline(0,'k--')
xline(sample,'k--')
xline(delay,'k--')
yline(size(dat,2)/2,'k--')
title('Null')
xlabel('Time (s) from go cue')
ylabel('Sessions')
cnull = colorbar;
cnull.Label.String = 'Activity (a.u.)';
ax.FontSize = 10;
% colormap(linspecer)

dat = cat(2, potentdat(:,:,1), potentdat(:,:,2));
f = figure;
f.Position = [318   531   348   415];
ax = gca;
hold on;
imagesc(obj(1).time, 1:size(dat,2), dat')
ax.XLim(1) = -2.5;
ax.XLim(2) = 2.5;
ylim([0.5 size(dat,2)+0.5])
xline(0,'k--')
xline(sample,'k--')
xline(delay,'k--')
yline(size(dat,2)/2 + 0.5,'k--')
title('Potent')
xlabel('Time (s) from go cue')
ylabel('Sessions')
cpotent = colorbar;
cpotent.Label.String = 'Activity (a.u.)';
ax.FontSize = 10;
% colormap(linspecer)
clim([cnull.Limits]);
% cpotent.Limits = cnull.Limits;
%%
end





















