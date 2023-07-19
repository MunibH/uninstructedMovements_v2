%% choice
close all


sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;
gc = 0;

nullsm = 11;
potentsm = 0;

% only show delay period
delayedges = [-2.39 0];
delayt = findTimeIX(obj(1).time,delayedges);
delayt = delayt(1):delayt(2);

edges = [-0.4 -0.02]; % motion energy sort times
t = findTimeIX(obj(1).time,edges);

cdix = 1;
cond2use = [8 9];
for sessix = 3%1:numel(meta)  %[2 3 ]
    trix{1} = params(sessix).trialid{cond2use(1)};
    trix{2} = params(sessix).trialid{cond2use(2)};
    % trix = cell2mat(params(sessix).trialid(cond2use)');

    dat = rez(sessix).recon.null;
    W = cd_null(sessix).cd_mode_orth(:,cdix);
    proj.null = tensorprod(dat,W,3,1);

    dat = rez(sessix).recon.potent;
    W = cd_potent(sessix).cd_mode_orth(:,cdix);
    proj.potent = tensorprod(dat,W,3,1);

    for j = 1:numel(trix)
        tempme{j} = me(sessix).data(:,trix{j});
        nullts{j} = proj.null(:,trix{j});
        potentts{j} = proj.potent(:,trix{j});

        % sort trials by avg late delay motion energy
        [~,sortix] = sort(mean(tempme{j}(t(1):t(2),:),1),'descend');

        tempme{j} = tempme{j}(delayt,sortix);
        nullts{j} = nullts{j}(delayt,sortix);
        potentts{j} = potentts{j}(delayt,sortix);

    end
    
    plt.me = cat(2,tempme{1},tempme{2});
    plt.null = cat(2,nullts{1},nullts{2});
    plt.potent = cat(2,potentts{1},potentts{2});

    plt.null = mySmooth(plt.null,nullsm,'reflect');
    plt.potent = mySmooth(plt.potent,potentsm,'reflect');

    % plt.null = zscore(plt.null);
    % plt.potent = zscore(plt.potent);

    
    f = figure;
    f.Position = [412   417   929   336];
    f.Renderer = 'painters';

    nTrials = size(plt.me,2);
    nCond = size(tempme{1},2);
    time_ = obj(1).time(delayt);

    ax = subplot(1,3,1);
    imagesc(time_,1:size(plt.me,2),plt.me'); colorbar; %clim([-5 10])
    ax = prettifyPlot(ax);
    colormap(ax,parula);
    line([sample sample], [1 nTrials],'Color','w','LineStyle','--');
    line([delay delay], [1 nTrials],'Color','w','LineStyle','--');
    % line([gc gc], [1 nTrials],'Color','w','LineStyle','--');
    line([time_(1) time_(end)], [nCond nCond],'Color','w','LineStyle','--');

    ax = subplot(1,3,2);
    imagesc(time_,1:size(plt.me,2),plt.potent'); colorbar; %clim([-5 5])
    ax = prettifyPlot(ax);
    colormap(ax,parula);
    line([sample sample], [1 nTrials],'Color','w','LineStyle','--');
    line([delay delay], [1 nTrials],'Color','w','LineStyle','--');
    % line([gc gc], [1 nTrials],'Color','w','LineStyle','--');
    line([time_(1) time_(end)], [nCond nCond],'Color','w','LineStyle','--');

    ax = subplot(1,3,3);
    % plt.null(:,80:end) = plt.null(:,80:end)-0.5;
    imagesc(time_,1:size(plt.me,2),plt.null'); colorbar; clim([-3 2])
    ax = prettifyPlot(ax);
    colormap(ax,parula);
    line([sample sample], [1 nTrials],'Color','w','LineStyle','--');
    line([delay delay], [1 nTrials],'Color','w','LineStyle','--');
    % line([gc gc], [1 nTrials],'Color','w','LineStyle','--');
    line([time_(1) time_(end)], [nCond nCond],'Color','w','LineStyle','--');

    sgtitle([meta(sessix).anm ' ' meta(sessix).date])


    cc = corrcoef(plt.me(:),plt.null(:));
    r2.null(sessix) = cc(1,2).^2;
    cc = corrcoef(plt.me(:),plt.potent(:));
    r2.potent(sessix) = cc(1,2).^2;

    % pth = 'C:\Users\munib\Desktop\NP_CDChoice_Heatmaps';
    % fn = [meta(sessix).anm ' ' meta(sessix).date '.svg'];
    % saveas(f,fullfile(pth,fn),'svg')

    break

end

sz = 40;

f = figure;
f.Position = [680   692   299   186];
f.Renderer = 'painters';
ax = gca;
ax = prettifyPlot(ax);
hold on;

s = scatter(r2.null,r2.potent,sz, ...
    'MarkerEdgeColor','w','MarkerFaceColor','k');

xlabel('R2(ME,Null)', 'Interpreter','none')
ylabel('R2(ME,Potent)', 'Interpreter','none')
ax.FontSize = 9;
axis(ax,'equal')

mins = min([ax.XLim(1) ax.YLim(1)]);
maxs = max([ax.XLim(2) ax.YLim(2)]);

xlim([0 maxs])
ylim([0 maxs])
ax = gca;
setLimits(ax,0.1);
plot(ax.XLim,ax.YLim,'k--','LineWidth',2)

%% ramp

close all


sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;
gc = 0;

nullsm = 11;
potentsm = 0;

% only show delay period
delayedges = [-2.39 0];
delayt = findTimeIX(obj(1).time,delayedges);
delayt = delayt(1):delayt(2);

edges = [-0.4 -0.02]; % motion energy sort times
t = findTimeIX(obj(1).time,edges);

cdix = 3;
cond2use = [8 9];
for sessix = 1:numel(meta)  %[2 3 ]
    trix = cell2mat(params(sessix).trialid(cond2use)');

    tempme = me(sessix).data(:,trix);

    dat = rez(sessix).recon.null;
    W = cd_null(sessix).cd_mode_orth(:,cdix);
    proj.null = tensorprod(dat,W,3,1);
    proj.null = proj.null(:,trix);

    dat = rez(sessix).recon.potent;
    W = cd_potent(sessix).cd_mode_orth(:,cdix);
    proj.potent = tensorprod(dat,W,3,1);
    proj.potent = proj.potent(:,trix);

    % sort trials by avg late delay motion energy
    [~,sortix] = sort(mean(tempme(t(1):t(2),:),1),'descend');

    plt.me = tempme(delayt,sortix);
    plt.null = proj.null(delayt,sortix);
    plt.potent = proj.potent(delayt,sortix);

    plt.null = mySmooth(plt.null,nullsm,'reflect');
    plt.potent = mySmooth(plt.potent,potentsm,'reflect');
    
    f = figure;
    f.Position = [412   417   929   336];
    f.Renderer = 'painters';

    nTrials = size(plt.me,2);
    time_ = obj(1).time(delayt);

    ax = subplot(1,3,1);
    imagesc(time_,1:size(plt.me,2),plt.me'); colorbar; %clim([-5 10])
    ax = prettifyPlot(ax);
    colormap(ax,parula);
    line([sample sample], [1 nTrials],'Color','w','LineStyle','--');
    line([delay delay], [1 nTrials],'Color','w','LineStyle','--');
    % line([gc gc], [1 nTrials],'Color','w','LineStyle','--');
    % line([time_(1) time_(end)], [nCond nCond],'Color','w','LineStyle','--');

    ax = subplot(1,3,2);
    imagesc(time_,1:size(plt.me,2),plt.potent'); colorbar; %clim([-5 5])
    ax = prettifyPlot(ax);
    colormap(ax,parula);
    line([sample sample], [1 nTrials],'Color','w','LineStyle','--');
    line([delay delay], [1 nTrials],'Color','w','LineStyle','--');
    % line([gc gc], [1 nTrials],'Color','w','LineStyle','--');
    % line([time_(1) time_(end)], [nCond nCond],'Color','w','LineStyle','--');

    ax = subplot(1,3,3);
    imagesc(time_,1:size(plt.me,2),plt.null'); colorbar; %clim([-5 10])
    ax = prettifyPlot(ax);
    colormap(ax,parula);
    line([sample sample], [1 nTrials],'Color','w','LineStyle','--');
    line([delay delay], [1 nTrials],'Color','w','LineStyle','--');
    % line([gc gc], [1 nTrials],'Color','w','LineStyle','--');
    % line([time_(1) time_(end)], [nCond nCond],'Color','w','LineStyle','--');

    sgtitle([meta(sessix).anm ' ' meta(sessix).date])


    cc = corrcoef(plt.me(:),plt.null(:));
    r2.null(sessix) = cc(1,2).^2;
    cc = corrcoef(plt.me(:),plt.potent(:));
    r2.potent(sessix) = cc(1,2).^2;

    pth = 'C:\Users\munib\Desktop\NP_CDRamp_Heatmaps';
    fn = [meta(sessix).anm ' ' meta(sessix).date '.svg'];
    saveas(f,fullfile(pth,fn),'svg')

    % break

end

sz = 40;

f = figure;
f.Position = [680   692   299   186];
f.Renderer = 'painters';
ax = gca;
ax = prettifyPlot(ax);
hold on;

s = scatter(r2.null,r2.potent,sz, ...
    'MarkerEdgeColor','w','MarkerFaceColor','k');

xlabel('R2(ME,Null)', 'Interpreter','none')
ylabel('R2(ME,Potent)', 'Interpreter','none')
ax.FontSize = 9;
axis(ax,'equal')

mins = min([ax.XLim(1) ax.YLim(1)]);
maxs = max([ax.XLim(2) ax.YLim(2)]);

xlim([0 maxs])
ylim([0 maxs])
ax = gca;
setLimits(ax,0.1);
plot(ax.XLim,ax.YLim,'k--','LineWidth',2)

