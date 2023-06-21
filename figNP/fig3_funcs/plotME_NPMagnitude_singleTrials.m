function [nmag,pmag] = plotME_NPMagnitude_singleTrials(meta,obj,params,me,rez,cond2use)

% plot magnitude of activity in N/Ps
% trials sorted by avg ME during late delay period


cols = linspecer(numel(rez));
defcols = getColors;

% xlims = [-2.4 2];

xlims = [-0.9 0];
% xlims = [params(1).tmin, params(1).tmax];


edges = xlims;
for i = 1:numel(edges)
    [~,tix(i)] = min(abs(obj(1).time - edges(i)));
end
tix(1) = 1;
tix(2) = numel(obj(1).time);

sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;
gc = 0;

sm = 0;

for sessix = 3%1:numel(rez) % 15
   f = figure;
    f.Renderer = 'painters';
    f.Position = [194   499   904   273];

    % get trials
    trix = cell2mat(params(sessix).trialid(cond2use)');
    %     trix = 1:obj(sessix).bp.Ntrials;

    % get data
    tempme = me(sessix).data(:,trix);

    temprez = rez(sessix);

    potent = temprez.N_potent(:,trix,:);
    potent = (mean(potent.^2,3));

    null = temprez.N_null(:,trix,:);
    null = (mean(null.^2,3));
    gct = [0 2];
    gcix = findTimeIX(obj(1).time,gct);
    gcix = gcix(1):gcix(2);
    null(gcix,:) = null(gcix,:)/1.5;

    % sort trials by avg late delay motion energy
    edges = [-0.4 -0.02];
    for i = 1:numel(edges)
        [~,t(i)] = min(abs(obj(sessix).time - edges(i)));
    end
    [~,sortix] = sort(mean(tempme(t(1):t(2),:),1),'descend');

    
    tempme = tempme(tix(1):tix(2),:);
    tempme = mySmooth(tempme,11,'reflect');
    potent = potent(tix(1):tix(2),:);
    % potent = normalize(potent(tix(1):tix(2),:));
    potent = mySmooth(potent,sm);
    null = null(tix(1):tix(2),:);
    % null = normalize(null(tix(1):tix(2),:));
    null = mySmooth(null,sm);
    time = obj(sessix).time(tix(1):tix(2));

    % PLOT MOTION ENERGY
    ax1 = subplot(1,3,1);
    imagesc(time,1:numel(trix),tempme(:,sortix)');
    colormap(ax1,parula)
    c = colorbar;
    line([sample sample], [1 numel(trix)],'Color','w','LineStyle','--');
    line([delay delay], [1 numel(trix)],'Color','w','LineStyle','--');
    line([gc gc], [1 numel(trix)],'Color','w','LineStyle','--');
    clim([0 80]);
    xlim(xlims)
    ax1 = prettifyPlot(ax1);
%     xlabel('Time (s) from go cue')
%     ylabel('Trials')
%     c.Label.String = 'Motion Energy (a.u.)';

    % PLOT POTENT
    ax2 = subplot(1,3,2);
    pmag = potent(:,sortix);
    imagesc(time,1:numel(trix),potent(:,sortix)');
    colormap(ax2,linspecer);
    % cm = colorbarpwn(0, 10, 'level', 20, 'colorN',defcols.potent, 'colorP', defcols.potent);
    c = colorbar;
    line([sample sample], [1 numel(trix)],'Color','w','LineStyle','--');
    line([delay delay], [1 numel(trix)],'Color','w','LineStyle','--');
    line([gc gc], [1 numel(trix)],'Color','w','LineStyle','--');
    lims = clim;
    % clim([lims(1)+0.5 lims(2)/4])
    xlim(xlims)
    ax2 = prettifyPlot(ax2);
%     c.Label.String = 'Potent - Magnitude (a.u.)';

    % PLOT NULL
    ax3 = subplot(1,3,3);
    nmag = null(:,sortix);
    imagesc(time,1:numel(trix), null(:,sortix)');
    colormap(ax3,linspecer)
    c = colorbar;
    line([sample sample], [1 numel(trix)],'Color','w','LineStyle','--');
    line([delay delay], [1 numel(trix)],'Color','w','LineStyle','--');
    line([gc gc], [1 numel(trix)],'Color','w','LineStyle','--');
    lims = clim;
    clim([lims(1) lims(2)/1.5])
    xlim(xlims)
%     c.Label.String = 'Null - Magnitude (a.u.)';
    ax3 = prettifyPlot(ax3);

    sgtitle([meta(sessix).anm ' ' meta(sessix).date])

    % break
end

end