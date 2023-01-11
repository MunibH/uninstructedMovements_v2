function plotME_NPMagnitude_singleTrials(meta,obj,params,me,rez,cond2use)

% plot magnitude of activity in N/Ps
% trials sorted by avg ME during late delay period


cols = linspecer(numel(rez));
xlims = [-2 2];


for sessix = 1:numel(rez)
    f = figure;
    f.Position = [282         397        1231         343];

    % get trials
    trix = cell2mat(params(sessix).trialid(cond2use)');
    %     trix = 1:obj(sessix).bp.Ntrials;

    % get data
    tempme = me(sessix).data(:,trix);

    temprez = rez(sessix);

    potent = temprez.N_potent(:,trix,:);
    potent = sum(potent.^2,3);

    null = temprez.N_null(:,trix,:);
    null = sum(null.^2,3);

    % sort trials by avg late delay motion energy
    edges = [-0.4 -0.02];
    for i = 1:numel(edges)
        [~,t(i)] = min(abs(obj(sessix).time - edges(i)));
    end
    [~,sortix] = sort(mean(tempme(t(1):t(2),:),1),'descend');

    % PLOT MOTION ENERGY
    ax = subplot(1,3,1);
    imagesc(obj(sessix).time,1:numel(trix),tempme(:,sortix)');
    colormap(linspecer)
    c = colorbar;
    lims = clim;
    %     clim([lims(1) lims(2) / 1])
    xlim(xlims)
    xlabel('Time (s) from go cue')
    ylabel('Trials')
    c.Label.String = 'Motion Energy (a.u.)';

    % PLOT POTENT
    potentax = subplot(1,3,2);
    imagesc(obj(sessix).time,1:numel(trix),potent(:,sortix)');
    colormap(linspecer);
    c = colorbar;
    lims = clim;
    %     clim([lims(1) lims(2) / 1.5])
    xlim(xlims)
    c.Label.String = 'Potent - Magnitude (a.u.)';

    % PLOT NULL
    ax = subplot(1,3,3);
    imagesc(obj(sessix).time,1:numel(trix), null(:,sortix)');
    colormap(linspecer)
    c = colorbar;
    clim_null = clim;
    %     clim([lims(1) lims(2) / 1.5])
    xlim(xlims)
    c.Label.String = 'Null - Magnitude (a.u.)';

    axes(potentax);
    clim(clim_null)

    sgtitle([meta(sessix).anm ' ' meta(sessix).date])


end

end