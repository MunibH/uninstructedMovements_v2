function plotME_NPMagnitude_singleTrials(meta,obj,params,me,rez,cond2use)

% plot magnitude of activity in N/Ps
% trials sorted by avg ME during late delay period


cols = linspecer(numel(rez));
% xlims = [-2 2];
% xlims = [-0.9 0];
% xlims = [-0.9 0];

% edges = [params(1).tmin, params(1).tmax];
edges = [-0.9 0];
for i = 1:numel(edges)
    [~,tix(i)] = min(abs(obj(1).time - edges(i)));
end
tix(1) = 1;
tix(2) = numel(obj(1).time);

for sessix = 1:numel(rez)
%     f = figure;
%     f.Position = [282         358        1105         382];

    % get trials
    trix = cell2mat(params(sessix).trialid(cond2use)');
    %     trix = 1:obj(sessix).bp.Ntrials;

    % get data
    tempme = me(sessix).data(:,trix);

    temprez = rez(sessix);

    potent = temprez.N_potent(:,trix,:);
    potent = mean(potent.^2,3);

    null = temprez.N_null(:,trix,:);
    null = mean(null.^2,3);

    % sort trials by avg late delay motion energy
    edges = [-0.4 -0.02];
    for i = 1:numel(edges)
        [~,t(i)] = min(abs(obj(sessix).time - edges(i)));
    end
    [~,sortix] = sort(mean(tempme(t(1):t(2),:),1),'descend');

    tempme = tempme(tix(1):tix(2),:);
    potent = potent(tix(1):tix(2),:);
    null = null(tix(1):tix(2),:);
    time = obj(sessix).time(tix(1):tix(2));

    % PLOT MOTION ENERGY
%     ax = subplot(1,3,1);
    f = figure;
    f.Position = [680   591   321   387];
    imagesc(time,1:numel(trix),tempme(:,sortix)');
    colormap(linspecer)
    c = colorbar;
    lims = clim;
    %     clim([lims(1) lims(2) / 1])
%     xlim(xlims)
    xlabel('Time (s) from go cue')
    ylabel('Trials')
%     c.Label.String = 'Motion Energy (a.u.)';

    % PLOT POTENT
    f = figure;
    f.Position = [680   591   321   387];
    potentax = gca;
    imagesc(time,1:numel(trix),potent(:,sortix)');
    colormap(linspecer);
    c = colorbar;
    lims = clim;
    %     clim([lims(1) lims(2) / 1.5])
%     xlim(xlims)
%     c.Label.String = 'Potent - Magnitude (a.u.)';

    % PLOT NULL
    f = figure;
    f.Position = [680   591   321   387];
%     ax = subplot(1,3,3);
    imagesc(time,1:numel(trix), null(:,sortix)');
    colormap(linspecer)
    c = colorbar;
    clim([0 10])
    clim_null = clim;
%     clim([clim_null(1) clim_null(2)/1.5])
%     xlim(xlims)
%     c.Label.String = 'Null - Magnitude (a.u.)';

    axes(potentax);
    clim(clim_null)

%     sgtitle([meta(sessix).anm ' ' meta(sessix).date])

    break
end

end