close all

edges = [-0.4 -0.02]; % relative to to go cue
for i = 1:numel(edges)
    [~,tix(i)] = min(abs(obj(1).time - edges(i)));
end

sm = 11;

for sessix = 1:numel(rez)
    temprez = rez(sessix);
    trix = cell2mat(params(sessix).trialid(2:3)');

    tempme = me(sessix).data(:,trix);
    
    null = rez(sessix).null_ssm;
    potent = rez(sessix).potent_ssm;
   
    [~,ix] = sort(mean(tempme(tix(1):tix(2),:),1),'descend');
    trix = ix;

    f = figure;
    f.Position = [510   453   831   336];

    ax1 = subplot(1,3,1);
    plotme = tempme(:,trix);
    imagesc(obj(sessix).time,1:numel(trix),plotme');
    colormap(linspecer)
    lims = clim;
%     clim([lims(1) lims(2) / 1])
    colorbar;
    ylabel('Trials')

    ax2 = subplot(1,3,2);
    potent = mySmooth(normalize(potent(:,trix)),sm);
    imagesc(obj(sessix).time,1:numel(trix),potent');
    colormap(linspecer);
    lims = clim;
%     clim([lims(1) lims(2) / 1.5])
    colorbar;
    xlabel('Time from go cue (s)')

    ax3 = subplot(1,3,3);
    null = mySmooth(normalize(null(:,trix)),sm);
    imagesc(obj(sessix).time,1:numel(trix), null');
    colormap(linspecer)
    lims = clim;
%     clim([lims(1) lims(2) / 1.5])
    colorbar;

    sgtitle([meta(sessix).anm ' ' meta(sessix).date])
    pause
    close all

end