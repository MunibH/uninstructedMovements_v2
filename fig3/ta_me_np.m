for sessix = 1:numel(meta)

    tempme = me(sessix);
    temprez = rez(sessix);
    tempobj = obj(sessix);
    tempparams = params(sessix);

    %

    % trix = randsample(tempobj.bp.Ntrials,1,false); % 95
    trix{1} = tempparams.trialid{2};
    trix{2} = tempparams.trialid{3};
    trix{3} = tempparams.trialid{4};
    trix{4} = tempparams.trialid{5};
    for i = 1:numel(trix)
        plotme{i} = tempme.data(:,trix{i});
        plotme{i} = normalize(mean(plotme{i},2),'range',[0 1]);

        potent{i} = squeeze(temprez.N_potent(:,trix{i},:));
        potent{i} = squeeze(sum(potent{i}.^2,2));
        potent{i} = normalize(mean(potent{i},2),'range',[0 1]);

        null{i} = squeeze(temprez.N_null(:,trix{i},:));
        null{i} = squeeze(sum(null{i}.^2,2));
        null{i} = normalize(mean(null{i},2),'range',[0 1]);
    end

    %
    dy = 1;
    cols = getColors;
    fns = fieldnames(cols);
    for i = 1:numel(fns)
        clrs{i} = cols.(fns{i});
    end
    div = 1.3;
    f = figure;
    ax = gca;
    hold on
    for i = 1:numel(trix)
        plot(tempobj.time,mySmooth(plotme{i},7),'Color',clrs{i},'LineWidth',2)
        plot(tempobj.time,mySmooth(potent{i},7)+dy,'Color',clrs{i},'LineWidth',2)
        plot(tempobj.time,mySmooth(null{i},7)+dy*2,'Color',clrs{i},'LineWidth',2)
    end

    xlim([tempobj.time(5), 2.5])

    xline(0,'k--')


end









