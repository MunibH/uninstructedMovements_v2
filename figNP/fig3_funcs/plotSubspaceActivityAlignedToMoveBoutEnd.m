function plotSubspaceActivityAlignedToMoveBoutEnd(dat,obj,me,rez,params,meta)

rng(pi)


% m2q
for sessix = 1:numel(me)
    f = figure;
    %     f.Position = [403   412   809   543];
    ax = gca;
    hold on

    dt = params(sessix).dt;

    % transition data
    mdat = dat.mdat{sessix};
    nTrials = numel(mdat);

    nbouts = 5; % approximately max number of bouts in a trial
    allnull = nan(numel(obj(sessix).time),nTrials*nbouts);
    allpotent = nan(numel(obj(sessix).time),nTrials*nbouts);
    allme = nan(numel(obj(sessix).time),nTrials*nbouts);

    temprez = rez(sessix);

    for trix = 1:nTrials
        mdat_trix = mdat{trix};
        if isempty(mdat_trix)
            continue
        end

        null = sum(rez(sessix).N_null(:,trix,:).^2,3); % (time,1)
        potent = sum(rez(sessix).N_potent(:,trix,:).^2,3); % (time,1)
        medat = me(sessix).data(:,trix);

        % for each transition
        for i = 1:size(mdat_trix,2)
            tix = mdat_trix(1,i):mdat_trix(end,i);
            ts = tix .* dt;
            ts = ts - mean(ts);

            nullplot = null(tix);
            potentplot = potent(tix);
            meplot = medat(tix);

            % get center of tix
            center = ceil(numel(tix)/2);
            try
                allnull(250-center:250+center-1,(trix+i-1)) = nullplot;
                allpotent(250-center:250+center-1,(trix+i)-1) = potentplot;
                allme(250-center:250+center-1,(trix+i)-1) = meplot;
            catch
                allnull(250-center:250+center-2,(trix+i-1)) = nullplot;
                allpotent(250-center:250+center-2,(trix+i)-1) = potentplot;
                allme(250-center:250+center-2,(trix+i)-1) = meplot;
            end
        end
    end

    % sumsqmag
    %     meannull = sum(allnull.^2,3,'omitnan');
    %     meanpotent = sum(allpotent.^2,3,'omitnan');
    meannull = allnull;
    meanpotent = allpotent;


    nullmeanplot = nanmean(meannull,2);
    nullerrplot = nanstd(meannull,[],2) ./ sqrt(size(allnull,2));
    potentmeanplot = nanmean(meanpotent,2);
    potenterrplot = nanstd(meanpotent,[],2) ./ sqrt(size(allnull,2));

    yyaxis(ax,'left')
    shadedErrorBar(obj(sessix).time,nullmeanplot,nullerrplot,{'Color',[62, 168, 105]./255,'LineWidth',2,'LineStyle','-'},0.5,ax)
    shadedErrorBar(obj(sessix).time,potentmeanplot,potenterrplot,{'Color',[255, 56, 140]./255,'LineWidth',2,'LineStyle','-'},0.5,ax)
    ax.YColor = [100 100 100 ]./255;
    xlabel('Time to quiet (s)')
    ylabel('SSM - all dims (a.u)')

    ylims = ax.YLim;

    %     meanme = squeeze(normalize(nanmean(allme,2),'range',[ylims(1)+0.1 ylims(2)-0.1]));
    meanme = squeeze(nanmean(allme,2));
    %     meanme = squeeze(nanmean(allme,2));
    meanme = meanme(:,1);
    yyaxis(ax,'right')
    plot(obj(sessix).time,meanme,'Color',[225, 144, 15]./255,'LineWidth',3);
    ax.YColor = [225, 144, 15]./255;
    ylabel('Motion Energy')

    xline(0,'k:','LineWidth',2)

    xlim([-0.5 0.5])
    title([meta(sessix).anm ' - ' meta(sessix).date])


    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'-','Color',[225, 144, 15]./255);
    h(2) = plot(NaN,NaN,'-','Color',[62, 168, 105]./255);
    h(3) = plot(NaN,NaN,'-','Color',[255, 56, 140]./255);
    legend(h, 'motion energy','null','potent');

    break

end




end