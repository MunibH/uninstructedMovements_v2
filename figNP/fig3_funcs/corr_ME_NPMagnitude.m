function [corr_me_null, corr_me_potent] = corr_ME_NPMagnitude(meta,obj,params,me,rez,cond2use)

% correlate magnitude of activity in N/Ps and ME
% plot by mouse and session

ms = {'o','<','^','v','>','square','diamond','o'};
sz = 50;

allanms = {meta(:).anm}; uAnm = unique(allanms); nAnm = numel(uAnm);
alldates = {meta(:).date};

f = figure;
ax = gca;
hold on;
for ianm = 1:nAnm
    % get sessions for current animal
    ix = ismember(allanms,uAnm{ianm});
    sessionDates = {alldates{ix}};

    % for each session, calculate avg feat value over trials
    for isess = 1:numel(sessionDates)
        sessix = ismember(allanms,uAnm{ianm}) & ismember(alldates,sessionDates{isess});

        % get data
        tempme = me(sessix).data;
        tempme = tempme(:);
        tempme = fillmissing(tempme,'nearest');

        temprez = rez(sessix);

        potent = temprez.N_potent;
        potent = sum(potent.^2,3);
        potent = potent(:);

        null = temprez.N_null;
        null = sum(null.^2,3);
        null = null(:);

        % correlate
        corr_me_potent{ianm}(isess) = corr(tempme,potent)^2;
        corr_me_null{ianm}(isess) = corr(tempme,null)^2;

        if ianm == 8
            s = scatter(corr_me_null{ianm}(isess),corr_me_potent{ianm}(isess),sz, ...
                'MarkerEdgeColor','k','MarkerFaceColor','none','Marker',ms{ianm});
        else
            s = scatter(corr_me_null{ianm}(isess),corr_me_potent{ianm}(isess),sz, ...
                'MarkerEdgeColor','w','MarkerFaceColor','k','Marker',ms{ianm});
        end

    end
end

xlabel('R2(ME,Null)', 'Interpreter','none')
ylabel('R2(ME,Potent)', 'Interpreter','none')
ax.FontSize = 9;
axis(ax,'equal')

mins = min([ax.XLim(1) ax.YLim(1)]);
maxs = max([ax.XLim(2) ax.YLim(2)]);

xlim([0 maxs])
ylim([0 maxs])
ax = gca;
plot(ax.XLim,ax.YLim,'k--','LineWidth',2)
% view([90 -90])


h = zeros(nAnm, 1);
for ianm = 1:numel(h)
    if ianm == 8
        h(ianm) = scatter(NaN,NaN,sz,'MarkerEdgeColor','k','MarkerFaceColor','none','Marker',ms{ianm});
    else
        h(ianm) = scatter(NaN,NaN,sz,'MarkerEdgeColor','none','MarkerFaceColor','k','Marker',ms{ianm});
    end

end
legString = uAnm;

leg = legend(h, legString);
leg.EdgeColor = 'none';
leg.Color = 'none';
leg.Location = 'best';



