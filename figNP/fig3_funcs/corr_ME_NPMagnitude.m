function [corr_me_null, corr_me_potent] = corr_ME_NPMagnitude(meta,obj,params,me,rez,cond2use)

% correlate magnitude of activity in N/Ps and ME
% plot by mouse and session

ms = {'o','<','^','v','>','square','diamond','o','<'};
sz = 50;

[objix,uAnm]  = groupSessionsByAnimal(meta);
nAnm = numel(uAnm);

f = figure;
f.Position = [680   692   299   186];
ax = gca;
hold on;

% edges = [params(1).tmin, params(1).tmax];
edges = [-0.9 0];
for i = 1:numel(edges)
    [~,tix(i)] = min(abs(obj(1).time - edges(i)));
end

%%


% for each session, calculate avg feat value over trials
for isess = 1:numel(meta)
    sessix = isess;

    % get data
    tempme = me(sessix).data;
    tempme = tempme(tix(1):tix(2),:);
    tempme = tempme(:);
    tempme = fillmissing(tempme,'nearest');

    temprez = rez(sessix);

    potent = temprez.N_potent;
    potent = sum(potent.^2,3);
    potent = (potent(tix(1):tix(2),:));
    potent = potent(:);

    null = temprez.N_null;
    null = sum(null.^2,3);
    null = (null(tix(1):tix(2),:));
    null = null(:);

    % correlate
    corr_me_potent(isess) = corr(tempme,potent)^2;
    corr_me_null(isess) = corr(tempme,null)^2;

   
end


s = scatter(corr_me_null,corr_me_potent,sz, ...
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
plot(ax.XLim,ax.YLim,'k--','LineWidth',2)
% view([90 -90])






%%
%
% for ianm = 1:nAnm
%     % get sessions for current animal
%     ix = find(objix{ianm});
%
%     % for each session, calculate avg feat value over trials
%     for isess = 1:numel(ix)
%         sessix = ix(isess);
%
%         % get data
%         tempme = me(sessix).data;
%         tempme = tempme(tix(1):tix(2),:);
%         tempme = tempme(:);
%         tempme = fillmissing(tempme,'nearest');
%
%         temprez = rez(sessix);
%
%         potent = temprez.N_potent;
%         potent = sum(potent.^2,3);
%         potent = (potent(tix(1):tix(2),:));
%         potent = potent(:);
%
%         null = temprez.N_null;
%         null = sum(null.^2,3);
%         null = (null(tix(1):tix(2),:));
%         null = null(:);
%
%         % correlate
%         corr_me_potent{ianm}(isess) = corr(tempme,potent)^2;
%         corr_me_null{ianm}(isess) = corr(tempme,null)^2;
%
%     end
%
%     mu.potent = nanmean(corr_me_potent{ianm});
%     mu.null = nanmean(corr_me_null{ianm});
%     if ianm > 7
%         s = scatter(mu.null,mu.potent,sz, ...
%             'MarkerEdgeColor','k','MarkerFaceColor','none','Marker',ms{ianm});
%     else
%         s = scatter(mu.null,mu.potent,sz, ...
%             'MarkerEdgeColor','w','MarkerFaceColor','k','Marker',ms{ianm});
%     end
% end
%
% xlabel('R2(ME,Null)', 'Interpreter','none')
% ylabel('R2(ME,Potent)', 'Interpreter','none')
% ax.FontSize = 9;
% axis(ax,'equal')
%
% mins = min([ax.XLim(1) ax.YLim(1)]);
% maxs = max([ax.XLim(2) ax.YLim(2)]);
%
% xlim([0 maxs])
% ylim([0 maxs])
% ax = gca;
% plot(ax.XLim,ax.YLim,'k--','LineWidth',2)
% % view([90 -90])
%
%
% h = zeros(nAnm, 1);
% for ianm = 1:numel(h)
%     if ianm > 7
%         h(ianm) = scatter(NaN,NaN,sz,'MarkerEdgeColor','k','MarkerFaceColor','none','Marker',ms{ianm});
%     else
%         h(ianm) = scatter(NaN,NaN,sz,'MarkerEdgeColor','none','MarkerFaceColor','k','Marker',ms{ianm});
%     end
%
% end
% legString = uAnm;
%
% leg = legend(h, legString);
% leg.EdgeColor = 'none';
% leg.Color = 'none';
% leg.Location = 'best';



