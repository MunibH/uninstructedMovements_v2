function plotVarianceExplained_DelayResponse(rez,meta)


% curate data

for sessix = 1:numel(rez)
    ve(1,sessix) = rez(sessix).ve.norm.null_delay;
    ve(2,sessix) = rez(sessix).ve.norm.null_resp;
    ve(3,sessix) = rez(sessix).ve.norm.potent_delay;
    ve(4,sessix) = rez(sessix).ve.norm.potent_resp;
    % ve.null_delay(sessix) = rez(sessix).ve.norm.null_delay;
    % ve.null_resp(sessix) = rez(sessix).ve.norm.null_resp;
    % ve.potent_delay(sessix) = rez(sessix).ve.norm.potent_delay;
    % ve.potent_resp(sessix) = rez(sessix).ve.norm.potent_resp;
end


% [objix,uAnm] = groupSessionsByAnimal(meta);
% nAnm = numel(uAnm);
% 
% for ianm = 1:nAnm
%     ix = find(objix{ianm});
%     anmVE.null_prep(ianm) = mean(ve.null_delay(ix));
%     anmVE.null_move(ianm) = mean(ve.null_resp(ix));
%     anmVE.potent_move(ianm) = mean(ve.potent_delay(ix));
%     anmVE.potent_prep(ianm) = mean(ve.potent_resp(ix));
% end



%% plot

defcols = getColors();

cols(1,:) = defcols.null;
cols(2,:) = defcols.potent;
cols(3,:) = defcols.null;
cols(4,:) = defcols.potent;

% fns = fieldnames(anmVE);
fns = {'null-prep','potent-prep','null-resp','potent-resp'};

f=figure;
f.Renderer = 'painters';
f.Position = [688   501   291   319];
ax = gca;
ax = prettifyPlot(ax);
hold on;
div = 1;
xs = [1 2 4 5];
for i = 1:numel(xs)
    temp = ve(i,:);
    n = numel(temp);
    h(i) = bar(xs(i),mean(temp)); % mean across dims and sessions
    cix = i;
    h(i).FaceColor = cols(cix,:);
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.5;
    %     scatter(xs(i)*ones(size(temp)),temp,10,'MarkerFaceColor','k', ...
    %         'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
    %         'MarkerFaceAlpha',1)
    e = errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);
    e.LineWidth = 0.5;
    e.CapSize = 2;

    % vs(i) = scatter(randn(n,1) * 0.1 + xs(i)*ones(n,1),temp,30,'MarkerFaceColor','k',...
    %     'MarkerEdgeColor','w','LineWidth',0.1);



    % xs_(:,i) = vs(i).XData';
    % ys_(:,i) = temp;
end

% for i = 1:size(xs_,1)
%     line(xs_(i,1:2),ys_(i,1:2),'LineWidth',0.5,'Color','k')
%     line(xs_(i,3:4),ys_(i,3:4),'LineWidth',0.5,'Color','k')
% end

% ylim([-0.001 ax.YLim(2)])
ax.XTick = xs;
xlabels  = strrep(fns,'_','-');
xticklabels(xlabels);

ylabel('Fraction of normVE')
% ax.FontSize = 12;

% title('Covariances')



end