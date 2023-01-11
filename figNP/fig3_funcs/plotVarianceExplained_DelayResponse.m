function plotVarianceExplained_DelayResponse(rez)


% curate data
ve.null_delay = zeros(numel(rez),1);
ve.null_resp = zeros(numel(rez),1);
ve.potent_delay = zeros(numel(rez),1);
ve.potent_resp = zeros(numel(rez),1);
for sessix = 1:numel(rez)
    ve.null_delay(sessix) = rez(sessix).ve.norm.null_delay;
    ve.null_resp(sessix) = rez(sessix).ve.norm.null_resp;
    ve.potent_delay(sessix) = rez(sessix).ve.norm.potent_delay;
    ve.potent_resp(sessix) = rez(sessix).ve.norm.potent_resp;
end


%% plot

defcols = getColors();

cols(1,:) = defcols.null;
cols(2,:) = defcols.potent;
cols(3,:) = defcols.null;
cols(4,:) = defcols.potent;

fns = fieldnames(ve);

f=figure; hold on;
ax = gca;
div = 1;
xs = [1 2 4 5];
for i = 1:numel(fns)
    temp = ve.(fns{i});
    h(i) = bar(xs(i),mean(temp)); % mean across dims and sessions
    cix = i;
    h(i).FaceColor = cols(cix,:);
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.5;
    scatter(xs(i)*ones(size(temp)),temp,40,'MarkerFaceColor',cols(cix,:)./div, ...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
        'MarkerFaceAlpha',0.7)
    errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);
end
% ylim([-0.001 ax.YLim(2)])
ax.XTick = xs;
xlabels  = strrep(fns,'_','-');
xticklabels(xlabels);

ylabel('Fraction of VE')
ax.FontSize = 10;



end