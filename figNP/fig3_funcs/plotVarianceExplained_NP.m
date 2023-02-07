function plotVarianceExplained_NP(rez,meta)


% curate data
ve.null_total = zeros(numel(rez),1);
ve.potent_total = zeros(numel(rez),1);
ve.null_prep = zeros(numel(rez),1);
ve.null_move = zeros(numel(rez),1);
ve.potent_move = zeros(numel(rez),1);
ve.potent_prep = zeros(numel(rez),1);
for sessix = 1:numel(rez)
    ve.null_total(sessix) = rez(sessix).ve.null_total;
    ve.potent_total(sessix) = rez(sessix).ve.potent_total;
    ve.null_prep(sessix) = rez(sessix).ve.norm.null_prep;
    ve.null_move(sessix) = rez(sessix).ve.norm.null_move;
    ve.potent_move(sessix) = rez(sessix).ve.norm.potent_move;
    ve.potent_prep(sessix) = rez(sessix).ve.norm.potent_prep;
end

[objix,uAnm] = groupSessionsByAnimal(meta);
nAnm = numel(uAnm);

for ianm = 1:nAnm
    ix = find(objix{ianm});
    anmVE.null_total(ianm) = mean(ve.null_total(ix));
    anmVE.potent_total(ianm) = mean(ve.potent_total(ix));
    anmVE.null_prep(ianm) = mean(ve.null_prep(ix));
    anmVE.null_move(ianm) = mean(ve.null_move(ix));
    anmVE.potent_move(ianm) = mean(ve.potent_move(ix));
    anmVE.potent_prep(ianm) = mean(ve.potent_prep(ix));
end




%% plot

defcols = getColors();

% cols(1,:) = defcols.null * 255 / 2;
% cols(2,:) = defcols.potent * 255 / 2;
% cols(3,:) = defcols.null * 255;
% cols(4,:) = defcols.potent * 255;
% cols(5,:) = cols(4,:);
% cols(6,:) = cols(3,:);
% cols = cols ./ 255;

cols(1,:) = defcols.null * 255;
cols(2,:) = defcols.potent * 255;
cols(3,:) = cols(2,:);
cols(4,:) = cols(1,:);
cols = cols ./ 255;

% cols = linspecer(10);

fns = fieldnames(anmVE);
fns = fns(3:end);

f=figure; 
f.Position = [688   501   291   319];
ax = gca;
hold on;
div = 1;
% xs = [1 2 4 5 7 8];
xs = [1 2 4 5];
for i = 1:numel(fns)
    temp = anmVE.(fns{i});
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

    vs(i) = scatter(randn(nAnm,1) * 0.1 + xs(i)*ones(nAnm,1),temp,10,'MarkerFaceColor','k',...
        'MarkerEdgeColor','k','LineWidth',1);



    xs_(:,i) = vs(i).XData';
    ys_(:,i) = temp;
end

for i = 1:size(xs_,1)
    patchline(xs_(i,1:2),ys_(i,1:2),'EdgeAlpha',0.4,'LineWidth',0.1)
    patchline(xs_(i,3:4),ys_(i,3:4),'EdgeAlpha',0.4,'LineWidth',0.1)
end

% ylim([-0.001 ax.YLim(2)])
ax.XTick = xs;
xlabels  = strrep(fns,'_','-');
xticklabels(xlabels);

ylabel('Fraction of normVE')
ax.FontSize = 12;

title('Covariances')

end