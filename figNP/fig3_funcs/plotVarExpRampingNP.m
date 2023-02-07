function plotVarExpRampingNP(cd_null,cd_potent,obj,rez,params)

rampix = find(strcmpi(cd_null(1).cd_labels,'ramping'));

[ve.null,r2.null] = varExpRampingNP(obj,rez,cd_null,params,rampix,'null');
[ve.potent,r2.potent] = varExpRampingNP(obj,rez,cd_potent,params,rampix,'potent');
% varExpRampingNP(obj,rez,cd_potent_all,params)

cols = getColors();

dat = r2;

fns = fieldnames(dat);

f=figure; hold on;
ax = gca;
xs = [1 2];
for i = 1:numel(fns)
    temp = dat.(fns{i});
    h(i) = bar(xs(i),mean(temp)); % mean across dims and sessions
    cix = i;
    h(i).FaceColor = cols.(fns{i});
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.5;
    scatter(xs(i)*ones(size(temp)),temp,60,'MarkerFaceColor',[0 0 0], ...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
        'MarkerFaceAlpha',0.7)
    errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);
end
ax.XTick = xs;
xlabels  = strrep(fns,'_','-');
xticklabels(xlabels);

ylabel('R2 b/w ramping mode recon and PSTH')
ax.FontSize = 12;