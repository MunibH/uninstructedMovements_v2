
temp = load('nmag_dr_wc.mat');
n.dr_wc = temp.nmag;
temp = load('nmag_wc.mat');
n.wc = temp.nmag;
temp = load('nmag_dr_wc_response.mat');
n.dr_wc_resp = temp.nmag;

temp = load('pmag_dr_wc.mat');
p.dr_wc = temp.pmag;
temp = load('pmag_wc.mat');
p.wc = temp.pmag;
temp = load('pmag_dr_wc_response.mat');
p.dr_wc_resp = temp.pmag;


%%

times = [-0.9 0];
ix = findTimeIX(obj(1).time,times);
ix = ix(1):ix(2);

%% null dr_wc and null wc resp
close all
cols = getColors;

f = figure;
f.Position = [680   678   233   200];
f.Renderer = 'painters';
ax = gca;
ax = prettifyPlot(ax);
hold on;
temp1 = mean(n.dr_wc(ix,:),1);
% temp1 = temp1(:);
temp2 = mean(n.dr_wc_resp(ix,:),1);
% temp2 = temp2(:);
mdl = fitlm(temp1,temp2);
pp = plot(mdl);
title(['R^2 = ' num2str(mdl.Rsquared.Ordinary)])
pp(1).Marker = 'o';
pp(1).MarkerSize = 5;
pp(1).MarkerFaceColor = 'none';
pp(1).MarkerEdgeColor = cols.null;
pp(2).Color = 'k';
pp(2).LineStyle = '--';

xs = pp(2).XData;
y1 = pp(2).YData;
y2 = pp(3).YData;
% ax = gca; hold on;
% shadedErrorBar(xs,y1,y2-y1,{'Color','k'},0.2,ax)
pp(3).Visible ='off';
pp(4).Visible = 'off';

xlabel('Null DR WC')
ylabel('Null DR WC Resp')
% ax.FontSize = 11;
ax.Legend.Visible = 'off';


%% null dr_wc and null dr_wc_resp

close all
cols = getColors;

f = figure;
f.Position = [680   678   233   200];
f.Renderer = 'painters';
ax = gca;
ax  = prettifyPlot(ax);
hold on;
temp1 = mean(n.dr_wc(ix,:),1);
% temp1 = temp1(:);
temp2 = mean(n.wc(ix,:),1);
% temp2 = temp2(:);
mdl = fitlm(temp1,temp2);
pp = plot(mdl);
title(['R^2 = ' num2str(mdl.Rsquared.Ordinary)])
pp(1).Marker = 'o';
pp(1).MarkerSize = 7;
pp(1).MarkerEdgeColor = cols.null;
% pp(1). = 'w';
pp(2).Color = 'k';
pp(2).LineStyle = '--';

xs = pp(2).XData;
y1 = pp(2).YData;
y2 = pp(3).YData;
% ax = gca; hold on;
% shadedErrorBar(xs,y1,y2-y1,{'Color','k'},0.2,ax)
pp(3).Visible ='off';
pp(4).Visible = 'off';

xlabel('Null DR WC')
ylabel('Null WC')
% ax.FontSize = 11;
ax.Legend.Visible = 'off';


%% potent dr_wc and potent wc
close all
cols = getColors;

f = figure;
f.Position = [680   678   233   200];
f.Renderer = 'painters';
ax = gca;
ax  = prettifyPlot(ax);
hold on;
temp1 = mean(p.dr_wc(ix,:),1);
% temp1 = temp1(:);
temp2 = mean(p.wc(ix,:),1);
% temp2 = temp2(:);
mdl = fitlm(temp1,temp2);
pp = plot(mdl);
title(['R^2 = ' num2str(mdl.Rsquared.Ordinary)])
pp(1).Marker = 'o';
pp(1).MarkerSize = 7;
pp(1).MarkerEdgeColor = cols.potent;
% pp(1). = 'w';
pp(2).Color = 'k';
pp(2).LineStyle = '--';

xs = pp(2).XData;
y1 = pp(2).YData;
y2 = pp(3).YData;
% ax = gca; hold on;
% shadedErrorBar(xs,y1,y2-y1,{'Color','k'},0.2,ax)
pp(3).Visible ='off';
pp(4).Visible = 'off';

xlabel('Potent DR WC')
ylabel('Potent WC')
% ax.FontSize = 11;
ax.Legend.Visible = 'off';


%% potent dr_wc and potent wc
close all
cols = getColors;

f = figure;
f.Position = [680   678   233   200];
f.Renderer = 'painters';
ax = gca;
ax  = prettifyPlot(ax);
hold on;
temp1 = mean(p.dr_wc(ix,:),1);
% temp1 = temp1(:);
temp2 = mean(p.dr_wc_resp(ix,:),1);
% temp2 = temp2(:);
mdl = fitlm(temp1,temp2);
pp = plot(mdl);
title(['R^2 = ' num2str(mdl.Rsquared.Ordinary)])
pp(1).Marker = 'o';
pp(1).MarkerSize = 5;
pp(1).MarkerFaceColor = 'none';
pp(1).MarkerEdgeColor = cols.potent;
pp(2).Color = 'k';
pp(2).LineStyle = '--';

xs = pp(2).XData;
y1 = pp(2).YData;
y2 = pp(3).YData;
% ax = gca; hold on;
% shadedErrorBar(xs,y1,y2-y1,{'Color','k'},0.2,ax)
pp(3).Visible ='off';
pp(4).Visible = 'off';

xlabel('Potent DR WC')
ylabel('Potent DR WC Resp')
% ax.FontSize = 11;
ax.Legend.Visible = 'off';



