function plotPerformanceAllMice(meta,obj,rez,dfparams,params,cond2use,connect)

perf = zeros(numel(rez),numel(cond2use)); % (anm,cond)
for i = 1:numel(rez)
    perf(i,:) = mean( rez(i).perf(:,cond2use) , 1 ); % mean across sessions
end
perf = perf .* 100;


%%
f = figure;
f.Position = [680   716   292   262];
ax = gca;
hold on;

xs = [1 2 4 5 7 8];
ct = 1;
for i = 1:size(perf,2)
    b(i) = bar(xs(i),mean(perf(:,i)));
    if mod(ct,2)~=0
        b(i).FaceColor = dfparams.plt.color{i};
        b(i).EdgeColor = 'none';
        b(i).FaceAlpha = 0.8;
    else
        b(i).FaceColor = 'none';
        b(i).EdgeColor = dfparams.plt.color{i};
        b(i).LineWidth = 2;
        b(i).FaceAlpha = 1;
    end

%     vs(i) = scatter(xs(i)*ones(size(perf(:,i))),perf(:,i),10,'MarkerFaceColor','k',...
%         'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25);
    
    xx = randn(size(perf(:,i))) * 0.01 + xs(i)*ones(size(perf(:,i)));
    vs(i) = scatter(xx,perf(:,i),10,'MarkerFaceColor','k',...
        'MarkerEdgeColor','k','LineWidth',1);

    xs_(:,i) = vs(i).XData';
    ys_(:,i) = perf(:,i);

    e = errorbar(b(i).XEndPoints,mean(perf(:,i)),std(perf(:,i)),'LineStyle','none','Color','k','LineWidth',1);
    e.LineWidth = 0.5;
    e.CapSize = 2;

    ct = ct + 1;
end

for i = 1:size(xs_,1)
    patchline(xs_(i,1:2),ys_(i,1:2),'EdgeAlpha',0.4,'LineWidth',0.1)
    patchline(xs_(i,3:4),ys_(i,3:4),'EdgeAlpha',0.4,'LineWidth',0.1)
    patchline(xs_(i,5:6),ys_(i,5:6),'EdgeAlpha',0.4,'LineWidth',0.1)
end

xticks(xs)
xticklabels(["All ctrl" "All stim"  "Right ctrl" "Right stim"  "Left ctrl" "Left stim"])
ylabel("Performance (%)")
ylim([0,100])
ax = gca;
ax.FontSize = 12;

anms = strjoin(unique({meta.anm}),' | ');
% stims = strjoin(unique({meta.stim}),' | ');
locs = strjoin(unique({meta.stimLoc}),' | ');
locs = strrep(locs,'_',' ');

title([anms ' - ' locs],'fontsize',9)

% if connect
%     % for each pair of conditions, plot connecting lines
%     hold on;
%     for i = 1:2:(numel(cond2use)-1)
%         conds = [i i+1];
%         xs = [vs(conds(1)).XData ; vs(conds(2)).XData];
%         ys = [vs(conds(1)).YData ; vs(conds(2)).YData];
%         patchline(xs, ys,'EdgeColor','k','EdgeAlpha',0.3)
%
% %         plot(xs,ys,'k-')
%     end
% end
%
%
%
%
% % %%
% %
% % if numel(meta) == 1
% %     f = figure;
% %     f.Position = [316          73        1232         905];
% %     ax = axes(f); hold on;
% %     for i = 1:numel(perf)
% %         b(i) = bar(categorical(dfparams.cond(cond2use(i))),perf(i));
% %         b(i).FaceColor = dfparams.plt.color{cond2use(i)};
% %         b(i).EdgeColor = 'none';
% %     end
% %     ylabel('Performance (%)')
% %     ylim([0,100])
% %     title('Performance on all sessions, all mice')
% %     ax.FontSize = 13;
% %     hold off;
% %     return
% % end
% %
% %
% % f = figure;
% % f.Position = [316          73        1232         905];
% % ax = axes(f);
% % violincols = reshape(cell2mat(dfparams.plt.color(cond2use)),3,numel(cond2use))';
% % vs = violinplot(perf,dfparams.cond(cond2use),...
% %     'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1}, 'ViolinColor', violincols,'Width',0.2);
% % ylabel('Performance (%)')
% % ylim([0,100])
% % title('Performance on all sessions, all mice')
% % ax.FontSize = 13;
% %
% %
% % if connect
% %     % for each pair of conditions, plot connecting lines
% %     hold on;
% %     for i = 1:2:(numel(cond2use)-1)
% %         conds = [i i+1];
% %         xs = [vs(conds(1)).ScatterPlot.XData ; vs(conds(2)).ScatterPlot.XData];
% %         ys = [vs(conds(1)).ScatterPlot.YData ; vs(conds(2)).ScatterPlot.YData];
% %         patchline(xs, ys,'EdgeColor','k','EdgeAlpha',0.3)
% %
% % %         plot(xs,ys,'k-')
% %     end
% % end
% %
% % t-test b/w groups
% ct = 1;
% for i = 1:2:(numel(cond2use)-1)
%     conds = [i i+1];
%     xs = [vs(conds(1)).XData ; vs(conds(2)).XData];
%     ys = [vs(conds(1)).YData ; vs(conds(2)).YData];
%     %     [h(ct),p(ct)] = ttest2(ys(1,:),ys(2,:));
%     [p(ct),h(ct)] = ranksum(ys(1,:),ys(2,:)); % if unsure both datasets have equal variance, use mann u whtiney instead of t test
%
%     lw = 0.4;
%     x = xs(i:i+1);
%
%     if h(ct)
%         plot(mean(x),95,'k*')
%     else
% %         text(mean(x),95,'n.s.','FontSize',9,'FontWeight','bold')
%     end
%
%     ct = ct + 1;
%
%
% end



end






