function plotPerformanceAllMice(meta,obj,rez,dfparams,params,cond2use,connect)

perf = zeros(numel(rez),numel(cond2use)); % (anm,cond)
for i = 1:numel(rez)
    perf(i,:) = mean( rez(i).perf(:,cond2use) , 1 ); % mean across sessions
end
perf = perf .* 100;


%%
f = figure;
f.Position = [680   716   292   262];
f.Renderer = 'painters';
ax = gca;
ax = prettifyPlot(ax);
hold on;

dfparams.plt.color{1}     = [0.0392 0.0392 0.0392];
dfparams.plt.color{2} = [0.4706 0.4588 0.4588];
dfparams.plt.color{3} = [0 0 1];
dfparams.plt.color{4} = [0 0 1];
dfparams.plt.color{5} = [1 0 0];
dfparams.plt.color{6} = [1 0 0];

xs = [1 2.5 4.5 6 8 9.5];
ct = 1;
for i = 1:size(perf,2)
    b(i) = bar(xs(i),mean(perf(:,i)));
    if mod(ct,2)~=0
        b(i).FaceColor = dfparams.plt.color{i};
        b(i).EdgeColor = dfparams.plt.color{i};
        b(i).FaceAlpha = 1;
    else
        b(i).FaceColor = 'none';
        b(i).EdgeColor = dfparams.plt.color{i};
        b(i).LineWidth = 1;
        b(i).FaceAlpha = 1;
    end

%     vs(i) = scatter(xs(i)*ones(size(perf(:,i))),perf(:,i),10,'MarkerFaceColor','k',...
%         'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25);
    
    xx = randn(size(perf(:,i))) * 0.01 + xs(i)*ones(size(perf(:,i)));
    vs(i) = scatter(xx,perf(:,i),30,'MarkerFaceColor','k',...
        'MarkerEdgeColor','w','LineWidth',1);

    xs_(:,i) = vs(i).XData';
    ys_(:,i) = perf(:,i);

    e = errorbar(b(i).XEndPoints,mean(perf(:,i)),std(perf(:,i)),'LineStyle','none','Color','k','LineWidth',1);
    e.LineWidth = 0.5;
    e.CapSize = 2;

    ct = ct + 1;
end

for i = 1:size(xs_,1)
    line(xs_(i,1:2),ys_(i,1:2),'LineWidth',0.1,'color','k')
    line(xs_(i,3:4),ys_(i,3:4),'LineWidth',0.1,'color','k')
    line(xs_(i,5:6),ys_(i,5:6),'LineWidth',0.1,'color','k')
    % patchline(xs_(i,1:2),ys_(i,1:2),'EdgeAlpha',0.4,'LineWidth',0.1)
    % patchline(xs_(i,3:4),ys_(i,3:4),'EdgeAlpha',0.4,'LineWidth',0.1)
    % patchline(xs_(i,5:6),ys_(i,5:6),'EdgeAlpha',0.4,'LineWidth',0.1)
end

xticks(xs)
xticklabels(["All ctrl" "All stim"  "Right ctrl" "Right stim"  "Left ctrl" "Left stim"])
ylabel("Performance (%)")
ylim([0,100])
ax = gca;
% ax.FontSize = 12;

anms = strjoin(unique({meta.anm}),' | ');
% stims = strjoin(unique({meta.stim}),' | ');
locs = strjoin(unique({meta.stimLoc}),' | ');
locs = strrep(locs,'_',' ');

% title([anms ' - ' locs],'fontsize',9)




end






