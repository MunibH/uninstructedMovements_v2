% rt = firstTongueRT(obj);
rt = firstJawRT(obj);

%%

close all
clear dat cols trialsBlock trix

for sessix = 1:numel(meta)

    t1 = -0.4; t2 = -0.01;
    [~,ix1] = min(abs(obj(sessix).time - t1));
    [~,ix2] = min(abs(obj(sessix).time - t2));

    trix{1} = params(sessix).trialid{2};
    trix{2} = params(sessix).trialid{3};
    alltrix = cell2mat(params(sessix).trialid(2:3)');

    for i = 1:numel(trix)
        dat.rt{i}{sessix} = rt{sessix}(trix{i});
    end
    dat.rt_alltrix{sessix} = rt{sessix}(alltrix)';
    dat.me{sessix} = mean(me(sessix).data(ix1:ix2,alltrix),1)';


end

for i = 1:numel(trix)
    dat.allrt{i} = cell2mat(dat.rt{i});
end

% correlate rt and late delay motion energy for each session
dat.rt_me_corr = cell2mat(cellfun(@corr, dat.me, dat.rt_alltrix, 'Uni',0));

%% trial type


f=figure; hold on;
ax = gca;
xs = [1 3];
clrs = getColors;
cols = {clrs.rhit, clrs.lhit};
div = 1.3;
for i = 1:numel(xs)
    temp = cell2mat(dat.rt{i});
    h(i) = bar(xs(i), nanmean(temp));
    h(i).FaceColor = cols{i};
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.5;
    mus = cell2mat(cellfun(@(x) nanmean(x), dat.rt{i}, 'UniformOutput',false));
    sem = cell2mat(cellfun(@(x) nanstd(x), dat.rt{i}, 'UniformOutput',false)) ./ numel(meta);
    scatter(xs(i)*ones(size(mus)),mus,30,'MarkerFaceColor',cols{i}./div, ...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',.35, ...
        'MarkerFaceAlpha',0.7)
    errorbar(h(i).XEndPoints,nanmean(temp),nanstd(temp)./sqrt(numel(meta)),'LineStyle','none','Color','k','LineWidth',1);
end
ax.XTick = xs;
ax.XTickLabel = {'Right Hit AFC', 'Left Hit AFC'};
ylabel('Reaction Time (s)')
ax.FontSize = 13;


%% motion energy and rt

% close all

figure;
for sessix = 1:numel(meta)
    
    x = dat.rt_alltrix{sessix};
    y = dat.me{sessix};

    mdl = fitlm(x,y);
    
    r2(sessix) = mdl.Rsquared.Ordinary;

%     plot(mdl)
%     title(['R2 = ' num2str(r2(sessix))])
%     xlabel('Reaction time (s)')
%     ylabel('Late delay motion energy (a.u.)');
%     pause
%     clf

end


f=figure; hold on;
ax = gca;
xs = [1];
div = 1.3;
for i = 1:numel(xs)
    temp = r2;
    h(i) = bar(xs(i), nanmean(temp));
    h(i).FaceColor = 'k';
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.5;
    scatter(xs(i)*ones(size(temp)),temp,30,'MarkerFaceColor','k', ...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',.35, ...
        'MarkerFaceAlpha',0.7)
    errorbar(h(i).XEndPoints,nanmean(temp),nanstd(temp)./sqrt(numel(meta)),'LineStyle','none','Color','k','LineWidth',1);
end
ax.XTick = xs;
ax.XTickLabel = {' '};
ylabel('R-squared (ME,RT)')
ax.FontSize = 13;





















