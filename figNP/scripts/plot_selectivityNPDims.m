close all

cols = getColors;

mainCond = 1; % right trials should always be higher

for sessix = 1:numel(meta)
    t1 = mode(obj(sessix).bp.ev.bitStart - 2.5);
    t2 = mode(obj(sessix).bp.ev.sample - 2.5);
%     t1 = mode(obj(sessix).bp.ev.sample - 2.5);
%     t2 = mode(obj(sessix).bp.ev.sample - 2.5 + 0.5);
    [~,ix1] = min(abs(obj(sessix).time - t1));
    [~,ix2] = min(abs(obj(sessix).time - t2));


    trix{1} = params(sessix).trialid{2};
    trix{2} = params(sessix).trialid{3};

    for i = 1:numel(trix)
        null{i} = mySmooth(squeeze(mean(rez(sessix).N_null(:,trix{i},:),2)),1,'reflect');

        potent{i} = mySmooth(squeeze(mean(rez(sessix).N_potent(:,trix{i},:),2)),1,'reflect');
    end
    %     null = cellfun(@(x) smoothdata(x,10,"movmean",1), null, 'UniformOutput', false);
    sample = mode(obj(sessix).bp.ev.sample) - 2.5;
    tstart = mode(obj(sessix).bp.ev.bitStart) - 2.5;
    for i = 1:size(null{1},2)
        temp = abs(null{1}(:,i) - null{2}(:,i));
        ps.null{sessix}(i) = mean(temp(ix1:ix2));
%         dim = i;
%         figure;
%         subplot(3,1,1)
%         plot(obj(sessix).time,null{1}(:,dim),'b','LineWidth',1); hold on; plot(obj(sessix).time,null{2}(:,dim),'r','LineWidth',1)
%         xline(sample,'k--')
%         xline(tstart,'k--')
%         subplot(3,1,2)
%         plot(obj(sessix).time,null{1}(:,dim) - null{2}(:,dim),'m','LineWidth',1); hold on; yline(0,'k--')
%         xline(tstart,'k--')
%         xline(sample,'k--')
%         subplot(3,1,3)
%         plot(obj(sessix).time, abs(null{1}(:,dim) - null{2}(:,dim)),'m','LineWidth',1); hold on; yline(0,'k--')
%         xline(sample,'k--')
%         xline(tstart,'k--')
%         pause
    end

    for i = 1:size(potent{1},2)
        temp = abs(potent{1}(:,i) - potent{2}(:,i));
        ps.potent{sessix}(i) = mean(temp(ix1:ix2));
    end

    selnull = abs(null{1} - null{2});
    selpotent = abs(potent{1} - potent{2});

%     sel.null(:,sessix) = sum(selnull,2);
        sel.null(:,sessix) = sel.null(:,sessix) - mean(sel.null(ix1:ix2,sessix));
%     sel.potent(:,sessix) = sum(selpotent,2);
        sel.potent(:,sessix) = sel.potent(:,sessix) - mean(sel.potent(ix1:ix2,sessix));

    %     f = figure; f.Position = [839   431   392   500];
    %     imagesc(obj(sessix).time,1:size(selpotent,2),selpotent'); cpotent = colorbar; colormap(linspecer) % caxis(cnull.Limits);
    %     xlabel('Time (s) from go cue'); ylabel('Potent Dims'); ax = gca; ax.FontSize = 12; cpotent.Label.String = 'Selectivity';
    %     f = figure; f.Position = [839   431   392   500];
    %     imagesc(obj(sessix).time,1:size(selnull,2),selnull'); cnull = colorbar; colormap(linspecer); %caxis(cpotent.Limits);
    %     xlabel('Time (s) from go cue'); ylabel('Null Dims'); ax = gca; ax.FontSize = 12; cnull.Label.String = 'Selectivity';

%     break
    %     pause
    %     close all





end

%%
close all
lw = 2;
alph = 0.2;

f = figure;
ax = gca;
hold on;

means = mean(sel.null,2);
errs = std(sel.null,[],2) ./ sqrt(size(sel.null,2));
% CI95 = tinv([0.025 0.975], size(sel.null,2)-1);
% errs = means*CI95(2);
shadedErrorBar(obj(1).time,means,errs,{'Color',cols.null,'LineWidth',lw},alph,ax)

means = mean(sel.potent,2);
errs = std(sel.potent,[],2) ./ sqrt(size(sel.potent,2));
shadedErrorBar(obj(1).time,means,errs,{'Color',cols.potent,'LineWidth',lw},alph,ax)

sample = mode(obj(sessix).bp.ev.sample) - 2.5;
delay = mode(obj(sessix).bp.ev.delay) - 2.5;
xline(0,'k--')
xline(sample,'k--')
xline(delay,'k--')
xlim([-2.4 obj(1).time(end)])

ylim([-0.1 ax.YLim(2)])
ax.FontSize = 12;
xlabel('Time (s) from go cue')
ylabel('Selectivity (a.u.)')


%%

nullpresamplesel = cell2mat(ps.null);
potentpresamplesel = cell2mat(ps.potent);

figure;
ax = gca;
hold on;
bar(1,mean(nullpresamplesel));
scatter(1*ones(size(nullpresamplesel)),nullpresamplesel,10,'MarkerFaceColor','k', ...
            'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',.35, ...
            'MarkerFaceAlpha',0.7)
bar(2,mean(potentpresamplesel));
scatter(2*ones(size(potentpresamplesel)),potentpresamplesel,10,'MarkerFaceColor','k', ...
            'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',.35, ...
            'MarkerFaceAlpha',0.7)


ylabel('selectivity in presample period')
ax.XTick = [1 2];
xticklabels({'null','potent'})
ax.FontSize = 13;



