function ax = plotCD_EpochSelectivity(meta,obj,allrez_null,rez_null,allrez_potent,rez_potent)

labels = {'latesample','latedelay'};
epochs = {'goCue','goCue'};
% tm = {[-1.47 -1.45], [-0.1 -0.01]}; % in seconds, relative to respective epochs
tm = {[-1.5 -1.4], [-0.1 -0.01]}; % in seconds, relative to respective epochs

for i = 1:numel(epochs)
    e1 = mode(obj(1).bp.ev.(epochs{i})) + tm{i}(1) - 2.5;
    e2 = mode(obj(1).bp.ev.(epochs{i})) + tm{i}(2) - 2.5;
    times.(labels{i}) = rez_null(1).time>e1 & rez_null(1).time<e2;
end



%%

% if strcmpi(spacename,'null')
%     c = [0.4 0.4 0.4];
% elseif strcmpi(spacename,'potent')
%     c = [0.7 0.7 0.7];
% end

[objix,uAnm] = groupSessionsByAnimal(meta);
nAnm = numel(uAnm);


cdix = find(ismember(rez_null(1).cd_labels,'late'));

sel = struct();
% Null
sel = getEpochSelectivity(sel,cdix,allrez_null,objix,nAnm,'null');

% Potent
sel = getEpochSelectivity(sel,cdix,allrez_potent,objix,nAnm,'potent');

%%

% tt = times.latesample;
% % tt = times.latedelay;
% % cc = [1 2]; % hits
% cc = [3 4]; % miss
% dat = squeeze(allrez_potent.cd_proj(:,cc,1,:));
% rpotent = mean(squeeze(dat(tt,1,:)),1);
% lpotent = mean(squeeze(dat(tt,2,:)),1);
% sel.prctile.potent = prctile(rpotent - lpotent, [5 95]);
% 
% 
% dat = squeeze(allrez_null.cd_proj(:,cc,1,:));
% rnull = mean(squeeze(dat(tt,1,:)),1);
% lnull = mean(squeeze(dat(tt,2,:)),1);
% sel.prctile.null = prctile(rnull - lnull, [5 95]);
% 
% % figure; hold on;
% % histogram(r,100,'facecolor','b')
% % histogram(l,100,'facecolor','r')
% 
% cols = getColors;
% 
% figure;
% hold on;
% nbins = 30;
% histogram(rpotent-lpotent,nbins,'facecolor',cols.potent,'facealpha',0.6,'edgecolor','none')
% xline(sel.prctile.potent(1),'--','color',cols.potent,'linewidth',2)
% xline(sel.prctile.potent(2),'--','color',cols.potent,'linewidth',2)
% histogram(rnull-lnull,nbins,'facecolor',cols.null,'facealpha',0.6,'edgecolor','none')
% xline(sel.prctile.null(1),'--','color',cols.null,'linewidth',2)
% xline(sel.prctile.null(2),'--','color',cols.null,'linewidth',2)
% title('late sample | miss')
% xline(0,'--','linewidth',2)
% 
% % pd = fitdist((r-l)','Normal');
% % xs = -2.5:0.01:1.5;
% % y = pdf(pd,xs);
% % figure;
% % plot(xs,y)

%%

fns = fieldnames(sel);
np = fieldnames(sel.hit);

xs = [1 4 ; 2 5];

ct = 1;
for i = 1:numel(fns) % hits/misses
    f = figure;
    title(fns{i})
    f.Position = [972   451   214   261];
    ax = gca;
    hold on;
    for j = 1:numel(np) % null/potent
        for k = 1:numel(labels)
            temp = mean(sel.(fns{i}).(np{j})(times.(labels{k}),:),1); % mean across time

            mu = nanmean(temp,2);
            b(i) = bar(xs(j,k), mu);

            if strcmpi(np{j},'null')
                c = [0.5 0.5 0.5];
            else
                c = [0.8 0.8 0.8];
            end

            
            if strcmpi(fns{i},'hit')
                b(i).FaceColor = c;
                b(i).EdgeColor = 'none';
                b(i).LineWidth = 0.01;
                b(i).LineStyle = '-';
            elseif strcmpi(fns{i},'miss')
                b(i).FaceColor = 'none';
                b(i).EdgeColor = c;
                b(i).LineWidth = 2;
                b(i).LineStyle = '-';
            end

            b(i).FaceAlpha = 0.85;
            
%             xx = nAnm;
            xx = size(sel.hit.null,2);
%             vs(i) = scatter(randn(xx,1) * 0.05 + xs(j,k)*ones(size(temp,2),1),temp,5,'MarkerFaceColor','k',...
%                 'MarkerEdgeColor','k','LineWidth',1);

%             e = errorbar(b(i).XEndPoints,mu,nanstd(temp)./sqrt(numel(obj)),'LineStyle','none','Color','k','LineWidth',1);
            e = errorbar(b(i).XEndPoints,mu,nanstd(temp),'LineStyle','none','Color','k','LineWidth',1);
            e.LineWidth = 0.5;
            e.CapSize = 2;

            %             xs_(:,i) = vs(i).XData';
            %             ys_(:,ct) = temp;
            %             ct = ct + 1;

%             ci = bootci(1000,@mean,temp); % 95 % confidence intervals
%             ci = paramci(fitdist(temp','Normal'));
%             ci = confidence_intervals( aa, 95 );


%             [h(ct),p(ct)] = ttest(temp,zeros(size(temp)));
%             [h(ct),p(ct)] = my_ttest(temp);
% %             [p(ct),h(ct)] = ranksum(temp,zeros(size(temp)));
% %             [p(ct),h(ct)] = signrank(temp,zeros(size(temp)));
%             if h(ct)
%                 text(xs(j,k),0,'*','FontSize',30);
%                 text(xs(j,k),0.1,num2str(p(ct)),'FontSize',10);
%             end
            ct = ct + 1;

        end
    end
    ylabel("Selectivity (a.u.)")
    ax.FontSize = 10;
    ax.XTick = [1.5 4.5];
    xticklabels([ "sample" "delay"])

end




%%
%
% f = figure;
% f.Position = [680   691   294   287];
% ax = gca;
% hold on;
%
% xs = [1 2 ; 4 5];
%
% fns = fieldnames(sel);
% ct = 1;
% for i = 1:numel(labels)
%     for j = 1:numel(fns)
%         temp = mean(sel.(fns{j}).(np{j})(times.(labels{i}),:),1); % mean across time
%
%         mu = nanmean(temp,2);
%
%         b(i) = bar(xs(i,j), mu);
%         if j == 2
%             b(i).FaceColor = 'none';
%             b(i).EdgeColor = c;
%             b(i).LineWidth = 2;
%         else
%             b(i).FaceColor = c;
%             b(i).EdgeColor = 'none';
%         end
%
%         b(i).FaceAlpha = 0.85;
%
%         vs(i) = scatter(randn(nAnm,1) * 0.05 + xs(i,j)*ones(size(temp,2),1),temp,5,'MarkerFaceColor','k',...
%             'MarkerEdgeColor','k','LineWidth',1);
%
%         b(i).XEndPoints
%         e = errorbar(b(i).XEndPoints,mu,nanstd(temp)./sqrt(numel(obj)),'LineStyle','none','Color','k','LineWidth',1);
%         e.LineWidth = 0.5;
%         e.CapSize = 2;
%
%         xs_(:,i) = vs(i).XData';
%         ys_(:,ct) = temp;
%         ct = ct + 1;
%
%
%         %         [h(ct),p(ct)] = ttest(temp);
%         [p(ct),h(ct)] = ranksum(temp,zeros(size(temp)));
%         if h(ct)
%             text(xs(i,j),-1,'*','FontSize',20)
%         end
%
%     end
% end
%
%
% % xs_ = repmat([1 2 4 5]',1,nAnm)';
% % for i = 1:size(xs_,1)
% %     patchline(xs_(i,1:2),ys_(i,1:2),'EdgeAlpha',0.4,'LineWidth',0.1)
% %     patchline(xs_(i,3:4),ys_(i,3:4),'EdgeAlpha',0.4,'LineWidth',0.1)
% % end
%
% ylabel("Selectivity (a.u.)")
% ax.FontSize = 10;
% ax.XTick = [1 2 4 5];
% xticklabels([ "sample-hit" "sample-miss" "delay-hit" "delay-miss"])
% title(spacename)
%





end


%%

function sel = getEpochSelectivity(sel,cdix,allrez,objix,nAnm,spacename)


cddat = squeeze(allrez.cd_proj(:,:,cdix,:));

sel.hit.(spacename) = squeeze(cddat(:,1,:) - cddat(:,2,:));
sel.miss.(spacename) = squeeze(cddat(:,3,:) - cddat(:,4,:));


% selfns = fieldnames(sel);
% for i = 1:numel(selfns)
%     for ianm = 1:nAnm
%         sessix = find(objix{ianm});
%         temp = nanmean(sel.(selfns{i}).(spacename)(:,sessix),2);
%         newsel.(selfns{i}).(spacename)(:,ianm) = temp;
%     end
% end
% 
% 
% sel.hit.(spacename) = newsel.hit.(spacename);
% sel.miss.(spacename) = newsel.miss.(spacename);


end