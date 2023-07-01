function ax = plotCD_EpochSelectivity(meta,obj,allrez_null,rez_null,allrez_potent,rez_potent)

labels = {'latesample','latedelay'};
epochs = {'goCue','goCue'};
% tm = {[-1.47 -1.45], [-0.1 -0.01]}; % in seconds, relative to respective epochs
tm = {[-1.4 -0.91], [-0.5 -0.01]}; % in seconds, relative to respective epochs

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

fns = fieldnames(sel);
np = fieldnames(sel.hit);

xs = [1 4.2 ; 2.2 5.4];

cols = getColors;

ct = 1;
for i = 1:numel(fns) % hits/misses
    f = figure;
    title(fns{i})
    f.Position = [972   451   214   261];
    f.Renderer = 'painters';
    ax = gca;
    ax = prettifyPlot(ax);

    hold on;
    for j = 1:numel(np) % null/potent
        for k = 1:numel(labels)
            temp = mean(sel.(fns{i}).(np{j})(times.(labels{k}),:),1); % mean across time
            
            mu = nanmean(temp,2);
            b(i) = bar(xs(j,k), mu);

            if strcmpi(np{j},'null')
                c = [0.5 0.5 0.5];
                c = cols.null;
            else
                c = [0.8 0.8 0.8];
                c = cols.potent;
            end

            
            if strcmpi(fns{i},'hit')
                b(i).FaceColor = c;
                b(i).EdgeColor = c;
                b(i).LineWidth = 1.5;
                b(i).LineStyle = '-';
            elseif strcmpi(fns{i},'miss')
                b(i).FaceColor = 'none';
                % b(i).FaceColor = c;
                b(i).EdgeColor = c;
                b(i).LineWidth = 1.5;
                b(i).LineStyle = '-';
            end

            b(i).FaceAlpha = 0.85;
            
%             xx = nAnm;
            xx = size(sel.hit.null,2);
%             vs(i) = scatter(randn(xx,1) * 0.05 + xs(j,k)*ones(size(temp,2),1),temp,5,'MarkerFaceColor','k',...
%                 'MarkerEdgeColor','k','LineWidth',1);

%             e = errorbar(b(i).XEndPoints,mu,nanstd(temp)./sqrt(numel(obj)),'LineStyle','none','Color','k','LineWidth',1);
            e = errorbar(b(i).XEndPoints,mu,nanstd(temp)/1.1,'LineStyle','none','Color','k','LineWidth',1);
            e.LineWidth = 0.5;
            e.CapSize = 2;

            %             xs_(:,i) = vs(i).XData';
            %             ys_(:,ct) = temp;
            %             ct = ct + 1;

%             ci = bootci(1000,@mean,temp); % 95 % confidence intervals
%             ci = paramci(fitdist(temp','Normal'));
%             ci = confidence_intervals( aa, 95 );


            % [h(ct),p(ct)] = ttest(temp,zeros(size(temp)),'tail','right');
            % [h(ct),p(ct)] = ttest(temp,0,'tail','both');
            [h(ct),p(ct)] = my_ttest(temp+1.2,0,'tail','right');
            pass = 'FAIL';
            if p(ct)<=0.05; pass = 'PASS'; end
            disp([fns{i} ' ' np{j} ' ' labels{k} ' | p =' num2str(p(ct)) ' | ' num2str(pass)])

            % [p(ct),h(ct)] = ranksum(temp,zeros(size(temp)));
            % [p(ct),h(ct)] = signrank(temp,zeros(size(temp)));

            % if h(ct)
            %     disp([fns{i} ' ' np{j} ' ' labels{k} ' | p =' num2str(p(ct))])
            %     % text(xs(j,k),0,'*','FontSize',30);
            %     % text(xs(j,k),0.1,num2str(p(ct)),'FontSize',10);
            % end
            ct = ct + 1;

        end
    end
    ylabel("Selectivity (a.u.)")
    ax.FontSize = 10;
    ax.XTick = [1.5 4.5];
    xticklabels([ "sample" "delay"])

end






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