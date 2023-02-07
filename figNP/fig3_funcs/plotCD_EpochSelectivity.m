function ax = plotCD_EpochSelectivity(meta,obj,allrez,rez,spacename)

labels = {'latesample','latedelay'};
epochs = {'delay','goCue'};
tm = {[-0.82 -0.42], [-0.22 -0.02]}; % in seconds, relative to respective epochs

for i = 1:numel(epochs)
    e1 = mode(obj(1).bp.ev.(epochs{i})) + tm{i}(1) - 2.5;
    e2 = mode(obj(1).bp.ev.(epochs{i})) + tm{i}(2) - 2.5;
    times.(labels{i}) = rez(1).time>e1 & rez(1).time<e2;
end



%%

lw = 2.5;
alph = 0.1;

sample = mode(rez(1).ev.sample - rez(1).align);
delay = mode(rez(1).ev.delay - rez(1).align);


if strcmpi(spacename,'null')
    c = [0.4 0.4 0.4];
elseif strcmpi(spacename,'potent')
    c = [0.7 0.7 0.7];
end

cdix = find(ismember(rez(1).cd_labels,'late'));

cddat = squeeze(allrez.cd_proj(:,:,cdix,:));

sel.hit = squeeze(cddat(:,1,:) - cddat(:,2,:));
sel.miss = squeeze(cddat(:,3,:) - cddat(:,4,:));


[objix,uAnm] = groupSessionsByAnimal(meta);
nAnm = numel(uAnm);

selfns = fieldnames(sel);
for i = 1:numel(selfns)
    for ianm = 1:nAnm
        sessix = find(objix{ianm});
        temp = nanmean(sel.(selfns{i})(:,sessix),2);
        newsel.(selfns{i})(:,ianm) = temp;
    end
end

clear sel;
sel = newsel;
clear newsel;

%%

f = figure;
f.Position = [680   691   294   287];
ax = gca;
hold on;

xs = [1 2 ; 4 5];

fns = fieldnames(sel);
ct = 1;
for i = 1:numel(labels)
    for j = 1:numel(fns)
        temp = mean(sel.(fns{j})(times.(labels{i}),:),1); % mean across time

        mu = nanmean(temp,2);

        b(i) = bar(xs(i,j), mu);
        if j == 2
            b(i).FaceColor = 'none';
            b(i).EdgeColor = c;
            b(i).LineWidth = 2;
        else
            b(i).FaceColor = c;
            b(i).EdgeColor = 'none';
        end

        b(i).FaceAlpha = 0.85;

        vs(i) = scatter(randn(nAnm,1) * 0.1 + xs(i,j)*ones(size(temp,2),1),temp,5,'MarkerFaceColor','k',...
            'MarkerEdgeColor','k','LineWidth',1);

        b(i).XEndPoints
        e = errorbar(b(i).XEndPoints,mu,nanstd(temp),'LineStyle','none','Color','k','LineWidth',1);
        e.LineWidth = 0.5;
        e.CapSize = 2;

        xs_(:,i) = vs(i).XData';
        ys_(:,ct) = temp;
        ct = ct + 1;


%         [h(ct),p(ct)] = ttest(temp);
        [p(ct),h(ct)] = ranksum(temp,zeros(size(temp)));
        if h(ct)
            text(xs(i,j),-1,'*','FontSize',20)
        end

    end
end


% xs_ = repmat([1 2 4 5]',1,nAnm)';
% for i = 1:size(xs_,1)
%     patchline(xs_(i,1:2),ys_(i,1:2),'EdgeAlpha',0.4,'LineWidth',0.1)
%     patchline(xs_(i,3:4),ys_(i,3:4),'EdgeAlpha',0.4,'LineWidth',0.1)
% end

ylabel("Selectivity (a.u.)")
ax.FontSize = 10;
ax.XTick = [1 2 4 5];
xticklabels([ "sample-hit" "sample-miss" "delay-hit" "delay-miss"])
title(spacename)






end