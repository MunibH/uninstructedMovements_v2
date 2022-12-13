
% rt = firstTongueRT(obj);
rt = firstJawRT(obj);

%%

close all
clear dat cols trialsBlock

nQuartiles = 4;
for sessix = 1:numel(meta)

    %     t1 = -0.5; t2 = -0.01;
    %     [~,ix1] = min(abs(obj(sessix).time - t1));
    %     [~,ix2] = min(abs(obj(sessix).time - t2));

    potent = sum(rez(sessix).N_potent.^2,3);
    cdlatepotent = cd_potent(sessix).trialdat(:,:,2);
    null = sum(rez(sessix).N_null.^2,3);
    cdlatenull = cd_null(sessix).trialdat(:,:,2);
    me_ = me(sessix).data;
    thisrt = rt{sessix};

    %     figure; subplot(2,1,1); histogram(thisrt,40); subplot(2,1,2); plot(thisrt); pause; close all; continue;

    % right and left hit trials
    trix.all = cell2mat(params(sessix).trialid(2:5)');
    thisrt = thisrt(trix.all);
    quartiles = quantile(thisrt,nQuartiles);
    trialsBlock{sessix} = getTrialsByBlock(thisrt,quartiles);

    trialsBlock{sessix} = trialsBlock{sessix}; % remove last quartile, verry slow RTs
%     trialsBlock{sessix} = trialsBlock{sessix}(1:4); % remove last quartile, verry slow RTs

    for i = 1:numel(trialsBlock{sessix})
        dat.rt{i}{sessix} = thisrt(trialsBlock{sessix}{i});
        dat.me{i}{sessix} = me_(:,trialsBlock{sessix}{i});
        dat.cdpotent{i}{sessix} = cdlatepotent(:,trialsBlock{sessix}{i});
        dat.cdnull{i}{sessix} = cdlatenull(:,trialsBlock{sessix}{i});
        dat.potent{i}{sessix} = potent(:,trialsBlock{sessix}{i});
        dat.null{i}{sessix} = null(:,trialsBlock{sessix}{i});
    end


end




%%

clear null potent me_ cdnull cdpotent

for i = 1:numel(dat.rt) % each rt group
    for sessix = 1:numel(dat.rt{i}) % each session
        cdnull(:,sessix,i) = mean(dat.cdnull{i}{sessix},2);
        cdpotent(:,sessix,i) = mean(dat.cdpotent{i}{sessix},2);
        null(:,sessix,i) = mean(dat.null{i}{sessix},2);
        potent(:,sessix,i) = mean(dat.potent{i}{sessix},2);
        me_(:,sessix,i) = mean(dat.me{i}{sessix},2);
    end
end

c = gray(nQuartiles*4);
for i = 1:nQuartiles
    cols(i,:) = c(1+(i-1)*2,:);
end
% cols(1,:) = c(1,:);
% cols(2,:) = c(3,:);
% cols(3,:) = c(5,:);
% cols(4,:) = c(7,:);

figure;
ax = gca;
hold on;
for i = 1:size(me_,3)
    temp = me_(:,:,i);
    plot(obj(1).time,mean(temp,2),'Color',cols(i,:),'LineWidth',1)
    %     shadedErrorBar(obj(1).time,mean(temp,2),std(temp,[],2)./sqrt(size(temp,2)),{'Color','k','LineWidth',1},0.01,ax);
end
xlim([-2.4,2.4])
ylabel('Motion Energy')


figure;
ax = gca;
hold on;
for i = 1:size(cdnull,3)
    temp = cdnull(:,:,i);
    plot(obj(1).time,mean(temp,2),'Color',cols(i,:),'LineWidth',1)
    %     shadedErrorBar(obj(1).time,mean(temp,2),std(temp,[],2)./sqrt(size(temp,2)),{'Color','k','LineWidth',1},0.01,ax);
end
xlim([-2.4,2.4])
ylabel('CD Late (Null)')

figure;
ax = gca;
hold on;
for i = 1:size(cdpotent,3)
    temp = cdpotent(:,:,i);
    plot(obj(1).time,mean(temp,2),'Color',cols(i,:),'LineWidth',1)
    %     shadedErrorBar(obj(1).time,mean(temp,2),std(temp,[],2)./sqrt(size(temp,2)),{'Color','k','LineWidth',1},0.01,ax);
end
xlim([-2.4,2.4])
ylabel('CD Late (Potent)')

figure;
ax = gca;
hold on;
for i = 1:size(null,3)
    temp = null(:,:,i);
    plot(obj(1).time,mean(temp,2),'Color',cols(i,:),'LineWidth',1)
    %     shadedErrorBar(obj(1).time,mean(temp,2),std(temp,[],2)./sqrt(size(temp,2)),{'Color','k','LineWidth',1},0.01,ax);
end
xlim([-2.4,2.4])
ylabel('Null')

figure;
ax = gca;
hold on;
for i = 1:size(potent,3)
    temp = potent(:,:,i);
    plot(obj(1).time,mean(temp,2),'Color',cols(i,:),'LineWidth',1)
    %     shadedErrorBar(obj(1).time,mean(temp,2),std(temp,[],2)./sqrt(size(temp,2)),{'Color','k','LineWidth',1},0.01,ax);
end
xlim([-2.4,2.4])
ylabel('Potent')

%%

f=figure; hold on;
ax = gca;
xs = 1:numel(quartiles);
for i = 1:numel(xs)
    temp = cell2mat(dat.rt{i});
    h(i) = bar(xs(i), nanmean(temp));
    h(i).FaceColor = 'k';
    h(i).EdgeColor = 'none';
    h(i).FaceAlpha = 0.5;
    mus = cell2mat(cellfun(@(x) nanmean(x), dat.rt{i}, 'UniformOutput',false));
    sem = cell2mat(cellfun(@(x) nanstd(x), dat.rt{i}, 'UniformOutput',false)) ./ numel(meta);
    scatter(xs(i)*ones(size(mus)),mus,30,'MarkerFaceColor','k', ...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',.35, ...
        'MarkerFaceAlpha',0.7)
    %     errorbar(h(i).XEndPoints,mus',sem','LineStyle','none','Color','k','LineWidth',1);
end
xlabel('Quantile')
ylabel('Reaction Time (s)')
ax.XTick = 1:numel(quartiles);
ax.FontSize = 13;


%% Helper functions

function trialsBlock = getTrialsByBlock(rt,quartiles)
for i = 1:numel(quartiles)
    if i == 1
        trialsBlock{i} = rt<=quartiles(i);
    elseif i == numel(quartiles)
        trialsBlock{i} = rt>quartiles(i);
    else
        trialsBlock{i} = rt>quartiles(i-1) & rt<=quartiles(i);
    end
end

% i = 1;
% trialsBlock{i}{sessix} = thisrt<=quartiles(i);
% i = 2;
% trialsBlock{i}{sessix} = thisrt>quartiles(i-1) & thisrt<=quartiles(i);
% i = 3;
% trialsBlock{i}{sessix} = thisrt>quartiles(i-1) & thisrt<=quartiles(i);
% i = 4;
% trialsBlock{i}{sessix} = thisrt>quartiles(i);

end
