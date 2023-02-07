function [h,p] = npSelectivityByRTQuartile(meta,obj,rez,params,nQuantiles,rt,me)

% reuse this code but loop over quartiles

cond2use_ta = [1 2]; % right and left hits, corresponding to trial-avg projs onto n/p
cond2use_st = [8 9]; % right and left hits, corresponding to single-trial projs onto n/p

modparams.subTrials = 35;

% Sort by delay selectivity
edges = [mode(obj(1).bp.ev.delay) mode(obj(1).bp.ev.goCue)] - 2.5;



cond2use = [2 3];
% cond2use = [2 3 4 5]; % right, left hits , including early
% cond2use = [8 9 10 11]; % right, left hits, right miss, left miss , excluding early
% cond2use = [8 9];

epochedges = {[mode(obj(1).bp.ev.delay)-0.4-2.5 mode(obj(1).bp.ev.delay)-2.5 ], ...
              [mode(obj(1).bp.ev.delay)-2.5 mode(obj(1).bp.ev.delay)+0.4-2.5 ], ...
              [mode(obj(1).bp.ev.goCue)-0.4-2.5 mode(obj(1).bp.ev.goCue)-2.5 ]};

for iq = 1:numel(epochedges)
    for j = 1:numel(epochedges{iq})
        [~,epochix{iq}(j)] = min(abs(obj(1).time - epochedges{iq}(j)));
    end
end

for sessix = 1:numel(obj)
    trix.all = cell2mat(params(sessix).trialid(cond2use)');

    if numel(cond2use) > 2
        trix.right = cell2mat(params(sessix).trialid(cond2use([1 4]))');
        trix.left = cell2mat(params(sessix).trialid(cond2use([2 3]))');
    else
        trix.right = params(sessix).trialid{cond2use(1)};
        trix.left = params(sessix).trialid{cond2use(2)};
    end

    dat.rt = rt{sessix}(trix.all);
    dat.me = me(sessix).data(:,trix.all);
    dat.null = rez(sessix).N_null(:,trix.all,:);
    dat.potent = rez(sessix).N_potent(:,trix.all,:);

    quantiles = quantile(dat.rt,nQuantiles-1);
    trialsBlockMask = getTrialsByRTBlock(dat.rt,quantiles);
    %     cellfun(@(x) sum(x), trialsBlock, 'uniformoutput',false)
    trialsBlock = cellfun(@(x) trix.all(x), trialsBlockMask ,'UniformOutput',false);

    for iq = 1:nQuantiles
        modparams.iqTrials{1}  = trialsBlock{iq}( ismember(trialsBlock{iq}, trix.right) ); % right trials for current quantile
        modparams.iqTrials{2}  = trialsBlock{iq}( ismember(trialsBlock{iq}, trix.left) );  % left trials

        [nullsel,nulldims] = getSortedSelectivityForRT(obj(sessix),params(sessix),rez(sessix),edges,cond2use_ta,cond2use_st,modparams,'null');
        [potentsel,potentdims] = getSortedSelectivityForRT(obj(sessix),params(sessix),rez(sessix),edges,cond2use_ta,cond2use_st,modparams,'potent');

        dat.sel.null(:,iq,sessix) = nanmean(nullsel(:,:),2); % mean/sum across dimensions
        dat.sel.potent(:,iq,sessix) = nanmean(potentsel(:,:),2);

        allrt(iq,sessix) = mean(dat.rt(trialsBlockMask{iq}));

        allme(:,iq,sessix) = mean(dat.me(:,trialsBlockMask{iq}),2);
    end

    for iepoch = 1:numel(epochix)
        ixs = epochix{iepoch}(1):epochix{iepoch}(2);
        dat.sel.epoch.null(:,iepoch,sessix) = nanmean(dat.sel.null(ixs,:,sessix),1); % (quantile,epoch,session)
        dat.sel.epoch.potent(:,iepoch,sessix) = nanmean(dat.sel.potent(ixs,:,sessix),1); % (quantile,epoch,session)

%         dat.sel.epoch.null_time(:,:,iepoch,sessix) = dat.sel.null(ixs,:,sessix); %(time,quantile,epoch,session)
%         dat.sel.epoch.potent_time(:,:,iepoch,sessix) = dat.sel.null(ixs,:,sessix); %(time,quantile,epoch,session)
    end
    
    
end

%% plot n/p quantile trajectories (session-averaged)

close all

sample = mode(obj(1).bp.ev.sample) - 2.5;
tstart = mode(obj(1).bp.ev.bitStart) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;

sm = 51;

xlims = [tstart 1.5];

% c = gray;
% c = c(ceil(linspace(100,150,nQuantiles-1)),:);
% col(1,:) = [0 0 0];
% col(2:nQuantiles,:) = c;

clear cols
cols = getColors();
cnull = repmat(cols.null,nQuantiles,1) .* linspace(1,1.5,nQuantiles)';
% c = flipud(c);

alph = 0.05;

% null space
f = figure;
f.Position = [130   435   366   448];
ax(1) = subplot(2,1,1);
hold on;
for iq = 1:nQuantiles
    temp = squeeze(dat.sel.null(:,iq,:));
    temp = mySmooth(temp,sm,'reflect');
    mu = nanmean(temp,2);
    sig = nanstd(temp,[],2) ./ numel(obj);
    shadedErrorBar(obj(1).time,mu,sig,{'Color',cnull(iq,:)},alph,ax(1))
    %     plot(obj(1).time,mu,'Color',col(iq,:))
end
title('Null')
xlim(xlims)
xline(sample,'k--')
xline(delay,'k--')
xline(0,'k--')

clear cols
cols = getColors();
c = repmat(cols.potent,nQuantiles,1) ./ linspace(1,2,nQuantiles)';
cpotent = flipud(c);

% potent space
ax(2) = subplot(2,1,2);
hold on;
for iq = 1:nQuantiles
    temp = squeeze(dat.sel.potent(:,iq,:));
    temp = mySmooth(temp,sm,'zeropad');
    mu = nanmean(temp,2);
    sig = nanstd(temp,[],2) ./ numel(obj);
    shadedErrorBar(obj(1).time,mu,sig,{'Color',cpotent(iq,:)},alph,ax(2))
    %     plot(obj(1).time,mu,'Color',col(iq,:))
end
title('Potent')
xlim(xlims)
xline(sample,'k--')
xline(delay,'k--')
xline(0,'k--')



% for iq = 1:numel(ax)
%     ys(iq,:) = ax(iq).YLim;
% end
% yy(1) = min(min(ys));
% yy(2) = max(max(ys));
% for iq = 1:numel(ax)
%     ax(iq).YLim = yy;
% end



%%  mean rt bar graph for each session and quantile

% group by animal
[objix,uAnm] = groupSessionsByAnimal(meta);
nAnm = numel(uAnm);

anmRT = zeros(nQuantiles,nAnm);
for ianm = 1:nAnm
    anmRT(:,ianm) = nanmean(allrt(:,objix{ianm}),2);
end


c = gray;
ix = [1 75 150];
c = c(ix,:);

xs = [1 2 3];
f = figure;
ax = gca;
hold on;
for i = 1:numel(xs)
    temp = anmRT(i,:);
    b(i) = bar(xs(i),nanmean(temp));
    b(i).FaceColor = c(i,:);
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.4;

    vs(i) = scatter(randn(numel(temp),1) * 0.1 + xs(i)*ones(numel(temp),1),temp,5,'MarkerFaceColor','k',...
        'MarkerEdgeColor','k','LineWidth',1);

    e = errorbar(b(i).XEndPoints,nanmean(temp),nanstd(temp)./sqrt(numel(obj)),'LineStyle','none','Color','k','LineWidth',1);
    e.LineWidth = 0.5;
    e.CapSize = 2;
end
ylabel('Reaction time (s)')
ylim([0 ax.YLim(2)])
ax.XTick = xs;
ax.FontSize = 12;


%% plot mean epoch selectivity for each quantile

% group by animal
[objix,uAnm] = groupSessionsByAnimal(meta);
nAnm = numel(uAnm);

null = zeros(nQuantiles,numel(epochix),nAnm); % (quantile,epoch,session)
potent = zeros(nQuantiles,numel(epochix),nAnm); % (quantile,epoch,session)
for ianm = 1:nAnm
    null(:,:,ianm) = nanmean(dat.sel.epoch.null(:,:,objix{ianm}),3);
    potent(:,:,ianm) = nanmean(dat.sel.epoch.potent(:,:,objix{ianm}),3);
end

clear col cols clrs
cols = getColors();

if nQuantiles == 3
xs = {[1 1.5 2], [3 3.5 4], [5 5.5 6]};
elseif nQuantiles == 4
    xs = {[1 1.5 2 2.5], [3.5 4 4.5 5], [6 6.5 7 7.5]};
end

f = figure;
f.Position = [680   522   355   456];

ax = subplot(2,1,1);
hold on;
for iepoch = 1:size(null,2)
    for iq = 1:nQuantiles
        % null
        temp = squeeze(null(iq,iepoch,:));
        b(iq) = bar(xs{iepoch}(iq),nanmean(temp));
        b(iq).FaceColor = cnull(iq,:);
        b(iq).EdgeColor = 'none';
        b(iq).FaceAlpha = 1;
        b(iq).BarWidth = 0.4;

        vs(iq) = scatter(randn(numel(temp),1) * 0.001 + xs{iepoch}(iq)*ones(size(temp)),temp,5,'MarkerFaceColor','k',...
            'MarkerEdgeColor','k','LineWidth',1);

        e = errorbar(b(iq).XEndPoints,nanmean(temp),nanstd(temp)./sqrt(numel(obj)),'LineStyle','none','Color','k','LineWidth',1);
        e.LineWidth = 0.5;
        e.CapSize = 2;

    end

end
ylabel('Selectivity')
title('Null')
ylim([0 ax.YLim(2)])
ax.XTick = [];
ax.FontSize = 12;

ax = subplot(2,1,2);
hold on;
for iepoch = 1:size(potent,2)
    for iq = 1:nQuantiles
        % null
        temp = squeeze(potent(iq,iepoch,:));
        b(iq) = bar(xs{iepoch}(iq),nanmean(temp));
        b(iq).FaceColor = cpotent(iq,:);
        b(iq).EdgeColor = 'none';
        b(iq).FaceAlpha = 1;
        b(iq).BarWidth = 0.4;

        vs(iq) = scatter(randn(numel(temp),1) * 0.001 + xs{iepoch}(iq)*ones(size(temp)),temp,5,'MarkerFaceColor','k',...
            'MarkerEdgeColor','k','LineWidth',1);

        e = errorbar(b(iq).XEndPoints,nanmean(temp),nanstd(temp)./sqrt(numel(obj)),'LineStyle','none','Color','k','LineWidth',1);
        e.LineWidth = 0.5;
        e.CapSize = 2;

    end

end
ylabel('Selectivity')
title('Potent')
ylim([0 ax.YLim(2)])
ax.XTick = [];
ax.FontSize = 12;





%% t-test for epoch means


clear h p
for iepoch = 1:numel(epochix)

    % null
    temp = squeeze(dat.sel.epoch.null(:,iepoch,:)); % (quantile,session)

%     temp = squeeze(dat.sel.epoch.null_time(:,:,iepoch,:)); % (time,quantile,session)
%     temp = permute(temp,[ 1 3 2]);
%     temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));

    [h.null{iepoch}(1),p.null{iepoch}(1)] = ttest(temp(1,:),temp(2,:)); % 1,2
    [h.null{iepoch}(2),p.null{iepoch}(2)] = ttest(temp(2,:),temp(3,:)); % 2,3
    [h.null{iepoch}(3),p.null{iepoch}(3)] = ttest(temp(1,:),temp(3,:)); % 1,3

%     [p.null{iepoch}(1),h.null{iepoch}(1)] = ranksum(temp(1,:),temp(2,:)); % 1,2
%     [p.null{iepoch}(2),h.null{iepoch}(2)] = ranksum(temp(2,:),temp(3,:)); % 2,3
%     [p.null{iepoch}(3),h.null{iepoch}(3)] = ranksum(temp(1,:),temp(3,:)); % 1,3

%     p.null_anova{iepoch}(1) = anova1([temp(1,:) ; temp(2,:)]',[1 2],'off');
%     p.null_anova{iepoch}(2) = anova1([temp(2,:) ; temp(3,:)]',[1 2],'off');
%     p.null_anova{iepoch}(3) = anova1([temp(1,:) ; temp(3,:)]',[1 2],'off');


    % potent
    temp = squeeze(dat.sel.epoch.potent(:,iepoch,:));

%     temp = squeeze(dat.sel.epoch.potent_time(:,:,iepoch,:)); % (time,quantile,session)
%     temp = permute(temp,[ 1 3 2]);
%     temp = reshape(temp,size(temp,1)*size(temp,2),size(temp,3));

    [h.potent{iepoch}(1),p.potent{iepoch}(1)] = ttest(temp(1,:),temp(2,:)); % 1,2
    [h.potent{iepoch}(2),p.potent{iepoch}(2)] = ttest(temp(2,:),temp(3,:)); % 2,3
    [h.potent{iepoch}(3),p.potent{iepoch}(3)] = ttest(temp(1,:),temp(3,:)); % 1,3

%     [p.potent{iepoch}(1),h.potent{iepoch}(1)] = ranksum(temp(1,:),temp(2,:)); % 1,2
%     [p.potent{iepoch}(2),h.potent{iepoch}(2)] = ranksum(temp(2,:),temp(3,:)); % 2,3
%     [p.potent{iepoch}(3),h.potent{iepoch}(3)] = ranksum(temp(1,:),temp(3,:)); % 1,3

%     p.potent_anova{iepoch}(1) = anova1([temp(1,:) ; temp(2,:)]',[1 2],'off');
%     p.potent_anova{iepoch}(2) = anova1([temp(2,:) ; temp(3,:)]',[1 2],'off');
%     p.potent_anova{iepoch}(3) = anova1([temp(1,:) ; temp(3,:)]',[1 2],'off');
end


%% motion energy plot


c = gray;
ix = [1 75 150];
c = c(ix,:);

f = figure;
ax = gca;
hold on;
alph = 0.2;
for i = 1:nQuantiles
    temp = mySmooth(squeeze(allme(:,i,:)),31,'reflect');
    mu = mean(temp,2);
    sig = std(temp,[],2) ./ sqrt(numel(obj));
    shadedErrorBar(obj(1).time,mu,sig,{'Color',c(i,:),'LineWidth',1},alph,ax)
end
xline(sample,'k--')
xline(delay,'k--')
xline(0,'k--')
xlim(xlims)
xlabel('Time from go cue (s)')
ylabel('Motion Energy')


end
