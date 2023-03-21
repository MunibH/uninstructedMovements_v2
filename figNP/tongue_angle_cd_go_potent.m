clear cond dat tdat
cond2use = [2 3 4 5];
 
feat = 'tongue_angle';
% feat = 'bottomleft_tongue_ydisp_view2';
featix = find(ismember(kin(1).featLeg,feat));
for isess = 1:numel(meta)
    for icond = 1:numel(cond2use)
        trix = params(isess).trialid{cond2use(icond)};
        dat(:,isess,icond) = squeeze(nanmean(kin(isess).dat_std(:,trix,featix),2));
        tdat{isess}{icond} = kin(isess).dat(:,trix,featix);
    end
end


avgdat = squeeze(nanmean(dat,2));
c = getColors;

% figure;
% ax = gca;
% hold on;
% plot(obj(1).time,avgdat(:,1),'color',c.rhit)
% plot(obj(1).time,avgdat(:,2),'color',c.lhit)
% plot(obj(1).time,avgdat(:,3),'color',c.rmiss)
% plot(obj(1).time,avgdat(:,4),'color',c.lmiss)
% xlim([0 1])

figure;
ax = gca;
alph = 0.2;
lw = 1.5;
hold on;
sm = 0;
bc = 'zeropad';
temp = mySmooth(dat(:,:,1),sm,bc); mu = nanmean(temp,2); sd = nanstd(temp,[],2) ./ sqrt(numel(meta))*1.96;
shadedErrorBar(obj(1).time,mu,sd,{'Color',c.rhit,'LineWidth',lw},alph,ax)
temp = mySmooth(dat(:,:,2),sm,bc); mu = nanmean(temp,2); sd = nanstd(temp,[],2) ./ sqrt(numel(meta))*1.96;
shadedErrorBar(obj(1).time,mu,sd,{'Color',c.lhit,'LineWidth',lw},alph,ax)
temp = mySmooth(dat(:,:,3),sm,bc); mu = nanmean(temp,2); sd = nanstd(temp,[],2) ./ sqrt(numel(meta))*1.96;
shadedErrorBar(obj(1).time,mu,sd,{'Color',c.rmiss,'LineWidth',lw,'LineStyle','--'},alph,ax)
temp = mySmooth(dat(:,:,4),sm,bc); mu = nanmean(temp,2); sd = nanstd(temp,[],2) ./ sqrt(numel(meta))*1.96;
shadedErrorBar(obj(1).time,mu,sd,{'Color',c.lmiss,'LineWidth',lw,'LineStyle','--'},alph,ax)
xlim([0 1])
xlabel('Time from go cue (s)')
ylabel('Tongue angle (z-score)')

%%
% 
% f = figure;
for isess = 9% 1:numel(meta)

f = figure;
%     ax = nexttile;
    temp = tdat{isess};
    for i = 1:numel(cond2use)
        ys(i) = size(temp{i},2);
    end
    ys = cumsum(ys(1:3));
    temp = cell2mat(tdat{isess});
    imagesc(obj(1).time,1:size(temp,2),temp')
    colormap(linspecer)
    for i = 1:numel(ys)
        yline(ys(i),'w--')
    end
    xlabel('Time from go cue (s)')
    ylabel('Trials')
    colorbar
    xlim([0 0.5])
    title([meta(isess).anm ' ' meta(isess).date])
end

%%

f = figure;
% ax = gca;
% hold on;
for isess = 1:numel(meta)
    ax = nexttile; hold on;
    temp = tdat{isess};
    plot(obj(1).time,temp{1},'color',c.rhit)
    plot(obj(1).time,temp{2},'color',c.lhit)
    plot(obj(1).time,temp{3},'color',c.rmiss)
    plot(obj(1).time,temp{4},'color',c.lmiss)
end








