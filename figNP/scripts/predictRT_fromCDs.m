% rt = firstTongueRT(obj);
rt = firstJawRT(obj);

%%

close all
clear dat cols trialsBlock r2


% params
opts.nFolds = 2; % number of iterations (bootstrap)

opts.binSize = 100; % ms
opts.dt = floor(opts.binSize / (params(1).dt*1000)); % samples
opts.tm = obj(1).time(1:opts.dt:numel(obj(1).time));
opts.numT = numel(opts.tm);

for sessix = 1:numel(meta)
    disp([num2str(sessix) '/' num2str(numel(meta))])
    for cdix = 1:numel(cd_potent(sessix).cd_labels)
        cd.potent = cd_potent(sessix).trialdat(:,:,cdix);
        cd.null = cd_null(sessix).trialdat(:,:,cdix);
        thisrt = rt{sessix};

        % right and left hit trials
        trix.all = cell2mat(params(sessix).trialid(2:3)');

        X.potent = cd.potent(:,trix.all);
        X.null = cd.null(:,trix.all);
        Y = thisrt(trix.all)';
        Y = zscore(Y);

        % predict rt from CDs in each subspace
        ix = 1;
        for i = 1:opts.dt:(size(X.null,1)-opts.dt) % each timepoint
            ixs = i:(i+opts.dt-1);

            dat = nanmean(X.null(ixs,:),1)';
            dat = cat(2,dat,ones(size(dat,1),1));

            CVMdl = fitrlinear(dat',Y,'KFold',opts.nFolds,'ObservationsIn','columns');
            pred = kfoldPredict(CVMdl);
            r2_ = corrcoef(Y,pred);
            r2.null(ix,cdix,sessix) = r2_(1,2).*2;

            dat = nanmean(X.potent(ixs,:),1)';
            dat = cat(2,dat,ones(size(dat,1),1));

            CVMdl = fitrlinear(dat',Y,'KFold',opts.nFolds,'ObservationsIn','columns');
            pred = kfoldPredict(CVMdl);
            r2_ = corrcoef(Y,pred);
            r2.potent(ix,cdix,sessix) = r2_(1,2).*2;

            ix = ix + 1;
        end

    end

end

%%

cols = getColors();

alph = 0.2;

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);
trialStart = mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue);

f = figure;
labels = cd_potent(1).cd_labels;
nCDs = numel(labels);
for i = 1:nCDs
    ax = subplot(nCDs,1,i);
    hold on;

    acc = squeeze(r2.null(:,i,:));
    shadedErrorBar(opts.tm(1:end-1),mean(acc,2),std(acc,[],2)./sqrt(numel(obj)),{'Color',cols.null,'LineWidth',2},alph,ax)
    acc = squeeze(r2.potent(:,i,:));
    shadedErrorBar(opts.tm(1:end-1),mean(acc,2),std(acc,[],2)./sqrt(numel(obj)),{'Color',cols.potent,'LineWidth',2},alph,ax)
    xline(0,'k--','LineWidth',1)
    xline(sample,'k--','LineWidth',1)
    xline(delay,'k--','LineWidth',1)
    xline(trialStart,'k--','LineWidth',1)
    xlim([trialStart, params(1).tmax-0.2])
    ylim([ax.YLim(1) 1])

    xlabel('Time (s) from go cue')
    ylabel([num2str(opts.nFolds) '-Fold R2'])
    title(labels{i},'FontSize',8)
end

