% rt = firstTongueRT(obj);
rt = firstJawRT(obj);

%%

close all
clear dat cols trialsBlock r2 temp trix acc


% params
opts.nFolds = 2; % number of CV folds
opts.binSize = 100; % ms
opts.dt = floor(opts.binSize / (params(1).dt*1000)); % samples
opts.tm = obj(1).time(1:opts.dt:numel(obj(1).time));
opts.numT = numel(opts.tm);

cond2use = 8:9; % right and left hits

nQuantiles = 4;
for sessix = 1:numel(meta)
    disp([num2str(sessix) '/' num2str(numel(meta))])
    for cdix = 1:numel(cd_potent(sessix).cd_labels)
        cd.potent = cd_potent(sessix).trialdat(:,:,cdix);
        cd.null = cd_null(sessix).trialdat(:,:,cdix);
        thisrt = rt{sessix};

        % right and left hit trials
        trix.all = cell2mat(params(sessix).trialid(cond2use)');

        X.potent = cd.potent(:,trix.all);
        X.null = cd.null(:,trix.all);
        thisrt = thisrt(trix.all)';

        quantiles = quantile(thisrt,nQuantiles-1); % returns n quantiles that divide the data set into evenly distributed n+1 segments.
        trialsBlock = getTrialsByRTBlock(thisrt,quantiles);
        y = zeros(size(trix.all));
        for i = 1:numel(trialsBlock)
            t = find(trialsBlock{i});
            y(t) = i;
        end

        Y = zeros(size(y,1),nQuantiles); % one-hot encoded labels (trials,rtBlock) -> 1st column is fastest rt-trial, last column is slowest rt-trial
        for i = 1:size(y,1)
            temp = zeros(nQuantiles,1);
            temp(y(i)) = 1;
            Y(i,:) = temp; 
        end

        % predict rt from CDs in each subspace
        ix = 1;
        for i = 1:opts.dt:(size(X.null,1)-opts.dt) % each timepoint
            ixs = i:(i+opts.dt-1);

            dat = nanmean(X.null(ixs,:),1)';

            CVMdl = fitcdiscr(dat,y,'KFold',opts.nFolds);
            pred = kfoldPredict(CVMdl);
            acc.null(ix,cdix,sessix) = sum(pred==y)/numel(y);

            dat = nanmean(X.potent(ixs,:),1)';

            CVMdl = fitcdiscr(dat,y,'KFold',opts.nFolds);
            pred = kfoldPredict(CVMdl);
            acc.potent(ix,cdix,sessix) = sum(pred==y)/numel(y);

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

    temp = squeeze(acc.null(:,i,:));
    shadedErrorBar(opts.tm(1:end-1),mean(temp,2),std(temp,[],2)./sqrt(numel(obj)),{'Color',cols.null,'LineWidth',2},alph,ax)
    temp = squeeze(acc.potent(:,i,:));
    shadedErrorBar(opts.tm(1:end-1),mean(temp,2),std(temp,[],2)./sqrt(numel(obj)),{'Color',cols.potent,'LineWidth',2},alph,ax)
    xline(0,'k--','LineWidth',1)
    xline(sample,'k--','LineWidth',1)
    xline(delay,'k--','LineWidth',1)
    xline(trialStart,'k--','LineWidth',1)
    xlim([trialStart, params(1).tmax-0.2])
%     ylim([ax.YLim(1) 1])

    xlabel('Time (s) from go cue')
    ylabel([num2str(opts.nFolds) '-Fold Accuracy'])
    title(labels{i},'FontSize',8)
end

