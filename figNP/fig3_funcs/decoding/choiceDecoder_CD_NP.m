function acc = choiceDecoder_CD_NP(obj,rez,null,potent,params,cv,cond2use,cdix)

% decode choice from N/P activity
acc.null = nan(cv.numT-1,numel(rez));
acc.potent = nan(cv.numT-1,numel(rez));
for sessix = 1:numel(rez)
    disp(['Decoding session ' num2str(sessix) ' / ' num2str(numel(rez))])

    % trials

    trials_cond = params(sessix).trialid(cond2use);

    minHitTrials = cellfun(@(x) numel(x),trials_cond, 'UniformOutput',false);
    nhits = min(cell2mat(minHitTrials));

    trials_hit = cellfun(@(x) randsample(x,nhits), trials_cond, 'UniformOutput', false);
    trialsHit = cell2mat(trials_hit);
    trialsHit = trialsHit(:);

    trials.all = trialsHit;

    % labels (1 for right choice, 0 for left choice)
    Y = [ones(nhits,1) ; -ones(nhits,1)]; % right hits, left hits

    % DECODING
    disp('null')
    % -- null space
    % input
    X = null(sessix).trialdat(:,trials.all,cdix);

    % train/test split

    [trials.train,trials.trainidx] = datasample(trials.all,round(numel(trials.all)*cv.train),'Replace',false);
    trials.testidx = find(~ismember(trials.all,trials.train));
    trials.test = trials.all(trials.testidx);

    in.train.y = Y(trials.trainidx);
    in.test.y  = Y(trials.testidx);
    in.train.X = X(:,trials.trainidx,:);
    in.test.X  = X(:,trials.testidx,:);

    % decoding
    acc.null(:,sessix) = NP_ChoiceDecoder(in,cv,trials);

    % -- potent space
    disp('potent')
    % input
    X = potent(sessix).trialdat(:,trials.all,cdix);

    % train/test split

    [trials.train,trials.trainidx] = datasample(trials.all,round(numel(trials.all)*cv.train),'Replace',false);
    trials.testidx = find(~ismember(trials.all,trials.train));
    trials.test = trials.all(trials.testidx);

    in.train.y = Y(trials.trainidx);
    in.test.y  = Y(trials.testidx);
    in.train.X = X(:,trials.trainidx,:);
    in.test.X  = X(:,trials.testidx,:);

    % decoding
    acc.potent(:,sessix) = NP_ChoiceDecoder(in,cv,trials);


    % %     % shuffle labels for a 'null' distribution
    % %     Y = randsample(Y,numel(Y));
    % %
    % %     % train/test split
    % %
    % %     in.train.y = Y(trials.trainidx);
    % %     in.test.y  = Y(trials.testidx);
    % %
    % %     for ishuf = 1:rez.nShuffles
    % %         acc_shuf(:,sessix,ishuf,ifeat) = DLC_ChoiceDecoder(in,rez,trials);
    % %     end



end

% acc_shuf_ = reshape(acc_shuf,size(acc_shuf,1),size(acc_shuf,2)*size(acc_shuf,3));



%% plot

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);
trialStart = mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue);

cols = getColors;
lw = 2;
alph = 0.2;

fns = fieldnames(acc);
f = figure;
ax = gca;
hold on;
for i = 1:numel(fns)
    mu = nanmean(acc.(fns{i}),2);
    sigma = nanstd(acc.(fns{i}),[],2) ./ sqrt(numel(rez));
    shadedErrorBar(cv.tm(1:end-1), mu, sigma, {'Color',cols.(fns{i}), 'LineWidth', lw}, alph, ax);
end
xline(0,'k--','LineWidth',1)
xline(sample,'k--','LineWidth',1)
xline(delay,'k--','LineWidth',1)
ylim([ax.YLim(1) 1])
xlim([trialStart obj(1).time(end-10)])

xlabel('Time (s) from go cue')
ylabel([num2str(cv.nFolds) '-Fold CV Accuracy'])
title('Choice decoding from CD_NP','FontSize',8)

end


%% Helper functions

function acc = NP_ChoiceDecoder(in,cv,trials)

% bootstrap, and sessions happen outside this function

X_train = in.train.X;
y_train = in.train.y;
X_test = in.test.X;
y_test = in.test.y;


ix = 1;
for i = 1:cv.dt:(size(X_train,1)-cv.dt) % each timepoint

    ixs = i:(i+cv.dt-1);

    % train
    x_train = squeeze(nanmean(X_train(ixs,:,:),1));
    x_train = fillmissing(x_train,'constant',0);
    if size(x_train,1) == 1
        x_train = x_train';
    end

    %     mdl = fitcecoc(x_train,y_train);
    %     mdl = fitclinear(x_train,y_train);
    mdl = fitcsvm(x_train,y_train,'Standardize',false);
    cv_mdl = crossval(mdl,'KFold',cv.nFolds);

    pred = kfoldPredict(cv_mdl);

    acc(ix) = sum(pred == y_train) / numel(y_train);
    ix = ix + 1;
end


end