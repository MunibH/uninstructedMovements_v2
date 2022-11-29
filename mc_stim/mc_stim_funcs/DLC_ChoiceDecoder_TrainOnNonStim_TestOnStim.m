function [acc_ctrl, acc_stim] = DLC_ChoiceDecoder_TrainOnNonStim_TestOnStim(in,rez)

% bootstrap, and sessions happen outside this function

X_train = in.train.X;
y_train = in.train.y;
X_test = in.test.X;
y_test = in.test.y;


%%

ix = 1;
for i = 1:rez.dt:(size(X_train,1)-rez.dt) % each timepoint

    ixs = i:(i+rez.dt-1);

    % train
    x_train = squeeze(nanmean(X_train(ixs,:,:),1));
    x_train = fillmissing(x_train,'constant',0);
    if size(x_train,1) == 1
        x_train = x_train';
    end

%     mdl = fitcecoc(x_train,y_train);
%     mdl = fitclinear(x_train,y_train);
    mdl = fitcsvm(x_train,y_train,'Standardize',true);
    cv_mdl = crossval(mdl,'KFold',rez.nFolds);
    pred = kfoldPredict(cv_mdl); % predict on heldout data from crossvalidated training
    acc_ctrl(ix) = sum(pred == y_train) / numel(y_train);

    % test
    x_test = squeeze(nanmean(X_test(ixs,:,:),1));
        x_test = fillmissing(x_test,'constant',0);
    if size(x_test,1) == 1
        x_test = x_test';
    end
    for iFold = 1:rez.nFolds
        test_pred = predict(cv_mdl.Trained{iFold},x_test);
        acc_stim(ix,iFold) = sum(test_pred == y_test) / numel(y_test);
    end

    ix = ix + 1;
end
acc_stim = mean(acc_stim,2); % mean over folds

end