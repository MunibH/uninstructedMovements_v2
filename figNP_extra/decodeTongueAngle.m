% decode tongue angle from np single trials / or go mode single trials (in
% np)

% data already loaded using st_elsayed_bootstrap_v2

%% single trial projs onto cd go

sm = 21;

for sessix = 1:numel(meta)
    cdix = 1;
    Wchoice = cd_null(sessix).cd_mode_orth(:,cdix);
    trialdat.null{sessix} = tensorprod(rez(sessix).recon.null,Wchoice,3,1); % (time,trials)
    trialdat.null{sessix} = mySmooth(trialdat.null{sessix},sm,'reflect');

    Wchoice = cd_potent(sessix).cd_mode_orth(:,cdix);
    trialdat.potent{sessix} = tensorprod(rez(sessix).recon.potent,Wchoice,3,1); % (time,trials)
    trialdat.potent{sessix} = mySmooth(trialdat.potent{sessix},sm,'reflect');
end

%% DECODING PARAMETERS

% input data = neural data (time*trials,neurons)
% output data = kin data   (time*trials,kin feats)

par.pre=6; % time bins prior to output used for decoding
par.post=0; % time bins after output used for decoding
par.dt = params(1).dt; % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using. Set params.dt and pre/post appropriately for you analysis
par.post_s = par.post .* params(1).dt;

% these parameters above are important for decoding accuracy. for actual
% analyses (like a final analysis to be included in a publication),
% you should vary these parameters only if you have a validation
% set that you won't test on until these parameters are set. Otherwise,
% there's a high risk of overfitting

% cross val folds
par.nFolds = 4;

% data sets
par.train = 1; % fraction of trials (just using cross-val here, matlab's kFoldPredict uses held out data for testing)
par.test = 1 - par.train;

% feature to decode using neural activity
par.feats = kin(1).featLeg;
feat = {'tongue_angle'};
par.featix = find(ismember(par.feats,feat));

% trials
par.cond2use = [5:8];

par.regularize = 1; % if 0, linear regression. if 1, ridge regression
par.regtype = 'ridge'; % ridge
par.learner = 'svm';

%% DECODING USING NP SINGLE TRIALS

close all

for sessix = 4%1:numel(meta) % 
    disp([num2str(sessix) ' / ' num2str(numel(meta))])

    % trials
    nTrials = cell2mat(cellfun(@(x) numel(x), params(sessix).trialid(par.cond2use) , 'UniformOutput',false));
    minTrials = min(nTrials(1:2));
    par.trials.train = [];
    for i = 1:numel(par.cond2use)
        trix = params(sessix).trialid{par.cond2use(i)};
        if i < 3
            k = min([50,minTrials]);
            % k = minTrials;
            % k = 50;
        else
            k = numel(trix);
        end
        par.trials.train = [par.trials.train ; randsample(trix,k,false)];
        par.trials.cond(i) = k;
    end

    % par.trials.train = cell2mat(params(sessix).trialid(par.cond2use)');

    % tongue angle
    Y.train = kin(sessix).dat(:,par.trials.train,par.featix); % (time,trials)
    Y.size = size(Y.train);
    Y.train = reshape(Y.train,Y.size(1)*Y.size(2),1);

    % NULL predictors
    X.train.null = rez(sessix).N_null(:,par.trials.train,:); % (time,trials,dims)
    X.size.null = size(X.train);
    X.train.null = reshape(X.train.null, size(X.train.null,1)*size(X.train.null,2),size(X.train.null,3));
    % reshape train and test data to account for prediction bin size
    X.train.null = reshapePredictors(X.train.null,par);
    % flatten inputs
    X.train.null = reshape(X.train.null,size(X.train.null,1),size(X.train.null,2)*size(X.train.null,3));

    % POTENT predictors
    X.train.potent = rez(sessix).N_potent(:,par.trials.train,:); % (time,trials,dims)
    X.size.potent = size(X.train);
    X.train.potent = reshape(X.train.potent, size(X.train.potent,1)*size(X.train.potent,2),size(X.train.potent,3));
    % reshape train and test data to account for prediction bin size
    X.train.potent = reshapePredictors(X.train.potent,par);
    % flatten inputs
    X.train.potent = reshape(X.train.potent,size(X.train.potent,1),size(X.train.potent,2)*size(X.train.potent,3));

    % FULL predictors
    X.train.full = permute(obj(sessix).trialdat(:,:,par.trials.train),[1 3 2]); % (time,trials,dims)
    X.size.full = size(X.train);
    X.train.full = reshape(X.train.full, size(X.train.full,1)*size(X.train.full,2),size(X.train.full,3));
    % reshape train and test data to account for prediction bin size
    X.train.full = reshapePredictors(X.train.full,par);
    % flatten inputs
    X.train.full = reshape(X.train.full,size(X.train.full,1),size(X.train.full,2)*size(X.train.full,3));

    X.train.null = mySmooth(X.train.null,21,'reflect');
    X.train.potent = mySmooth(X.train.potent,21,'reflect');

    % standardize data
    X.mu.null = mean(X.train.null,1,'omitnan');
    X.sigma.null = std(X.train.null,[],1,'omitnan');
    X.train.null = (X.train.null - X.mu.null) ./ X.sigma.null;

    X.mu.potent = mean(X.train.potent,1,'omitnan');
    X.sigma.potent = std(X.train.potent,[],1,'omitnan');
    X.train.potent = (X.train.potent - X.mu.potent) ./ X.sigma.potent;

    X.mu.full = mean(X.train.full,1,'omitnan');
    X.sigma.full = std(X.train.full,[],1,'omitnan');
    X.train.full = (X.train.full - X.mu.full) ./ X.sigma.full;
    
    Y.mu = mean(Y.train,1,'omitnan');
    Y.sigma = std(Y.train,[],1,'omitnan');
    Y.train = (Y.train - Y.mu);% ./ Y.sigma;

    % fill missing values in kinematics
    Y.nanmask = isnan(Y.train);
    Y.train = fillmissing(Y.train,'nearest'); % will replace 0s with nans after prediction
    X.train.null = fillmissing(X.train.null,'constant',0);
    X.train.potent = fillmissing(X.train.potent,'constant',0);
    X.train.full = fillmissing(X.train.full,'constant',0);

    if par.regularize
        mdl.null = fitrlinear(X.train.null,Y.train,'Learner',par.learner,'KFold',par.nFolds,'Regularization',par.regtype);
        mdl.potent = fitrlinear(X.train.potent,Y.train,'Learner',par.learner,'KFold',par.nFolds,'Regularization',par.regtype);
        mdl.full = fitrlinear(X.train.full,Y.train,'Learner',par.learner,'KFold',par.nFolds,'Regularization',par.regtype);
    else
        mdl.null = fitrlinear(X.train.null,Y.train,'Learner',par.learner,'KFold',par.nFolds);
        mdl.potent = fitrlinear(X.train.potent,Y.train,'Learner',par.learner,'KFold',par.nFolds);
        mdl.full = fitrlinear(X.train.full,Y.train,'Learner',par.learner,'KFold',par.nFolds);
    end

    pred.null = kfoldPredict(mdl.null);
    pred.potent = kfoldPredict(mdl.potent);
    pred.full = kfoldPredict(mdl.full);

    y = reshape(Y.train,Y.size(1),Y.size(2)); % original input data (centered)
    yhat.null = reshape(pred.null,Y.size(1),Y.size(2)); % prediction
    yhat.potent = reshape(pred.potent,Y.size(1),Y.size(2)); % prediction
    yhat.full = reshape(pred.full,Y.size(1),Y.size(2)); % prediction

    timeix = [0 1];
    ix = findTimeIX(obj(1).time,timeix);
    ix = ix(1):ix(2);

    rsq = corrcoef(y(ix,:),yhat.null(ix,:));
    r2.null(sessix) = rsq(1,2).^2;
    rsq = corrcoef(y(ix,:),yhat.potent(ix,:));
    r2.potent(sessix) = rsq(1,2).^2;
    rsq = corrcoef(y(ix,:),yhat.full(ix,:));
    r2.full(sessix) = rsq(1,2).^2;

    Y.train(Y.nanmask) = nan;
    pred.null(Y.nanmask) = nan;
    pred.potent(Y.nanmask) = nan;
    pred.full(Y.nanmask) = nan;

    y = reshape(Y.train,Y.size(1),Y.size(2)); % original input data (centered)
    yhat.null = reshape(pred.null,Y.size(1),Y.size(2)); % prediction
    yhat.potent = reshape(pred.potent,Y.size(1),Y.size(2)); % prediction
    yhat.full = reshape(pred.full,Y.size(1),Y.size(2)); % prediction

    %%
    xlims = [0 1];
    cmap = flipud(redblue);
    cmap(end,:) = 0;
    
    f = figure;
    f.Renderer = 'painters';
    f.Position = [557   599   817   276];

    ax1 = subplot(1,4,1); 
    temp = y;
    temp(isnan(temp)) = max(max(temp)) ;
    imagesc(obj(1).time, 1:size(y,2), temp')
    colormap(cmap); c1 = colorbar;
    clim(ax1,[-1 1])
    ax1 = prettifyPlot(ax1);
    yy = sum(par.trials.cond(1:2));
    yline(yy,'c--');

    ax2 = subplot(1,4,2); 
    temp = yhat.full;
    temp(isnan(temp)) = max(max(temp)) ;
    imagesc(obj(1).time, 1:size(y,2), temp')
    % title(num2str(r2.full(sessix)))
    colormap(cmap); c2 = colorbar;
    % clim(ax2,[c1.Limits])
    clim(ax2,[-1.5 1.5])
    ax2 = prettifyPlot(ax2);
    yy = sum(par.trials.cond(1:2));
    yline(yy,'c--');

    ax3 = subplot(1,4,3); 
    temp = yhat.null;
    temp(isnan(temp)) = max(max(temp)) +1;
    imagesc(obj(1).time, 1:size(y,2), temp')
    % title(num2str(r2.null(sessix)))
    colormap(cmap); c2 = colorbar;
    % clim(ax3,[c1.Limits])
    clim(ax3,[-2.3 2.3])
    ax3 = prettifyPlot(ax3);
    yy = sum(par.trials.cond(1:2));
    yline(yy,'c--');

    ax4 = subplot(1,4,4); 
    temp = yhat.potent;
    temp(isnan(temp)) = max(max(temp)) ;
    imagesc(obj(1).time, 1:size(y,2), temp')
    % title(num2str(r2.potent(sessix)))
    colormap(cmap); c3 = colorbar;
    % clim(ax4,[c1.Limits])
    clim(ax4,[-2 2])
    ax4 = prettifyPlot(ax4);
    yy = sum(par.trials.cond(1:2));
    yline(yy,'c--');

    % sgtitle([meta(sessix).anm ' ' meta(sessix).date])
    xlim(ax1,xlims)
    xlim(ax2,xlims)
    xlim(ax3,xlims)
    xlim(ax4,xlims)
    % 0.31 0.031 0.29
    %%

end



%% r2 bar plot

close all


f = figure;
f.Position = [748   394   202   272];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

dat = [r2.null' r2.potent'];

cols = getColors;
col{1} = cols.null;
col{2} = cols.potent;

xs = [1 2];
for i = 1:numel(xs)
    temp = dat(:,i);
    b(i) = bar(xs(i),mean(temp,'omitmissing'));
    b(i).FaceColor = col{i};
    b(i).EdgeColor = 'none';
    % b(i).FaceAlpha = 0.5;
    b(i).BarWidth = 0.5;
    vs(i) = scatter(xs(i)*ones(numel(temp),1),temp,15,'MarkerFaceColor','k',...
        'MarkerEdgeColor','w','XJitter','randn','XJitterWidth',0.001);
    errorbar(b(i).XEndPoints,nanmean(temp),nanstd(temp),'LineStyle','none','Color','k','LineWidth',1)
end
xx = [ones(numel(temp),1)*xs(1)  ones(numel(temp),1)*xs(2) ];
ll = line(xx',dat','color','k');
% ax.XTick = [];
ylabel("R^2")
ax = gca;
% ax.FontSize = 12;
xlim([xs(1)-0.5 xs(2)+0.5])



[h,p] = ttest(dat(:,1),dat(:,2))



%% DECODING USING CD

close all

for sessix = 11%1:numel(meta) % 11,12 example sessions
    disp([num2str(sessix) ' / ' num2str(numel(meta))])

    % trials
    par.trials.train = cell2mat(params(sessix).trialid(par.cond2use)');
    par.trials.cond = cell2mat(cellfun(@(x) numel(x) , params(sessix).trialid(par.cond2use), 'UniformOutput',false));

    % tongue angle
    Y.train = kin(sessix).dat(:,par.trials.train,par.featix); % (time,trials)
    Y.size = size(Y.train);
    Y.train = reshape(Y.train,Y.size(1)*Y.size(2),1);

    % NULL predictors
    X.train.null = trialdat.null{sessix}(:,par.trials.train); % (time,trials)
    X.size.null = size(X.train);
    X.train.null = reshape(X.train.null, size(X.train.null,1)*size(X.train.null,2),size(X.train.null,3));
    % reshape train and test data to account for prediction bin size
    X.train.null = reshapePredictors(X.train.null,par);
    % flatten inputs
    X.train.null = reshape(X.train.null,size(X.train.null,1),size(X.train.null,2)*size(X.train.null,3));

    % POTENT predictors
    X.train.potent = trialdat.potent{sessix}(:,par.trials.train); % (time,trials)
    X.size.potent = size(X.train);
    X.train.potent = reshape(X.train.potent, size(X.train.potent,1)*size(X.train.potent,2),size(X.train.potent,3));
    % reshape train and test data to account for prediction bin size
    X.train.potent = reshapePredictors(X.train.potent,par);
    % flatten inputs
    X.train.potent = reshape(X.train.potent,size(X.train.potent,1),size(X.train.potent,2)*size(X.train.potent,3));

    % standardize data
    X.mu.null = mean(X.train.null,1,'omitnan');
    X.sigma.null = std(X.train.null,[],1,'omitnan');
    X.train.null = (X.train.null - X.mu.null) ./ X.sigma.null;

    X.mu.potent = mean(X.train.potent,1,'omitnan');
    X.sigma.potent = std(X.train.potent,[],1,'omitnan');
    X.train.potent = (X.train.potent - X.mu.potent) ./ X.sigma.potent;
    
    Y.mu = mean(Y.train,1,'omitnan');
    Y.sigma = std(Y.train,[],1,'omitnan');
    Y.train = (Y.train - Y.mu);% ./ Y.sigma;

    % fill missing values in kinematics
    Y.nanmask = isnan(Y.train);
    Y.train = fillmissing(Y.train,'nearest'); % will replace 0s with nans after prediction
    X.train.null = fillmissing(X.train.null,'constant',0);
    X.train.potent = fillmissing(X.train.potent,'constant',0);

    if par.regularize
        mdl.null = fitrlinear(X.train.null,Y.train,'Learner','leastsquares','KFold',par.nFolds,'Regularization',par.regtype);
        mdl.potent = fitrlinear(X.train.potent,Y.train,'Learner','leastsquares','KFold',par.nFolds,'Regularization',par.regtype);
    else
        mdl.null = fitrlinear(X.train.null,Y.train,'Learner','leastsquares','KFold',par.nFolds);
        mdl.potent = fitrlinear(X.train.potent,Y.train,'Learner','leastsquares','KFold',par.nFolds);
    end

    pred.null = kfoldPredict(mdl.null);
    pred.potent = kfoldPredict(mdl.potent);

    y = reshape(Y.train,Y.size(1),Y.size(2)); % original input data (centered)
    yhat.null = reshape(pred.null,Y.size(1),Y.size(2)); % prediction
    yhat.potent = reshape(pred.potent,Y.size(1),Y.size(2)); % prediction

    timeix = [0 1];
    ix = findTimeIX(obj(1).time,timeix);
    ix = ix(1):ix(2);

    rsq = corrcoef(y(ix,:),yhat.null(ix,:));
    r2.null(sessix) = rsq(1,2).^2;
    rsq = corrcoef(y(ix,:),yhat.potent(ix,:));
    r2.potent(sessix) = rsq(1,2).^2;

    Y.train(Y.nanmask) = nan;
    pred.null(Y.nanmask) = nan;
    pred.potent(Y.nanmask) = nan;

    y = reshape(Y.train,Y.size(1),Y.size(2)); % original input data (centered)
    yhat.null = reshape(pred.null,Y.size(1),Y.size(2)); % prediction
    yhat.potent = reshape(pred.potent,Y.size(1),Y.size(2)); % prediction

    %%
    xlims = [0 1];
    cmap = flipud(linspecer);
    
    f = figure;
    f.Renderer = 'painters';
    f.Position = [557   599   817   276];

    ax1 = subplot(1,4,1); 
    temp = y;
    temp(isnan(temp)) = max(max(temp)) ;
    imagesc(obj(1).time, 1:size(y,2), temp')
    colormap(cmap); c1 = colorbar;
    ax1 = prettifyPlot(ax1);
    yy = sum(par.trials.cond(1:2));
    yline(yy,'w--');

    ax2 = subplot(1,4,2); 
    temp = y;
    temp(isnan(temp)) = max(max(temp)) ;
    imagesc(obj(1).time, 1:size(y,2), temp')
    % title(num2str(r2.null(sessix)))
    colormap(cmap); c2 = colorbar;
    % clim(ax2,[c1.Limits])
    % clim(ax2,[-2 2])
    ax2 = prettifyPlot(ax2);
    yy = sum(par.trials.cond(1:2));
    yline(yy,'w--');

    ax3 = subplot(1,4,3); 
    temp = yhat.null;
    temp(isnan(temp)) = max(max(temp)) ;
    imagesc(obj(1).time, 1:size(y,2), temp')
    % title(num2str(r2.null(sessix)))
    colormap(cmap); c2 = colorbar;
    % clim(ax3,[c1.Limits])
    % clim(ax3,[-2 2])
    ax3 = prettifyPlot(ax3);
    yy = sum(par.trials.cond(1:2));
    yline(yy,'w--');

    ax4 = subplot(1,4,4); 
    temp = yhat.potent;
    temp(isnan(temp)) = max(max(temp)) ;
    imagesc(obj(1).time, 1:size(y,2), temp')
    % title(num2str(r2.potent(sessix)))
    colormap(cmap); c3 = colorbar;
    % clim(ax4,[c1.Limits])
    % clim(ax4,[-2 2])
    ax4 = prettifyPlot(ax4);
    yy = sum(par.trials.cond(1:2));
    yline(yy,'w--');

    sgtitle([meta(sessix).anm ' ' meta(sessix).date])
    xlim(ax1,xlims)
    xlim(ax2,xlims)
    xlim(ax3,xlims)
    xlim(ax4,xlims)
    % 0.21 0.43
    %%
end



%% r2 bar plot

close all


f = figure;
f.Position = [748   414   249   252];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

dat = [r2.null' r2.potent'];

cols = getColors;
col{1} = cols.null;
col{2} = cols.potent;

xs = [1 2];
for i = 1:numel(xs)
    temp = dat(:,i);
    b(i) = bar(xs(i),mean(temp,'omitmissing'));
    b(i).FaceColor = col{i};
    b(i).EdgeColor = 'none';
    % b(i).FaceAlpha = 0.5;
    b(i).BarWidth = 0.7;
    vs(i) = scatter(xs(i)*ones(numel(temp),1),temp,20,'MarkerFaceColor','k',...
        'MarkerEdgeColor','w','XJitter','randn','XJitterWidth',0.25);
    errorbar(b(i).XEndPoints,nanmean(temp),nanstd(temp)./sqrt(numel(meta)),'LineStyle','none','Color','k','LineWidth',1)
end
% ax.XTick = [];
ylabel("R^2")
ax = gca;
% ax.FontSize = 12;




[h,p] = ttest(dat(:,1),dat(:,2))




























