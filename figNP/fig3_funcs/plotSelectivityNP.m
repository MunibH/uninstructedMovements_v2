function plotSelectivityNP(meta,obj,params,rez,cond2use_ta,cond2use_st)

cols = getColors;
lw = 2;
alph = 0.2;

mainCond = 1; % right trials should always be higher

%%

% Sort by late-delay selectivity
% edges = [-0.6 0]; % (s) relative to go cue

% sort by late-sample selectivity
edges = [mode(obj(1).bp.ev.delay)-0.5 mode(obj(1).bp.ev.delay)] - 2.5;

% % sort by presample-selectivity
% edges = [mode(obj(1).bp.ev.sample)-0.3 mode(obj(1).bp.ev.sample)] - 2.5;

modparams.subTrials = 35;

[plotsel.null, plotsel.null_sep] = getSortedSelectivity(obj,params,rez,edges,cond2use_ta,cond2use_st,modparams,'null');
[plotsel.potent, plotsel.potent_sep] = getSortedSelectivity(obj,params,rez,edges,cond2use_ta,cond2use_st,modparams,'potent');

%% plot

xlims = [-2.3 2.5];

f = figure; f.Position = [794   516   302   384];
ax = gca;
hold on;
imagesc(obj(1).time, 1:size(plotsel.null,2),plotsel.null'); 
c = colorbar;
colormap(flipud(redblue))
xlabel('Time (s) from go cue')
ylabel('Dimensions')
c.Label.String = 'Selectivity (a.u.)';
title('null')
ax.YDir = 'reverse';
xlim(xlims)
ylim([1 size(plotsel.null,2)])
yline(plotsel.null_sep(1),'k--','LineWidth',2)
yline(plotsel.null_sep(2),'k--','LineWidth',2)


f = figure; f.Position = [794   516   302   384];
ax = gca;
hold on;
imagesc(obj(1).time, 1:size(plotsel.potent,2),plotsel.potent'); 
c = colorbar;
colormap(flipud(redblue))
xlabel('Time (s) from go cue')
ylabel('Dimensions')
c.Label.String = 'Selectivity (a.u.)';
title('potent')
ax.YDir = 'reverse';
xlim(xlims)
ylim([1 size(plotsel.potent,2)])
yline(plotsel.potent_sep(1),'k--','LineWidth',2)
yline(plotsel.potent_sep(2),'k--','LineWidth',2)



%% mean selectivity

col = getColors;
sm = 31;
lw = 2;
alph = 0.1;

sample = mode(obj(1).bp.ev.sample - 2.5);
delay = mode(obj(1).bp.ev.delay - 2.5);

null = mySmooth(plotsel.null(:,1:plotsel.null_sep(2)),sm,'reflect');
potent = mySmooth(plotsel.potent(:,1:plotsel.potent_sep(2)),sm,'reflect');

mu.null = mean(null,2);
mu.potent = mean(potent,2);

sig.null = std(null,[],2) ./ sqrt(size(null,2));
sig.potent = std(potent,[],2) ./ sqrt(size(potent,2));

f = figure; 
ax = gca;
hold on;
shadedErrorBar(obj(1).time,mu.null,sig.null,{'Color',col.null,'LineWidth',2,'LineStyle','-'},alph,ax);
shadedErrorBar(obj(1).time,mu.potent,sig.potent,{'Color',col.potent,'LineWidth',2,'LineStyle','-'},alph,ax);
xlim(xlims)
xline(0,'k--')
xline(delay,'k--')
xline(sample,'k--')

xlabel('Time (s) from go cue')
ylabel('Selectivity (a.u.)')


%% selectivity correlation matrix

sel = plotsel.null;

sel_corr_mat.null = zeros(size(sel,1),size(sel,1));
for i = 1:size(sel_corr_mat.null,1)
    for j = 1:size(sel_corr_mat.null,1)
        temp = corrcoef(sel(i,:),sel(j,:));
        sel_corr_mat.null(i,j) = temp(1,2);
    end
end

sel = plotsel.potent;

sel_corr_mat.potent = zeros(size(sel,1),size(sel,1));
for i = 1:size(sel_corr_mat.potent,1)
    for j = 1:size(sel_corr_mat.potent,1)
        temp = corrcoef(sel(i,:),sel(j,:));
        sel_corr_mat.potent(i,j) = temp(1,2);
    end
end



f = figure;
f.Position = [400   400   332   250];
ax = gca;
hold on;
imagesc(obj(1).time,obj(1).time,sel_corr_mat.null)
title('Null')
c = colorbar;
xlim(xlims)
ylim(xlims)
colormap(summer)
xline(0,'k--')
xline(delay,'k--')
xline(sample,'k--')
yline(0,'k--')
yline(delay,'k--')
yline(sample,'k--')

f = figure;
f.Position = [680   400   332   250];
ax = gca;
hold on;
imagesc(obj(1).time,obj(1).time,sel_corr_mat.potent)
title('Potent')
c = colorbar;
xlim(xlims)
ylim(xlims)
colormap(spring)
xline(0,'k--')
xline(delay,'k--')
xline(sample,'k--')
yline(0,'k--')
yline(delay,'k--')
yline(sample,'k--')



end


%% Helper functions

function [plotsel , sep] = getSortedSelectivity(obj,params,rez,edges,cond2use_ta,cond2use_st,modparams,space)
% Find cells that are significantly modulated by context (ranksum test, p-value = 0.01)
% As in Inagaki et al., Cell, 2022 ('A midbrain-thalamus...')
% modulatedCells = (1 x nDims) array where 1 means the dim is latedelay-selective and 0 means it is not

[~,ix1] = min(abs(obj(1).time - edges(1)));
[~,ix2] = min(abs(obj(1).time - edges(2)));

selectiveDims = [];
alldat = [];
for sessix = 1:length(obj)                                           % For every session...
    if strcmpi(space,'null')
        trialdat = rez(sessix).N_null;                                 % Get the trial PSTH (time x trials x dims)
    elseif strcmpi(space,'potent')
        trialdat = rez(sessix).N_potent;
    end
    trialdat_mean = mean(trialdat(ix1:ix2,:,:),1);                 % Take the average FR for all cells during the presamp period
    temp = squeeze(trialdat_mean);                                % (cells x trials)

    for cond = 1:numel(cond2use_st)
        trix = params(sessix).trialid{cond2use_st(cond)};
        trix2use = randsample(trix,modparams.subTrials,false);
        epochAvg{cond} = temp(trix,:); % (trials,dims)
    end


    % The p-value that you want to perform the ranksum test at
    sig = 0.01;

    nDims = size(epochAvg{1},2);                   % Get the number of dims
    pvals = zeros(1,nDims);                     % Store p-values for each cell
    hyp = zeros(1,nDims);                       % Store hyp test results for each cell
    for c = 1:nDims
        [pvals(c),hyp(c)] = ranksum(epochAvg{1}(:,c),epochAvg{2}(:,c),'alpha',sig);
    end

    selectiveDims{sessix} = hyp;

    if strcmpi(space,'null')
        dat = rez(sessix).N_null_psth(:,:,cond2use_ta);
    elseif strcmpi(space,'potent')
        dat = rez(sessix).N_potent_psth(:,:,cond2use_ta);
    end

    
    alldat = cat(2,alldat, dat);



end

% selectivity
selectivity = (alldat(:,:,1)-alldat(:,:,2));     % (time,dims), trial-averaged activity

% Do normalization (as in Li et al, 2015 'A motor cortex circuit...')
maxSelect =  max(abs(selectivity),[],1);                    % Max magnitude of selectivity for each dim
selectNorm = selectivity ./ maxSelect;

temp = mean(selectNorm(ix1:ix2,:),1,'omitnan');          % Find the average late delay sel for each dim
[~,sortix] = sort(temp,'descend');                          % Sort the selectivity  in descending order accordingly
sorted_select = selectNorm(:,sortix);

selmask = logical(cell2mat(selectiveDims)); % indices correspond to selectNorm
selmask = selmask(sortix); % corresponds to sorted_select
pctSelective = sum(selmask) / numel(selmask) * 100;

sep(1) = sum(selmask(1:floor(numel(selmask)/2))); % number of right-selective cells

sep(2) = sep(1) + sum(selmask(floor(numel(selmask)/2)+1:end)); % number of left-selective cells

plotsel = cat(2,sorted_select(:,selmask), sorted_select(:,~selmask)); % (time,(selective dims , nonselective dims))

end



















