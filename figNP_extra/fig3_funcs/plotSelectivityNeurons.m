function plotSelectivityNeurons(meta,obj,params,cond2use_ta,cond2use_st)

cols = getColors;
lw = 2;
alph = 0.2;

mainCond = 1; % right trials should always be higher

%%

% % Sort by delay selectivity
% edges = [mode(obj(1).bp.ev.delay) mode(obj(1).bp.ev.goCue)] - 2.5;
% edges = [-0.8 0]; % (s) relative to go cue

% % sort by late-sample selectivity
% edges = [mode(obj(1).bp.ev.delay)-0.5 mode(obj(1).bp.ev.delay)] - 2.5;

% % sort by presample-selectivity
% edges = [mode(obj(1).bp.ev.sample)-0.3 mode(obj(1).bp.ev.sample)] - 2.5;

% % sort by response-selectivity
edges = [mode(obj(1).bp.ev.goCue)+0.02 mode(obj(1).bp.ev.goCue)+0.5] - 2.5;

modparams.subTrials = 35;

[plotsel.dat, plotsel.sep, plotsel.sel_dims] = getSortedSelectivity(obj,params,edges,cond2use_ta,cond2use_st,modparams);

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


end




function [plotsel , sep, sel_dims] = getSortedSelectivity(obj,params,edges,cond2use_ta,cond2use_st,modparams)
% Find cells that are significantly modulated by context (ranksum test, p-value = 0.01)
% As in Inagaki et al., Cell, 2022 ('A midbrain-thalamus...')
% modulatedCells = (1 x nDims) array where 1 means the dim is latedelay-selective and 0 means it is not

[~,ix1] = min(abs(obj(1).time - edges(1)));
[~,ix2] = min(abs(obj(1).time - edges(2)));

selectiveDims = [];
alldat = [];
for sessix = 1:length(obj)                                           % For every session...
    trialdat = obj(sessix).trialdat;

    trialdat_mean = mean(trialdat(ix1:ix2,:,:),1);                 % Take the average FR for all cells during the presamp period
    temp = squeeze(trialdat_mean);                                % (cells x trials)

    for cond = 1:numel(cond2use_st)
        trix = params(sessix).trialid{cond2use_st(cond)};
        trix2use = randsample(trix,modparams.subTrials,false);
        epochAvg{cond} = temp(:,trix); % (cells,trials)
    end


    % The p-value that you want to perform the ranksum test at
    sig = 0.01;

    nCells = size(epochAvg{1},1);                   % Get the number of dims
    pvals = zeros(1,nCells);                     % Store p-values for each cell
    hyp = zeros(1,nCells);                       % Store hyp test results for each cell
    for c = 1:nCells
        [pvals(c),hyp(c)] = ranksum(epochAvg{1}(c,:),epochAvg{2}(c,:),'alpha',sig);
    end

    selectiveDims{sessix} = hyp;

    dat = obj(sessix).psth(:,:,cond2use_ta);


    alldat = cat(2,alldat, dat);

    clear pvals hyp

end

% selectivity
selectivity = (alldat(:,:,1)-alldat(:,:,2));     % (time,dims), trial-averaged activity

% Do normalization (as in Li et al, 2015 'A motor cortex circuit...')
maxSelect =  max(abs(selectivity),[],1);                    % Max magnitude of selectivity for each dim
selectNorm = selectivity ./ maxSelect;

temp = mean(selectNorm(ix1:ix2,:),1,'omitnan');          % Find the average late delay sel for each dim
[~,sortix] = sort(temp,'descend');                          % Sort the selectivity  in descending order accordingly
sorted_select = selectNorm(:,sortix);

sel_dims = cell2mat(selectiveDims);

selmask = logical(sel_dims); % indices correspond to selectNorm
selmask = selmask(sortix); % corresponds to sorted_select
pctSelective = sum(selmask) / numel(selmask) * 100;

sep(1) = sum(selmask(1:floor(numel(selmask)/2))); % number of right-selective cells

sep(2) = sep(1) + sum(selmask(floor(numel(selmask)/2)+1:end)); % number of left-selective cells

plotsel = cat(2,sorted_select(:,selmask), sorted_select(:,~selmask)); % (time,(selective dims , nonselective dims))

end







