function [plotsel , sep, sel_dims] = getSortedSelectivity(obj,params,rez,edges,cond2use_ta,cond2use_st,modparams,space)
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

sel_dims = cell2mat(selectiveDims);

selmask = logical(sel_dims); % indices correspond to selectNorm
selmask = selmask(sortix); % corresponds to sorted_select
pctSelective = sum(selmask) / numel(selmask) * 100;

sep(1) = sum(selmask(1:floor(numel(selmask)/2))); % number of right-selective cells

sep(2) = sep(1) + sum(selmask(floor(numel(selmask)/2)+1:end)); % number of left-selective cells

plotsel = cat(2,sorted_select(:,selmask), sorted_select(:,~selmask)); % (time,(selective dims , nonselective dims))

end