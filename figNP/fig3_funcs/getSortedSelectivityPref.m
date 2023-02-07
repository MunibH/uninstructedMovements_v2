function [plotsel , selDims] = getSortedSelectivityPref(obj,params,rez,edges,cond2use_ta,cond2use_st,modparams,space)
% Find cells that are significantly modulated by context (ranksum test, p-value = 0.01)
% As in Inagaki et al., Cell, 2022 ('A midbrain-thalamus...')
% modulatedCells = (1 x nDims) array where 1 means the dim is latedelay-selective and 0 means it is not

[~,ix1] = min(abs(obj(1).time - edges(1)));
[~,ix2] = min(abs(obj(1).time - edges(2)));

selectiveDims = [];
alldat = [];
allpref = [];
allve = [];
for sessix = 1:length(obj)                                           % For every session...

    if strcmpi(space,'null')
        trialdat = rez(sessix).N_null;                                 % Get the trial PSTH (time x trials x dims)
    elseif strcmpi(space,'potent')
        trialdat = rez(sessix).N_potent;
    end
    trialdat_mean = mean(trialdat(ix1:ix2,:,:),1);                 % Take the average FR for all cells during the presamp period
    temp = squeeze(trialdat_mean);                                % (trials x dims)

    if strcmpi(space,'null')
        dat = rez(sessix).N_null_psth(:,:,cond2use_ta);
        thisve = rez(sessix).ve.null;
    elseif strcmpi(space,'potent')
        dat = rez(sessix).N_potent_psth(:,:,cond2use_ta);
        thisve = rez(sessix).ve.potent;
    end
    
    if any(isnan(dat(:))) % skip session if no trials of one of the types -> EKH3 08-11 has no left misses
        continue
    end

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
    pref = zeros(1,nDims);
    for c = 1:nDims
        [pvals(c),hyp(c)] = ranksum(epochAvg{1}(:,c),epochAvg{2}(:,c),'alpha',sig);
        % get preferred direction for each dimension
        mus(1) = mean(epochAvg{1}(:,c));
        mus(2) = mean(epochAvg{2}(:,c));
        [~,pref(c)] = max(mus); % preferred direction is just direction with higher epochAvg
    end

    selectiveDims{sessix} = hyp;



    alldat = cat(2,alldat, dat);
    allpref = cat(2,allpref,pref);

    % vexp to scale by later
    allve = cat(2,allve, thisve);

end

% get nonpreferred directions for each dimension
allnonpref = allpref;
temp = allnonpref;
allnonpref(temp==1) = 2;
allnonpref(temp==2) = 1;

% find and remove dimensions that are all nans, unsure why?

% selectivity
for c = 1:size(alldat,2)
    selectivity(:,c) = alldat(:,c,allpref(c)) - alldat(:,c,allnonpref(c));     % (time,dims), trial-averaged activity
end

% scale by variance explained
selectivity = selectivity .* (allve.*100);


% Do normalization (as in Li et al, 2015 'A motor cortex circuit...')
maxSelect =  max(abs(selectivity),[],1);                    % Max magnitude of selectivity for each dim
% selectNorm = selectivity ./ maxSelect;
selectNorm = (selectivity);

selDims = logical(cell2mat(selectiveDims));

temp = mean(selectNorm(ix1:ix2,:),1,'omitnan');          % Find the average late delay sel for each dim
[~,sortix] = sort(temp,'descend');                          % Sort the selectivity  in descending order accordingly
sorted_select = selectNorm(:,sortix);
selmask = selDims(sortix);


plotsel = cat(2,sorted_select(:,selmask), sorted_select(:,~selmask)); % (time,(selective dims , nonselective dims))
selDims = 1:sum(selmask);

end