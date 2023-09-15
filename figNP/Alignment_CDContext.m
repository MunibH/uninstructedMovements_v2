% run this script after estimating subspaces and projecting neural data on
% to them (after running the main script)

%% find L/R selective cells per session
% only using cells with significant selectivity in this analysis

edges = [-0.8 0];
cond2use = [5 6];
for i = 1:numel(obj)
    cluix{i} = findSelectiveCells(obj(i),params(i),edges,cond2use);
end

%% reconstruct single cell activity from CDchoice projs in either null or potent space

cdix = 1; % index of cdchoice in the projs
cond2use = {'hit|miss'}; % specifying which trials 
for sessix = 1:numel(meta)
    clear trialdat W proj 
    disp(['Session ' num2str(sessix) '/' num2str(numel(meta))])
    trix = findTrials(obj(sessix), cond2use);
    trix = trix{1};
    clus = find(cluix{sessix});

    % get full single trial data (will compare reconstructed against this)
    % trialdat.full = zscore_singleTrialNeuralData(obj(sessix));
    % trialdat.full = permute(obj(sessix).trialdat,[1 3 2]);
    trialdat.full = trialdat.full(:,trix,:); % (time,trials,neurons);

    % get CDs
    W.null = cd_null(sessix).cd_mode_orth(:,cdix);
    W.potent = cd_potent(sessix).cd_mode_orth(:,cdix);
    fns = {'null','potent'};
    for j = 1:numel(fns)
        % single trials neural activity reconstructed from n/p
        trialdat.(fns{j}) = rez(sessix).recon.(fns{j})(:,trix,:); % (time,trials,dims)
        % project onto Wchoice
        proj.(fns{j}) = tensorprod(trialdat.(fns{j}),W.(fns{j}),3,1);
        % reconstruct data from CD choice proj
        trialdat.recon.(fns{j}) = tensorprod(proj.(fns{j}),W.(fns{j}),3,2);

        % for each cell, get R^2 b/w it's original data and reconstructed
        for k = 1:numel(clus) % for each cell
            thisclu = clus(k);
            % calculcate variance explained by CD choice
            orig = trialdat.full(:,:,thisclu); % (time,trials)

            fr = mean(mean(orig)); % subspace contribution method
            % weight = norm(W.(fns{j})(k));
            % r2.(fns{j}){sessix}(k) = fr*weight;

            recon = trialdat.recon.(fns{j})(:,:,thisclu); % (time,trials) % ve by recon method
            mdl = fitlm(orig(:),recon(:));
            r2.(fns{j}){sessix}(k) = mdl.Rsquared.Ordinary;
        end
    end
end


%% plot
close all

% concatenate R^2s from null and potent spaces into two vectors
alln = [];
allp = [];
for sessix = 1:numel(meta)
    n = r2.null{sessix};
    alln = [alln n];
    p = r2.potent{sessix};
    allp = [allp p];
    % scatter(n,p,10,'filled','MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
end
alignment = (alln - allp) ./ (allp + alln); % calculate alignment index

% histogram
f = figure;
f.Position = [644   483   338   231];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

h = histogram(alignment,30,'edgecolor','none','Normalization','count','Visible','off');






