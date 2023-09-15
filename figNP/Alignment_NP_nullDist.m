% null distributinos for alignment indicies
% see Elsayed 2016 supplements
% (https://static-content.springer.com/esm/art%3A10.1038%2Fncomms13239/MediaObjects/41467_2016_BFncomms13239_MOESM852_ESM.pdf)

iters = 50;

r2.null = [];
r2.potent = [];
for i = 1:iters
    disp(['Iteration ' num2str(i) '/' num2str(iters)])

    % choose a random session

    sessix = randsample(numel(meta),1);

    % get covariance matrices
    Cnull = rez(sessix).covNull;
    Cpot  = rez(sessix).covPotent;
    % Cnull = rez(sessix).Cprep;
    % Cpot  = rez(sessix).Cmove;

    % eigenvecs and vals of C
    [Vnull,Dnull] = eig(Cnull);
    [Vpot,Dpot] = eig(Cpot);

    % sample a random vector (each element drawn from N(0,1))
    vnull = randn(size(Cnull,1),rez(sessix).dPrep);
    vpot = randn(size(Cpot,1),rez(sessix).dMove);

    % sample random vector aligned to covariance structures as in Elsayed
    % 2016
    a = Vnull * Dnull * vnull;
    b = a ./ norm(a,2);

    % get orthonormal basis of b, defined by left singular vecs
    % can just use matlab's orth() function, since it uses svd's 'U' to
    % obtain orth basis
    valign.null = orth(b);

    % repeat for potent
    a = Vpot * Dpot * vpot;
    b = a ./ norm(a,2);
    valign.pot = orth(b);

    clear trialdat cluix trix clus
    % find L/R selective cells per session
    edges = [-0.8 0];
    cond2use = [2 3]; %[5 6];
    cluix = 1:numel(params(sessix).cluid);
    % cluix = findSelectiveCells(obj(sessix),params(sessix),edges,cond2use);

    % reconstruct single cell activity from n/p
    cond2use = {'hit|miss'};

    trix = findTrials(obj(sessix), cond2use);
    trix = trix{1};
    clus = find(cluix);


    % trialdat.full = zscore_singleTrialNeuralData(obj(sessix));
    % trialdat.full = trialdat.full(:,trix,:); % (time,trials,neurons);

    trialdat.full = permute(obj(sessix).trialdat,[1 3 2]);
    trialdat.full = trialdat.full(:,trix,:);

    % project trialdat onto valign
    trialdat.null = tensorprod(trialdat.full,valign.null,3,1);
    trialdat.potent = tensorprod(trialdat.full,valign.pot,3,1);

    % reconstruct full
    rez(sessix).recon.null = tensorprod(trialdat.null,valign.null,3,2);
    rez(sessix).recon.potent = tensorprod(trialdat.potent,valign.pot,3,2);

    fns = {'null','potent'};
    for j = 1:numel(fns)
        % single trials neural activity reconstructed from n/p
        % trialdat.recon.(fns{j}) = rez(sessix).recon.(fns{j})(:,trix,:); % (time,trials,dims)
        trialdat.recon.(fns{j}) = rez(sessix).recon.(fns{j}); % (time,trials,dims)

        % % % getting r2 between all neurons at once, rather than
        % % % individually
        % % x = trialdat.recon.(fns{j});
        % % y = trialdat.full;
        % % mdl = fitlm(x(:),y(:));
        % % r2.(fns{j}) = [r2.(fns{j}) ; mdl.Rsquared.Ordinary];

        for k = 1:numel(clus) % for each cell
            thisclu = clus(k);
            % calculcate variance explained by CD choice
            orig = trialdat.full(:,:,thisclu); % (time,trials)

            fr = mean(mean(orig)); % subspace contribution method
            % weight = norm(W.(fns{j})(k));
            % r2.(fns{j}){sessix}(k) = fr*weight;

            recon = trialdat.recon.(fns{j})(:,:,thisclu); % (time,trials) % ve by recon method
            mdl = fitlm(orig(:),recon(:));
            r2.(fns{j}) = [r2.(fns{j}) ; mdl.Rsquared.Ordinary];
        end
    end
end



%%

close all

% r2.potent = r2.potent(1:numel(r2.null));

alignment = (r2.null - r2.potent) ./ (r2.potent + r2.null);

cols = getColors;

f = figure;
f.Position = [644   483   338   231];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

% h = histogram(alignment,30,'edgecolor','none','Normalization','count','facecolor',[0.3 0.3 0.3]);
h = histogram(alignment,30,'edgecolor','none','Normalization','probability','facecolor',[0.3 0.3 0.3]);
% bars = h.Values;
% binedges = h.BinEdges;
%
% x = find(binedges<0);
% b = bar(binedges(x),bars(x));
% b.BarWidth = 1;
% b.EdgeColor = 'none';
% b.FaceColor = cols.potent;
%
% x = find(binedges>0);
% x = x(1:end-1);
% b = bar(binedges(x),bars(x));
% b.BarWidth = 1;
% b.EdgeColor = 'none';
% b.FaceColor = cols.null;




















