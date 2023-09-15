% var exp for reconstructions in null vs potent of cd choice

%% find L/R selective cells per session
edges = [-0.8 0];
cond2use = [5 6];
for i = 1:numel(obj)
    cluix{i} = findSelectiveCells(obj(i),params(i),edges,cond2use);
    % cluix{i} = logical(ones(numel(params(i).cluid),1));
end

%% reconstruct single cell activity from cd choice in n/p

clear r2

cdix = 3;
cond2use = {'hit|miss'};
for sessix = 1:numel(meta)
    clear trialdat W proj 
    disp(['Session ' num2str(sessix) '/' num2str(numel(meta))])
    trix = findTrials(obj(sessix), cond2use);
    trix = trix{1};
    clus = find(cluix{sessix});


    trialdat.full = zscore_singleTrialNeuralData(obj(sessix));
    % trialdat.full = permute(obj(sessix).trialdat,[1 3 2]);
    trialdat.full = trialdat.full(:,trix,:); % (time,trials,neurons);

    W.null = cd_null(sessix).cd_mode_orth(:,cdix);
    W.potent = cd_potent(sessix).cd_mode_orth(:,cdix);
    fns = {'null','potent'};
    for j = 1:numel(fns)
        % single trials neural activity reconstructed from n/p
        trialdat.(fns{j}) = rez(sessix).recon.(fns{j})(:,trix,:); % (time,trials,dims)
        % project onto Wchoice
        proj.(fns{j}) = tensorprod(trialdat.(fns{j}),W.(fns{j}),3,1);
        % reconstruct n/p reconstructions from CD choice proj
        trialdat.recon.(fns{j}) = tensorprod(proj.(fns{j}),W.(fns{j}),3,2);

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

f = figure;
f.Position = [644   483   338   231];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

alln = [];
allp = [];
for sessix = 1:numel(meta)
    n = r2.null{sessix};
    alln = [alln n];
    p = r2.potent{sessix};
    allp = [allp p];
    % scatter(n,p,10,'filled','MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
end
% scatterhist(alln,allp)
scatter(alln,allp,10,'MarkerEdgeColor','k');%'MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('VE,null')
ylabel('VE,potent')

axis square
% r = refline(1,0);
% r.Color = 'k';

mdl = fitlm(alln,allp);
allr2 = mdl.Rsquared.Ordinary;
W = mdl.Coefficients.Estimate; % intercept, x1
xs = ax.XLim;
xs = linspace(xs(1),xs(2),100);
rr =  xs.* W(2) + W(1);
plot(xs,rr,'k--')
text(0.2,0.87,['r2=' num2str(allr2)])
text(0.2,0.77,['nCells=' num2str(numel(allp))])

ax = setLimits(ax,0.1);

%%

close all

alignment = (alln - allp) ./ (allp + alln);

cols = getColors;

f = figure;
f.Position = [644   483   338   231];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

h = histogram(alignment,30,'edgecolor','none','Normalization','count','Visible','off');
bars = h.Values;
binedges = h.BinEdges;

x = find(binedges<0);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.potent;

x = find(binedges>0);
x = x(1:end-1);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.null;

mu = mean(alignment);
sd = std(alignment);
n = numel(alignment);
nulldist = normrnd(0,sd,1,n);

%% hartigans dip test

[x,f] = ecdf(alignment);

[dip, p_value] = HartigansDipSignifTest(f,10000)
% [dip, p_value] = HartigansDipSignifTest(alignment,1000)






