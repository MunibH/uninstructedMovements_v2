%% number of sessions and total single units

nSessions = numel(meta);
disp(['# sessions = ' num2str(nSessions)])

nSingleCells = 0;
for i = 1:numel(meta)
    nSingleCells = nSingleCells + numel(params(i).cluid);
end
disp(['# single units = ' num2str(nSingleCells)])

%% sample selectivity (DR L/R)

edges = [-1.5  -0.9];
cond2use = [8 9];
nSampleSelective = 0;
for i = 1:numel(meta)
    cluix{i} = findSelectiveCells(obj(i),params(i),edges,cond2use);
    cluix{i}(isnan(cluix{i})) = 0;
    nSampleSelective = nSampleSelective + sum(cluix{i});
end
disp(['% sample selective = ' num2str(nSampleSelective/nSingleCells)])



%% delay selectivity

edges = [-0.8 0];
cond2use = [8 9];
nDelaySelective = 0;
for i = 1:numel(meta)
    cluix{i} = findSelectiveCells(obj(i),params(i),edges,cond2use);
    cluix{i}(isnan(cluix{i})) = 0;
    nDelaySelective = nDelaySelective + sum(cluix{i});
end
disp(['% delay selective = ' num2str(nDelaySelective/nSingleCells)])

%% response selectivity

edges = [0 0.5];
cond2use = [8 9];
nRespSelective = 0;
for i = 1:numel(meta)
    cluix{i} = findSelectiveCells(obj(i),params(i),edges,cond2use);
    cluix{i}(isnan(cluix{i})) = 0;
    nRespSelective = nRespSelective + sum(cluix{i});
end
disp(['% response selective = ' num2str(nRespSelective/nSingleCells)])





