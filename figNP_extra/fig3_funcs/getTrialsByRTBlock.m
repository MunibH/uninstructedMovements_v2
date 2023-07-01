function trialsBlock = getTrialsByRTBlock(rt,quantiles)
for i = 1:numel(quantiles)
    if i == 1
        trialsBlock{i} = rt<=quantiles(i);
%     elseif i == numel(quantiles)
%         trialsBlock{i} = rt>=quantiles(i);
    else
        trialsBlock{i} = rt>quantiles(i-1) & rt<=quantiles(i);
    end
    trialsBlock{i+1} = rt >= quantiles(i);
end