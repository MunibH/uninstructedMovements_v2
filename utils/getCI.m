function CI = getCI(ts)



%%


% CI for a time series (time,dim2)
SEM = nanstd(ts, [], 2) ./ sqrt(size(ts,2));               % Standard Error
t = tinv([0.025  0.975],size(ts,2)-1);      % T-Score
% CI = 1.96 * SEM;
CI = t(1) * SEM;
% CI(:,1) = mean(ts,2) + t(2)*SEM;                      % Confidence Intervals
% CI(:,2) = mean(ts,2) + t(1)*SEM;                      % Confidence Intervals

% f = figure;
% ax = gca;
% shadedErrorBar(1:(size(ts,1)),mean(ts,2),CI',{'Color','k','LineWidth',2},0.15,ax)

end