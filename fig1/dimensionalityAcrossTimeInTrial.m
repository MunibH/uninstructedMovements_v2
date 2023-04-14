
close all
% pca

% psth = obj.psth(:,:,1); 
% % psth = psth - mean(psth,2);
% psth = normalize(psth);
% [coeff, score, eig, tsquared, explained] = pca(psth);

trialdat = permute(obj.trialdat,[1 3 2]);
dims = size(trialdat);
trialdat = reshape(trialdat,size(trialdat,1)*size(trialdat,2),size(trialdat,3));
trialdat = mySmooth(trialdat,21);
trialdat = normalize(trialdat);
[coeff, score, eig, tsquared, explained] = pca(trialdat);
score = reshape(score,dims(1),dims(2),dims(3));
score = squeeze(mean(score,2));

f = figure;
for i = 1:size(psth,2)
    ax = nexttile;
    plot(obj.time,score(:,i),'k','LineWidth',2)
    title(num2str(explained(i)))
end
ax = nexttile;
hold on;
plot(1:size(psth,2),cumsum(explained),'b.','MarkerSize',10)
plot(1:size(psth,2),cumsum(explained),'b-')
d = numComponentsToExplainVariance(explained,90);
yline(90,'k--')
xline(d,'k--')
title(['d = ' num2str(d)])
% dimensionality is 8

% estimate dimensionailty using 100 ms bins

binSize = 0.1; % 100 ms window
nSamp = binSize / params.dt; 

% for each 100 ms window, concatenate data across trials -> (time*trials,neurons)
clear d_sub d_sub_trial
edges = 1:nSamp:numel(obj.time);
for i = 1:numel(edges)
    clear tix
    if edges(i)+nSamp > numel(obj.time)
        tix = edges(i) : ( edges(i) + (numel(obj.time) - edges(i)) );
    else
        tix = edges(i) : ( edges(i) + nSamp );
    end
    dat = psth(tix,:);
    [coeff, score, eig, tsquared, explained] = pca(dat);
    d_sub(i) = numComponentsToExplainVariance(explained,90);

    dat = permute(obj.trialdat(tix,:,:),[1 3 2]);
    dat = reshape(dat,size(dat,1)*size(dat,2),size(dat,3));
    dat = normalize(dat);
    [coeff, score, eig, tsquared, explained] = pca(dat);
    d_sub_trial(i) = numComponentsToExplainVariance(explained,90);

end



ax = nexttile;
hold on;
plot(obj.time(edges),d_sub,'b.','MarkerSize',10)
plot(obj.time(edges),d_sub,'b-')
% ylim([0,ax.YLim(2)]);
title('PSTH')

ax = nexttile;
hold on;
plot(obj.time(edges),d_sub_trial,'b.','MarkerSize',10)
plot(obj.time(edges),d_sub_trial,'b-')
% ylim([0,ax.YLim(2)]);
title('single-trial')
















