function trialdat_zscored = zscore_singleTrialNeuralData(dat,obj)
% zscore obj.trialdat, which is single trial binned neural data
% (time,neurons,trials)
% The zscoring takes place on the matrix of (time*trials,neurons) for all
% neurons
% trialdat_zscored is of size(time,trials,neurons)



t1 = mode(obj.bp.ev.bitStart - 2.5);
t2 = mode(obj.bp.ev.sample - 2.5);
[~,ix1] = min(abs(obj.time - t1));
[~,ix2] = min(abs(obj.time - t2));


temp = permute(dat,[1 3 2]);
dims = size(temp); % (time,trials,neurons)

% method 1 (pressample)
% presampleMean = mean(temp(ix1:ix2,:,:),1);
% presampleStd = std(temp(ix1:ix2,:,:), [], 1);
% presampleMean = obj.presampleFR(:,1)';
% presampleStd = obj.presampleSigma(:,1)';
% presampleMean = reshape(presampleMean,1,1,size(presampleMean,2));
% presampleStd = reshape(presampleStd,1,1,size(presampleStd,2));
% trialdat_zscored = (temp - presampleMean) ./ presampleStd;

% % method 2, zscore across all time*trials
temp = reshape(temp,dims(1)*dims(2),dims(3));
temp2 = zscore(temp);
temp2 = temp ./ max(temp);
trialdat_zscored = reshape(temp2,dims(1),dims(2),dims(3));




end