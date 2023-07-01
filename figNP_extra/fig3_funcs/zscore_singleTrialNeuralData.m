function trialdat_zscored = zscore_singleTrialNeuralData(obj)
% zscore obj.trialdat, which is single trial binned neural data
% (time,trials,neurons)
% The zscoring takes place on the matrix of (time*trials,neurons) for all
% neurons
% trialdat_zscored is of size(time,trials,neurons)

% trialdat_zscored = permute(dat,[1 3 2]); % (time,trials,neurons)

t1 = mode(obj.bp.ev.bitStart - 2.5);
t2 = mode(obj.bp.ev.sample - 2.5);
[~,ix1] = min(abs(obj.time - t1));
[~,ix2] = min(abs(obj.time - t2));

dat = obj.trialdat; % (time,neurons,trials)

temp = permute(dat,[1 3 2]);
dims = size(temp); % (time,trials,neurons)

% method 1 (pressample)
% psth_ = squeeze(mean(dat,3,'omitmissing'));
% presampleMean = mean(psth_(ix1:ix2,:),1);
% presampleStd = std(psth_(ix1:ix2,:),[],1);
% % presampleMean = mean(temp(ix1:ix2,:,:),1,'omitmissing');
% % presampleStd = std(temp(ix1:ix2,:,:), [], 1,'omitmissing');
% % presampleMean = obj.presampleFR(:,1)';
% % presampleStd = obj.presampleSigma(:,1)';
% % presampleMean = reshape(presampleMean,1,1,size(presampleMean,2));
% % presampleStd = reshape(presampleStd,1,1,size(presampleStd,2));
% % mu = repmat(presampleMean,size(temp,1),1,size(temp,2));
% % sd = repmat(presampleStd,size(temp,1),1,size(temp,2));
% trialdat_zscored = (dat - presampleMean) ./ presampleStd;
% trialdat_zscored = permute(trialdat_zscored, [1 3 2]);
% for i = 1:size(trialdat_zscored,3)
%     trialdat_zscored(:,:,i) = fillmissing(trialdat_zscored(:,:,i),'constant',presampleMean(i));
% end

% trialdat_zscored = (temp - permute(mu,[1 3 2]) ) ./ permute(sd,[1 3 2]);

% % method 2, zscore across all time*trials
temp = reshape(temp,dims(1)*dims(2),dims(3));
temp2 = zscore(temp);
trialdat_zscored = reshape(temp2,dims(1),dims(2),dims(3));

% % method 3, subtract appropriate PSTH (Left or Right trial)
% trialdat_zscored = dat;
% % for cluix = 1:size(dat,2) % loop through clusters
% for trix = 1:size(dat,3) % loop through trials
%     if obj.bp.R(trix)
%         cond = psth_cond(1);
%     elseif obj.bp.L(trix)
%         cond = psth_cond(2);
%     end
% 
%     mu = obj.psth(:,:,cond);
%     trialdat_zscored(:,:,trix) = trialdat_zscored(:,:,trix) - mu; 
%     sigma = std(dat(:,:,params.trialid{cond}),[],3);
%     % smooth sigma to avoid dividing by zero
%     % division by zeros happen b/c for low firing neurons and a small bin
%     % size, there are zero variance time points across trials
%     sigma = mySmooth(sigma,21,'reflect');
%     trialdat_zscored(:,:,trix) = trialdat_zscored(:,:,trix) ./ sigma; 
%     if any(isnan(trialdat_zscored))
%         'a'
%     end
% end
% % end



end







