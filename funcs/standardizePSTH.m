function psth = standardizePSTH(obj)
% standardize by subtracting baseline firing rate, dividing by baseline std dev

psth = zeros(size(obj.psth));
for i = 1:size(obj.psth,3)
    
    temp = obj.psth(:,:,i);

    psth(:,:,i) = temp - mean(temp,2);
    
    % standardize using presample stats
    % mu = mean(temp(1:20,:),1);
    % sd = std(temp(1:20,:),[],1);
    % psth(:,:,i) = (temp - mu) ./ sd;
    
end

end