function [ve,r2] = varExpRampingNP(obj,rez,cd,params,rampix,type)


for sessix = 1:numel(obj)
    
%     ramp = cd
    orig = rez(sessix).recon_psth.(type);
    orig = cat(1, orig(:,:,1), orig(:,:,2));
    C = cov(orig);
    s = eig(C);

    rampmode = cd(sessix).cd_mode_orth(:,rampix);
    
    ve(sessix) = var_proj(rampmode,C,sum(s));


    scores = cd(sessix).trialdat(:,:,rampix);
    recon = tensorprod(scores,rampmode,3,2); % (time, trials, neurons)
    recon = squeeze(mean(recon,2));

    dat = zscore_singleTrialNeuralData(obj(sessix)); % (time, trials, neurons)
    dat = squeeze(mean(dat,2));

    r2(sessix) = corr(dat(:),recon(:)).^2;
    
end
