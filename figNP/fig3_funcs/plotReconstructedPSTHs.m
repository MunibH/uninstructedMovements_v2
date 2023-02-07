function plotReconstructedPSTHs(obj,rez,params)

cols = getColors;

f = figure;
for sessix = 1:numel(obj)
%     dat = nanmean(obj(sessix).trialdat,3);
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix));
    dat = squeeze(nanmean(trialdat_zscored,2));
    null = squeeze(nanmean(rez(sessix).recon.null,2));
    potent = squeeze(nanmean(rez(sessix).recon.potent,2));
    
    for cluix = 1:size(dat,2)
        ax = nexttile;
        hold on;
        plot(obj(sessix).time,dat(:,cluix),'k')
        plot(obj(sessix).time,null(:,cluix),'Color',cols.null)
        plot(obj(sessix).time,potent(:,cluix),'Color',cols.potent)
%         title(['NullVE=' num2str(rez(sessix).ve_recon.null(cluix))  '| PotentVE=' num2str(rez(sessix).ve_recon.potent(cluix))])
        if cluix > 40
            break
        end
    end

    pause
    clf
end


end