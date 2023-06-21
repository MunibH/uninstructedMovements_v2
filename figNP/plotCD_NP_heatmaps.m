%% choice
close all


cdix = 1;
cond2use = [8 9];
for sessix = 2:numel(meta)
    trix = cell2mat(params(sessix).trialid(cond2use)');

    dat = rez(sessix).recon.null;
    W = cd_null(sessix).cd_mode_orth(:,cdix);
    proj.null = tensorprod(dat,W,3,1);
    proj.null = zscore(proj.null(:,trix));

    dat = rez(sessix).recon.potent;
    W = cd_potent(sessix).cd_mode_orth(:,cdix);
    proj.potent = tensorprod(dat,W,3,1);
    proj.potent = zscore(proj.potent(:,trix));
    
    f = figure;
    f.Renderer = 'painters';
    
    ax = subplot(1,2,1);
    imagesc(proj.null'); colorbar; %clim([-5 10])
    ax = prettifyPlot(ax);

    ax = subplot(1,2,2);
    imagesc(proj.potent'); colorbar; %clim([-5 5])
    ax = prettifyPlot(ax);

    sgtitle([meta(sessix).anm ' ' meta(sessix).date])

    % break

end


%%


cdix = 3;

