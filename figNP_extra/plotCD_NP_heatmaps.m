%% choice
close all


cdix = 1;
cond2use = [8 9]; 
for sessix = 20%[5 15 16 17 20 21]%1:numel(meta) % [5 15 16 17 20 21]
    trix = cell2mat(params(sessix).trialid(cond2use)');

    dat = rez(sessix).recon.null;
    W = cd_null(sessix).cd_mode_orth(:,cdix);
    proj.null = tensorprod(dat,W,3,1);
    proj.null = zscore(proj.null(:,trix));
    proj.null = mySmooth(proj.null,41,'zeropad');

    dat = rez(sessix).recon.potent;
    W = cd_potent(sessix).cd_mode_orth(:,cdix);
    proj.potent = tensorprod(dat,W,3,1);
    proj.potent = zscore(proj.potent(:,trix));
    proj.potent = mySmooth(proj.potent,21,'zeropad');

    % baseline subtract
    bix = 1:10;
    mu = mean(proj.null(bix,:),1);
    proj.null = proj.null - mu;
    mu = mean(proj.potent(bix,:),1);
    proj.potent = proj.potent - mu;
    
    f = figure;
    f.Renderer = 'painters';
    
    ax = subplot(1,2,1);
    imagesc(obj(1).time,1:size(proj.null,2),proj.null'); colorbar; clim([-3 3])
    colormap(redblue)
    xl = numel(params(sessix).trialid{cond2use(1)});
    line(ax.XLim, [xl xl],'linestyle','--')
    ax = prettifyPlot(ax);

    ax = subplot(1,2,2);
    imagesc(obj(1).time,1:size(proj.null,2),proj.potent'); colorbar; clim([-4 4])
    colormap(redblue)
    line(ax.XLim, [xl xl],'linestyle','--')
    ax = prettifyPlot(ax);

    sgtitle([meta(sessix).anm ' ' meta(sessix).date])

    % break

end


%% ramping

close all

cdix = 3;
cond2use = [8 9];
for sessix = 10:numel(meta)
    trix = cell2mat(params(sessix).trialid(cond2use)');
    trix = sort(trix,'ascend');

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
    imagesc(obj(1).time,1:size(proj.null,2),mySmooth(proj.null,11)'); colorbar; clim([-3 3])
    colormap(parula)
    ax = prettifyPlot(ax);

    ax = subplot(1,2,2);
    imagesc(obj(1).time,1:size(proj.null,2),proj.potent'); colorbar; %clim([-5 5])
    colormap(parula)
    ax = prettifyPlot(ax);

    sgtitle([meta(sessix).anm ' ' meta(sessix).date])

    break

end
