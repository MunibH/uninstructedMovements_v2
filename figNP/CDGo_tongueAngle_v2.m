%% single trial heatmaps of cd go in n/p and tongue angle
close all

cols = getColors;

feat = 'tongue_angle';
featix = find(ismember(kin(1).featLeg,feat));
cond2use = [2 3 4 5];
for isess = 1:numel(meta)
    W = cd_null(isess).cd_mode_orth(:,2);
    trialdat.null{isess} = tensorprod(obj(isess).trialdat,W,2,1);
    W = cd_potent(isess).cd_mode_orth(:,2);
    trialdat.potent{isess} = tensorprod(obj(isess).trialdat,W,2,1);
    trialdat.angle{isess} =  kin(isess).dat(:,:,featix);

    trix = params(isess).trialid{cond2use};

    f = figure;
    f.Position = [1           1        1496         854];
    t = tiledlayout('flow');

    ax = nexttile;
    imagesc(obj(1).time,1:numel(trix),trialdat.angle{isess}(:,trix)')
    colormap(linspecer)
    title('tongue angle')

    ax = nexttile;
    imagesc(obj(1).time,1:numel(trix),trialdat.potent{isess}(:,trix)')
    colormap(linspecer)
    title('potent')


    ax = nexttile;
    imagesc(obj(1).time,1:numel(trix),trialdat.null{isess}(:,trix)')
    colormap(linspecer)
    title('null')

    c{1} = cols.rhit;
    c{2} = cols.lhit;
    c{3} = cols.rmiss;
    c{4} = cols.lmiss;
    alph = 0.2;

    ax = nexttile;
    hold on;
    for i = 1:numel(cond2use)
        trix = params(isess).trialid{cond2use(i)};
        tempangle = trialdat.angle{isess}(:,trix);
        mu.angle = mean(tempangle,2,'omitmissing');
        sd.angle = getCI(tempangle,1);
        plot(obj(isess).time,mu.angle,'Color',c{i},'LineWidth',2);
        % shadedErrorBar(obj(isess).time,mu.angle,sd.angle,{'Color',c{i},'LineWidth',2},alph,ax)
    end
    title('tongue angle')


    ax = nexttile;
    hold on;
    for i = 1:numel(cond2use)
        trix = params(isess).trialid{cond2use(i)};
        temppotent = trialdat.potent{isess}(:,trix);
        mu.potent = mean(temppotent,2);
        sd.potent = getCI(temppotent,1);
        shadedErrorBar(obj(isess).time,mu.potent,sd.potent,{'Color',c{i},'LineWidth',2},alph,ax)
    end
    title('Potent')


    ax = nexttile;
    hold on;
    for i = 1:numel(cond2use)
        trix = params(isess).trialid{cond2use(i)};
        tempnull = trialdat.null{isess}(:,trix);
        mu.null = mean(tempnull,2);
        sd.null = getCI(tempnull,1);
        shadedErrorBar(obj(isess).time,mu.null,sd.null,{'Color',c{i},'LineWidth',2},alph,ax)
    end
    title('Null')

    
   

    % break

    xlabel(t,'Time from go cue (s)')
    ylabel(t,'Trials')


end



%%