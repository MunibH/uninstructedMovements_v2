

%%
sessix = 3;

cond2use_ta = [1 2]; % right and left hits, corresponding to trial-avg projs onto n/p
cond2use_st = [2 3]; % right and left hits, corresponding to single-trial projs onto n/p

modparams.subTrials = 35;
edges = [mode(obj(1).bp.ev.delay) mode(obj(1).bp.ev.goCue)] - 2.5;

[plotsel.null, plotsel.null_sep, plotsel.null_sel_dims] = getSortedSelectivity(obj(sessix),params(sessix),rez(sessix),edges,cond2use_ta,cond2use_st,modparams,'null');
[plotsel.potent, plotsel.potent_sep, plotsel.potent_sel_dims] = getSortedSelectivity(obj(sessix),params(sessix),rez(sessix),edges,cond2use_ta,cond2use_st,modparams,'potent');

%%

xlims = [-2.3 2.5];

f = figure; f.Position = [794   516   302   384];
ax = gca;
hold on;
imagesc(obj(1).time, 1:size(plotsel.null,2),plotsel.null');
c = colorbar;
colormap(flipud(redblue))
xlabel('Time (s) from go cue')
ylabel('Dimensions')
c.Label.String = 'Selectivity (a.u.)';
title('null')
ax.YDir = 'reverse';
xlim(xlims)
ylim([1 size(plotsel.null,2)])
yline(plotsel.null_sep(1),'k--','LineWidth',2)
yline(plotsel.null_sep(2),'k--','LineWidth',2)


f = figure; f.Position = [794   516   302   384];
ax = gca;
hold on;
imagesc(obj(1).time, 1:size(plotsel.potent,2),plotsel.potent');
c = colorbar;
colormap(flipud(redblue))
xlabel('Time (s) from go cue')
ylabel('Dimensions')
c.Label.String = 'Selectivity (a.u.)';
title('potent')
ax.YDir = 'reverse';
xlim(xlims)
ylim([1 size(plotsel.potent,2)])
yline(plotsel.potent_sep(1),'k--','LineWidth',2)
yline(plotsel.potent_sep(2),'k--','LineWidth',2)



%%

scores = rez(sessix).N_null_psth;
Q = rez(sessix).Qnull;

dims = logical(plotsel.null_sel_dims);
recon.null = tensorprod(scores(:,dims,:),Q(:,dims),2);
recon.null = permute(recon.null,[1 3 2]); % (time, neurons, cond)


scores = rez(sessix).N_potent_psth;
Q = rez(sessix).Qpotent;

dims = logical(plotsel.potent_sel_dims);
recon.potent = tensorprod(scores(:,dims,:),Q(:,dims),2);
recon.potent = permute(recon.potent,[1 3 2]); % (time, neurons, cond)

%%
figure; imagesc(cat(2,recon.null(:,:,1),recon.null(:,:,2))')
title('null')

figure; imagesc(cat(2,recon.potent(:,:,1),recon.potent(:,:,2))')
title('potent')

figure; imagesc(cat(2,obj(sessix).psth(:,:,2),obj(sessix).psth(:,:,3))')
title('original')

%% variance explained by reconstructing data

for sessix = 1:numel(rez)

    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat,obj(sessix)); % (time, trials, neurons)

    scores = rez(sessix).N_null;
    Q = rez(sessix).Qnull;
    recon.null = tensorprod(scores,Q,3,2); % (time, trials, neurons)

    scores = rez(sessix).N_potent;
    Q = rez(sessix).Qpotent;

    recon.potent = tensorprod(scores,Q,3,2); % (time, trials, neurons)

    for cluix = 1:size(trialdat_zscored,3)
        dat = trialdat_zscored(:,:,cluix);
        pred = recon.null(:,:,cluix);
        ve.null{sessix}(cluix) = getR2(dat(:),pred(:));
    end
    for cluix = 1:size(trialdat_zscored,3)
        dat = trialdat_zscored(:,:,cluix);
        pred = recon.potent(:,:,cluix);
        ve.potent{sessix}(cluix) = getR2(dat(:),pred(:));
    end
end









