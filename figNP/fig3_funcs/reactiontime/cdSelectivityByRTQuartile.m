function cdSelectivityByRTQuartile(obj,cd_null,cd_potent,params,nQuantiles, rt)

%%

cd_labels = cd_null(1).cd_labels;
nCDs = numel(cd_labels);

tm_mask = cd_null(1).cd_times;
fns = fieldnames(tm_mask);
for i = 1:numel(fns)
    ix{i} = tm_mask.(fns{i});
end

cond2use = [8 9]; % right, left hits , excluding early

for sessix = 1:numel(obj)
    trix.all = cell2mat(params(sessix).trialid(cond2use)');
    trix.right = ismember(trix.all,params(sessix).trialid{cond2use(1)});
    trix.left = ismember(trix.all,params(sessix).trialid{cond2use(2)});

    dat.rt = rt{sessix}(trix.all);
%     dat.me = me(sessix).data(:,trix.all);

    quantiles = quantile(dat.rt,nQuantiles-1);
    trialsBlock = getTrialsByRTBlock(dat.rt,quantiles);
    %     cellfun(@(x) sum(x), trialsBlock, 'uniformoutput',false)

    for icd = 1:nCDs
        dat.null{1} = cd_null(sessix).trialdat(:,trix.right,icd);
        dat.potent{1}= cd_potent(sessix).trialdat(:,trix.right,icd);
        dat.null{2} = cd_null(sessix).trialdat(:,trix.left,icd);
        dat.potent{2} = cd_potent(sessix).trialdat(:,trix.left,icd);

        % null preferred
        [~,pref] = max( [mean(mean(dat.null{1}(ix{icd},:))) , mean(mean(dat.null{2}(ix{icd},:)))] );
        if pref==1; nonpref=2; else; nonpref=1; end
        dat.null_pref = pref; dat.null_nonpref = nonpref;
        % potent preferred
        [~,pref] = max( [mean(mean(dat.potent{1}(ix{icd},:))) , mean(mean(dat.potent{2}(ix{icd},:)))] );
        if pref==1; nonpref=2; else; nonpref=1; end
        dat.potent_pref = pref; dat.potent_nonpref = nonpref;

        for iq = 1:nQuantiles
            iqTrials{1}  = trialsBlock{iq}' & trix.right;
            iqTrials{2}  = trialsBlock{iq}' & trix.left;

            dat.sel.null(:,iq,icd,sessix) = mean(cd_null(sessix).trialdat(:,iqTrials{1},icd),2) ... % (time,quartile,cd,sessix)
                - mean(cd_null(sessix).trialdat(:,iqTrials{2},icd),2);

            dat.sel.potent(:,iq,icd,sessix) = mean(cd_potent(sessix).trialdat(:,iqTrials{1},icd),2) ...
                - mean(cd_potent(sessix).trialdat(:,iqTrials{2},icd),2);
        end
    end


end

%%
close all

sample = mode(obj(1).bp.ev.sample) - 2.5;
tstart = mode(obj(1).bp.ev.bitStart) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;

sm = 21;

xlims = [tstart 1.5];

c = gray(nQuantiles*4);
c = linspace(0.15,0.8,nQuantiles);
for i = 1:nQuantiles
    col(i,:) = [c(i) c(i) c(i)];
end

alph = 0.05;

% null space
f = figure;
t = tiledlayout('flow');
for icd = 1:nCDs
    ax = nexttile;
    hold on;
    for iq = 1:nQuantiles
        temp = squeeze(dat.sel.null(:,iq,icd,:));
        temp = mySmooth(temp,sm,'zeropad');
        mu = mean(temp,2);
        sig = std(temp,[],2) ./ numel(obj);
        %         shadedErrorBar(obj(1).time,mu,sig,{'Color',col(iq,:)},alph,ax)
        plot(obj(1).time,mu,'Color',col(iq,:))
    end
    title(cd_labels{icd})
    xlim(xlims)
    xline(sample,'k--')
    xline(delay,'k--')
    xline(0,'k--')
end
sgtitle('Null')


% potent space
f = figure;
t = tiledlayout('flow');
for icd = 1:nCDs
    ax = nexttile;
    hold on;
    for iq = 1:nQuantiles
        temp = squeeze(dat.sel.potent(:,iq,icd,:));
        temp = mySmooth(temp,sm,'zeropad');
        mu = mean(temp,2);
        sig = std(temp,[],2) ./ numel(obj);
        %         shadedErrorBar(obj(1).time,mu,sig,{'Color',col(iq,:)},alph,ax)
        plot(obj(1).time,mu,'Color',col(iq,:))
    end
    title(cd_labels{icd})
    xlim(xlims)
    xline(sample,'k--')
    xline(delay,'k--')
    xline(0,'k--')
end
sgtitle('Potent')

end